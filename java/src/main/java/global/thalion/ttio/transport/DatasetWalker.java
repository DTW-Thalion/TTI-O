/*
 * TTI-O Java Implementation.
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.transport;

import global.thalion.ttio.AcquisitionRun;
import global.thalion.ttio.Enums;
import global.thalion.ttio.IRImage;
import global.thalion.ttio.MSImage;
import global.thalion.ttio.RamanImage;
import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.ReferenceImport;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Walks a {@link SpectralDataset}, applying an optional
 * {@link AUFilter} and dispatching transport-stream events to an
 * {@link AccessUnitVisitor}. Stateless once constructed; safe to
 * reuse across walks.
 *
 * <p>Iteration order matches
 * {@link TransportWriter#writeDataset(SpectralDataset)} so the
 * sequence of visitor calls is byte-equivalent to canonical dataset
 * emission when the visitor is one that encodes events as packets.
 *
 * <p>Cross-language equivalents:
 * <ul>
 *   <li>Objective-C: {@code TTIODatasetWalker}</li>
 *   <li>Python:       {@code ttio.transport.walker.walk_dataset}</li>
 * </ul>
 */
public final class DatasetWalker {

    /** Walks {@code dataset} event-by-event. Returns silently when
     *  the dataset or visitor is null. */
    public void walk(SpectralDataset dataset,
                      AUFilter filter,
                      AccessUnitVisitor visitor) {
        if (dataset == null || visitor == null) return;
        AUFilter flt = (filter == null) ? new AUFilter() : filter;

        Map<String, AcquisitionRun> runs = dataset.msRuns();
        Map<String, GenomicRun> genomicRuns = dataset.genomicRuns();
        List<String> features = new ArrayList<>();
        for (String f : dataset.featureFlags().features()) features.add(f);

        // 1. StreamHeader — n_datasets is spectral + genomic, matching
        // TransportWriter.writeDataset.
        visitor.visitStreamHeader(this, "1.2",
                                    dataset.title() == null
                                        ? "" : dataset.title(),
                                    dataset.isaInvestigationId() == null
                                        ? "" : dataset.isaInvestigationId(),
                                    features,
                                    runs.size() + genomicRuns.size());

        // 1a. v0.11 §5.4 prelude — match TransportWriter.writeDataset
        // ordering verbatim. Each gate uses the same "populated?" check
        // the writer uses, so the visitor sees events in the same order
        // (and only when present) as the on-wire packets the writer
        // emits. ObjC parity: TTIODatasetWalker.m v0.11 prelude block.
        // §5.4.1 ENCRYPTION_ALGORITHM
        if (dataset.isEncrypted()
            && dataset.encryptedAlgorithm() != null
            && !dataset.encryptedAlgorithm().isEmpty()) {
            visitor.visitEncryptionAlgorithm(this,
                dataset.encryptedAlgorithm());
        }
        // §5.4.2 DATASET_PROVENANCE
        if (dataset.provenanceRecords() != null
            && !dataset.provenanceRecords().isEmpty()) {
            visitor.visitDatasetProvenance(this, dataset.provenanceRecords());
        }
        // §5.4.3 SUBJECT_METADATA → SAMPLE_METADATA (subjects first so
        // a soft-FK target is visible ahead of any sample row that
        // references it).
        if (dataset.subjects() != null && !dataset.subjects().isEmpty()) {
            visitor.visitSubjectMetadata(this, dataset.subjects());
        }
        if (dataset.samples() != null && !dataset.samples().isEmpty()) {
            visitor.visitSampleMetadata(this, dataset.samples());
        }
        // §5.4.4 reference groups — sorted by URI key for determinism
        // (matches ObjC walker which sorts the NSDictionary keys).
        Map<String, ReferenceImport> refs = dataset.references();
        if (refs != null && !refs.isEmpty()) {
            for (Map.Entry<String, ReferenceImport> e
                 : new TreeMap<>(refs).entrySet()) {
                visitor.visitReferenceGroup(this, e.getValue());
            }
        }
        // §5.4.5 image cubes — MS → Raman → IR.
        MSImage ms = (MSImage) dataset.imageForKind(Enums.ImageKind.MS);
        if (ms != null) {
            visitor.visitImage(this, ms);
        }
        RamanImage raman = (RamanImage) dataset.imageForKind(Enums.ImageKind.RAMAN);
        if (raman != null) {
            visitor.visitRamanImage(this, raman);
        }
        IRImage ir = (IRImage) dataset.imageForKind(Enums.ImageKind.IR);
        if (ir != null) {
            visitor.visitIRImage(this, ir);
        }
        // §5.4.6 IDENTIFICATIONS_TABLE → QUANTIFICATIONS_TABLE
        if (dataset.identifications() != null
            && !dataset.identifications().isEmpty()) {
            visitor.visitIdentificationsTable(this, dataset.identifications());
        }
        if (dataset.quantifications() != null
            && !dataset.quantifications().isEmpty()) {
            visitor.visitQuantificationsTable(this,
                dataset.quantifications());
        }

        // 2. DatasetHeaders — spectral runs (ids 1..N), then genomic
        // runs (ids N+1..N+M). Matches TransportWriter.writeDataset.
        int did = 1;
        for (Map.Entry<String, AcquisitionRun> e : runs.entrySet()) {
            if (flt.datasetId != null && did != flt.datasetId) {
                did++;
                continue;
            }
            AcquisitionRun run = e.getValue();
            List<String> channelNames =
                new ArrayList<>(run.channels().keySet());
            String instrumentJson =
                TransportWriter.instrumentConfigJson(run.instrumentConfig());
            visitor.visitDatasetHeader(this, did, e.getKey(),
                                         run.acquisitionMode().ordinal(),
                                         run.spectrumClassName(),
                                         channelNames,
                                         instrumentJson,
                                         run.spectrumCount());
            did++;
        }
        for (Map.Entry<String, GenomicRun> e : genomicRuns.entrySet()) {
            if (flt.datasetId != null && did != flt.datasetId) {
                did++;
                continue;
            }
            GenomicRun grun = e.getValue();
            visitor.visitDatasetHeader(this, did, e.getKey(),
                                         grun.acquisitionMode().ordinal(),
                                         "TTIOGenomicRead",
                                         TransportWriter.genomicChannelNames(grun),
                                         TransportWriter.genomicRunMetadataJson(grun),
                                         grun.readCount());
            did++;
        }

        // 3. AccessUnits + 4. EndOfDataset, interleaved per dataset to
        // match TransportWriter.writeDataset: AUs(ds_k) → EOD(ds_k).
        int emitted = 0;
        Integer maxAu = flt.maxAu;
        did = 1;
        outer:
        for (Map.Entry<String, AcquisitionRun> e : runs.entrySet()) {
            if (flt.datasetId != null && did != flt.datasetId) {
                did++;
                continue;
            }
            AcquisitionRun run = e.getValue();
            int count = run.spectrumCount();
            List<String> channelNames =
                new ArrayList<>(run.channels().keySet());
            for (int i = 0; i < count; i++) {
                if (maxAu != null && emitted >= maxAu) {
                    visitor.visitEndOfDataset(this, did, count);
                    break outer;
                }
                AccessUnit au = TransportWriter.spectrumToAccessUnit(
                    run, i, channelNames);
                if (!flt.matches(au, did)) continue;
                visitor.visitAccessUnit(this, au, did, i);
                emitted++;
            }
            visitor.visitEndOfDataset(this, did, count);
            did++;
        }
        outerGenomic:
        for (Map.Entry<String, GenomicRun> e : genomicRuns.entrySet()) {
            if (flt.datasetId != null && did != flt.datasetId) {
                did++;
                continue;
            }
            GenomicRun grun = e.getValue();
            List<AccessUnit> aus =
                TransportWriter.genomicRunAccessUnits(grun);
            for (int i = 0; i < aus.size(); i++) {
                if (maxAu != null && emitted >= maxAu) {
                    visitor.visitEndOfDataset(this, did, grun.readCount());
                    break outerGenomic;
                }
                AccessUnit au = aus.get(i);
                if (!flt.matches(au, did)) continue;
                visitor.visitAccessUnit(this, au, did, i);
                emitted++;
            }
            visitor.visitEndOfDataset(this, did, grun.readCount());
            did++;
        }

        // 5. EndOfStream
        visitor.visitEndOfStream(this);
    }
}
