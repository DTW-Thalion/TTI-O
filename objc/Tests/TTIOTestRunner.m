#import <Foundation/Foundation.h>
#import "Testing.h"

extern void testValueRange(void);
extern void testIsolationWindow(void);
extern void testActivationMethodEnum(void);
extern void testSpectrumKind(void);
extern void testMassSpectrumActivationAndIsolationFields(void);
extern void testSpectrumIndexM74RoundTrip(void);
extern void testSpectralDatasetM74FeatureFlag(void);
extern void testEncodingSpec(void);
extern void testAxisDescriptor(void);
extern void testCVParam(void);
extern void testHDF5(void);
extern void testSignalArray(void);
extern void testSpectra(void);
extern void testAcquisitionRun(void);
extern void testSpectralDataset(void);
extern void testMSImage(void);
extern void testEncryption(void);
extern void testQueryAndStreaming(void);
extern void testMzMLReader(void);
extern void testMzMLReaderM74(void);
extern void testMzMLReaderHugeBinaryNode(void);
extern void testMilestone10(void);
extern void testMilestone11(void);
extern void testMilestone12(void);
extern void testNmrMLReaderHugeFidNode(void);
extern void testMilestone13(void);
extern void testMilestone14(void);
extern void testMilestone17(void);
extern void testMilestone18(void);
extern void testMilestone19(void);
extern void testMilestone19M74(void);
extern void testMilestone21(void);
extern void testMilestone23(void);
extern void testMilestone24(void);
extern void testMilestone25(void);
extern void testMilestone27(void);
extern void testMilestone28(void);
extern void testMilestone29(void);
extern void testNmrMLReaderParity(void);
extern void testMilestone39(void);
extern void testMilestone49(void);
extern void testMilestone52(void);
extern void testMilestone53(void);
extern void testImzMLReader(void);
extern void testMzTabReader(void);
extern void testStress(void);
extern void testWatersMassLynxReader(void);
extern void testWriteMinimal(void);
extern void testImzMLWriter(void);
extern void testMzTabWriter(void);
extern void testSqliteProvider(void);
extern void testCloudAccess(void);
extern void testCanonicalBytesCrossBackend(void);
extern void testCipherSuite(void);
extern void testNdDatasetCrossBackend(void);
extern void testTransportCodec(void);
extern void testTransportWriterSink(void);
extern void testDatasetWalker(void);
extern void testDatasetWalkerV011Prelude(void);
extern void testTransportIngest(void);
extern void testTransportClient(void);
extern void testAcquisitionSimulator(void);
extern void testTransportServer(void);
extern void testTransportConformance(void);
extern void testSelectiveAccess(void);
extern void testPerAUEncryption(void);
extern void testPerAUFile(void);
extern void testPerAUFileDecryptInPlace(void);
extern void testEncryptedTransport(void);
extern void testMilestone73(void);
extern void testMilestone73_1(void);
extern void testM76JcampConformance(void);
extern void testMilestone77(void);
extern void testMilestone78(void);
extern void testV11EncryptionParity(void);
extern void testV111DecryptInPlace(void);
extern void testM79GenomicEnums(void);
extern void testM82GenomicRun(void);
extern void testM83Rans(void);
extern void testCodecRegistry(void);
extern void testImportedDataset(void);
extern void testRunSelection(void);
extern void testBrukerReaderDataset(void);
extern void testFormatReaders(void);
extern void testFormatWriters(void);
extern void testImporterExporterRegistry(void);
extern void testRegistryParity(void);
extern void testEncodeExportCli(void);
extern void testImageBase(void);
extern void testM84BasePack(void);
extern void testPackedReference(void);
extern void testFloatDeltaZstdPlainTransform(void);
extern void testFloatDeltaZstd(void);
extern void testQualitiesV5(void);
extern void testM85Quality(void);
extern void testM86GenomicCodecWiring(void);
extern void testM94ZFqzcompPerf(void);
extern void testTtioRansNative(void);
extern void testM94ZV4Dispatch(void);
extern void testM94ZV4ByteExact(void);
extern void testMateInfoV2(void);
extern void testMateInfoV2Dispatch(void);
extern void testRefDiffV2(void);
extern void testRefDiffV2Dispatch(void);
extern void testM97LongReadProfile(void);
extern void testM98AssemblyGraph(void);
extern void testM99BlocksV1PerAU(void);
extern void testM101SamTags(void);
extern void testNameTokenizerV2(void);
extern void testNameTokenizedV2Dispatch(void);
extern void testOffsetsCumsum(void);
extern void testGenomicIndexRegionQuery(void);
extern void testM95DeltaRansUnit(void);
extern void testM87BamImporter(void);
extern void testM88CramBamRoundTrip(void);
extern void testM89GenomicTransport(void);
extern void testM90GenomicProtection(void);
extern void testM90Parity(void);
extern void testM90Final(void);
extern void testV4EdgeCases(void);
extern void testV8Hdf5Corruption(void);
extern void testC1ToolsCli(void);
extern void testC2HDF5ErrorPaths(void);
extern void testC3ProvidersErrorPaths(void);
extern void testC5ProtectionGap(void);
extern void testC3bProvidersWritePaths(void);
extern void testProviderExtendable(void);
extern void testGenomicBlocks(void);
extern void testGenomicStreamWriter(void);
extern void testGenomicStreamWriterThreads(void);
extern void testSpectralStreamWriterThreads(void);
extern void testSpectralIterBlocks(void);
extern void testFastqRecordScanner(void);
extern void testThreads(void);
extern void testBlocksV1Golden(void);
extern void testGenomicBlocksReader(void);
extern void testReferenceResolverExternal(void);
extern void testFloatDeltaZstdBlocks(void);
extern void testSpectralStreamWriter(void);
extern void testStreamingImporters(void);
extern void testStreamingExporters(void);
extern void testC2bHDF5CompoundType(void);
extern void testC3cCanonicalBytes(void);
extern void testPhase12RunProtocol(void);
extern void testTask30MSProviderURL(void);
extern void testTask31InstanceWriterParity(void);
extern void testFastaFastqIo(void);
extern void testFastaFastqTioRoundTrip(void);
extern void testProgressSink(void);
extern void testProgressSinkStageC(void);
extern void testProgressSinkStageD(void);
extern void testReferencesAccessor(void);
extern void testReferenceImportWriteToDataset(void);
extern void testSpectralDatasetWritableOpen(void);
extern void testMSImageMzAxis(void);
extern void testMSImageLegacyMzAxisAbsent(void);
extern void testSpectralDatasetMsImageAccessor(void);
extern void testPixelSpectraRaisesWhenMzAxisAbsent(void);
extern void testJcampVibrationalRoundTrip(void);
extern void testFD1MultiRecipient(void);
extern void testMultiRecipientXLang(void);
extern void testDekWrappedXLang(void);
extern void testC0StandaloneKeyWrap(void);
extern void testC2aServerKekId(void);
extern void testTransportPacketTypeV011(void);
extern void testTransportReaderSkipUnknown(void);
extern void testTransportWriterReference(void);
extern void testTransportReaderReference(void);
extern void testArrowIpcCodec(void);
extern void testTransportEncryptionAlgorithm(void);
extern void testTransportDatasetProvenance(void);
extern void testTransportImage(void);
extern void testTransportImageProcessed(void);
extern void testSpectralDatasetIRImage(void);
extern void testSpectralDatasetImages(void);
extern void testIRImageWriteReadTypedAttrs(void);
extern void testTransportRamanImage(void);
extern void testTransportIRImage(void);
extern void testTransportIdentificationsQuantifications(void);
extern void testTransportWriteDatasetPrelude(void);
extern void testSubjectSample(void);
extern void testSpectralDatasetSubjectsSamples(void);
extern void testArrowIpcCodecSubjectsSamples(void);
extern void testTransportSubjectsSamples(void);
extern void testAccessorMatrixConformance(void);
extern void testCoverageGapWatchdog(void);

int main(int argc, const char *argv[])
{
    @autoreleasepool {
        START_SET("TTIOValueRange")
            testValueRange();
        END_SET("TTIOValueRange")

        START_SET("TTIOIsolationWindow")
            testIsolationWindow();
        END_SET("TTIOIsolationWindow")

        START_SET("TTIOActivationMethod")
            testActivationMethodEnum();
        END_SET("TTIOActivationMethod")

        START_SET("TTIOSpectrumKind")
            testSpectrumKind();
        END_SET("TTIOSpectrumKind")

        START_SET("TTIOMassSpectrum M74 fields")
            testMassSpectrumActivationAndIsolationFields();
        END_SET("TTIOMassSpectrum M74 fields")

        START_SET("TTIOSpectrumIndex M74 round-trip")
            testSpectrumIndexM74RoundTrip();
        END_SET("TTIOSpectrumIndex M74 round-trip")

        START_SET("M74 Slice E: opt_ms2_activation_detail + format 1.3")
            testSpectralDatasetM74FeatureFlag();
        END_SET("M74 Slice E: opt_ms2_activation_detail + format 1.3")

        START_SET("TTIOEncodingSpec")
            testEncodingSpec();
        END_SET("TTIOEncodingSpec")

        START_SET("TTIOAxisDescriptor")
            testAxisDescriptor();
        END_SET("TTIOAxisDescriptor")

        START_SET("TTIOCVParam")
            testCVParam();
        END_SET("TTIOCVParam")

        START_SET("TTIOHDF5 wrappers")
            testHDF5();
        END_SET("TTIOHDF5 wrappers")

        START_SET("TTIOSignalArray")
            testSignalArray();
        END_SET("TTIOSignalArray")

        START_SET("Spectrum classes")
            testSpectra();
        END_SET("Spectrum classes")

        START_SET("AcquisitionRun + SpectrumIndex")
            testAcquisitionRun();
        END_SET("AcquisitionRun + SpectrumIndex")

        START_SET("SpectralDataset (multi-run)")
            testSpectralDataset();
        END_SET("SpectralDataset (multi-run)")

        START_SET("MSImage")
            testMSImage();
        END_SET("MSImage")

        START_SET("Encryption (AES-256-GCM)")
            testEncryption();
        END_SET("Encryption (AES-256-GCM)")

        START_SET("Query + Streaming")
            testQueryAndStreaming();
        END_SET("Query + Streaming")

        START_SET("mzML Reader (Milestone 9)")
            testMzMLReader();
        END_SET("mzML Reader (Milestone 9)")

        START_SET("mzML Reader M74 (activation + isolation)")
            testMzMLReaderM74();
        END_SET("mzML Reader M74 (activation + isolation)")

        START_SET("mzML Reader: 10 MB+ binary node, streamed from a path")
            testMzMLReaderHugeBinaryNode();
        END_SET("mzML Reader: 10 MB+ binary node, streamed from a path")

        START_SET("Milestone 10: protocols + modality-agnostic")
            testMilestone10();
        END_SET("Milestone 10: protocols + modality-agnostic")

        START_SET("Milestone 11: compound types + dataset encryption")
            testMilestone11();
        END_SET("Milestone 11: compound types + dataset encryption")

        START_SET("Milestone 12: MSImage inheritance + native 2D NMR")
            testMilestone12();
        END_SET("Milestone 12: MSImage inheritance + native 2D NMR")

        START_SET("Milestone 13: nmrML reader")
            testMilestone13();
        END_SET("Milestone 13: nmrML reader")

        START_SET("nmrML reader: 10 MB+ FID node, streamed from a path")
            testNmrMLReaderHugeFidNode();
        END_SET("nmrML reader: 10 MB+ FID node, streamed from a path")

        START_SET("Milestone 14: digital signatures + verification")
            testMilestone14();
        END_SET("Milestone 14: digital signatures + verification")

        START_SET("Milestone 17: compound per-run provenance")
            testMilestone17();
        END_SET("Milestone 17: compound per-run provenance")

        START_SET("Milestone 18: canonical byte-order signatures")
            testMilestone18();
        END_SET("Milestone 18: canonical byte-order signatures")

        START_SET("Milestone 19: mzML writer + indexedmzML")
            testMilestone19();
        END_SET("Milestone 19: mzML writer + indexedmzML")

        START_SET("Milestone 19 M74: writer activation + isolation")
            testMilestone19M74();
        END_SET("Milestone 19 M74: writer activation + isolation")

        START_SET("Milestone 21: LZ4 + Numpress-delta codecs")
            testMilestone21();
        END_SET("Milestone 21: LZ4 + Numpress-delta codecs")

        START_SET("Milestone 23: thread-safe TTIOHDF5File")
            testMilestone23();
        END_SET("Milestone 23: thread-safe TTIOHDF5File")

        START_SET("Milestone 24: chromatogram API + mzML writer completion")
            testMilestone24();
        END_SET("Milestone 24: chromatogram API + mzML writer completion")

        START_SET("Milestone 25: envelope encryption + key rotation")
            testMilestone25();
        END_SET("Milestone 25: envelope encryption + key rotation")

        START_SET("Milestone 27: ISA-Tab / ISA-JSON exporter")
            testMilestone27();
        END_SET("Milestone 27: ISA-Tab / ISA-JSON exporter")

        START_SET("Milestone 28: spectral anonymization")
            testMilestone28();
        END_SET("Milestone 28: spectral anonymization")

        START_SET("Milestone 29: nmrML writer + Thermo RAW stub")
            testMilestone29();
            testNmrMLReaderParity();
        END_SET("Milestone 29: nmrML writer + Thermo RAW stub")

        START_SET("Milestone 39: provider abstraction")
            testMilestone39();
        END_SET("Milestone 39: provider abstraction")

        START_SET("Milestone 41: SQLite storage provider")
            testSqliteProvider();
        END_SET("Milestone 41: SQLite storage provider")

        START_SET("Cloud access (ROS3 / S3)")
            testCloudAccess();
        END_SET("Cloud access (ROS3 / S3)")

        START_SET("M43: canonical bytes cross-backend")
            testCanonicalBytesCrossBackend();
        END_SET("M43: canonical bytes cross-backend")

        START_SET("M48: cipher suite catalog")
            testCipherSuite();
        END_SET("M48: cipher suite catalog")

        START_SET("M45: N-D dataset cross-backend")
            testNdDatasetCrossBackend();
        END_SET("M45: N-D dataset cross-backend")

        START_SET("M49: post-quantum crypto (liboqs)")
            testMilestone49();
        END_SET("M49: post-quantum crypto (liboqs)")

        START_SET("M52: ObjC ZarrProvider")
            testMilestone52();
        END_SET("M52: ObjC ZarrProvider")

        START_SET("M53: ObjC Bruker TDF reader")
            testMilestone53();
        END_SET("M53: ObjC Bruker TDF reader")

        START_SET("M59: ObjC imzML reader")
            testImzMLReader();
        END_SET("M59: ObjC imzML reader")

        START_SET("M60: ObjC mzTab reader")
            testMzTabReader();
        END_SET("M60: ObjC mzTab reader")

        START_SET("M62: ObjC stress + concurrency")
            testStress();
        END_SET("M62: ObjC stress + concurrency")

        START_SET("M63: ObjC Waters MassLynx reader")
            testWatersMassLynxReader();
        END_SET("M63: ObjC Waters MassLynx reader")

        START_SET("writeMinimal flat-buffer fast path")
            testWriteMinimal();
        END_SET("writeMinimal flat-buffer fast path")

        START_SET("imzML writer (v0.9+)")
            testImzMLWriter();
        END_SET("imzML writer (v0.9+)")

        START_SET("mzTab writer (v0.9+)")
            testMzTabWriter();
        END_SET("mzTab writer (v0.9+)")

        START_SET("M67: transport codec (v0.10)")
            testTransportCodec();
        END_SET("M67: transport codec (v0.10)")

        START_SET("S3: transport writer streaming sink protocol")
            testTransportWriterSink();
        END_SET("S3: transport writer streaming sink protocol")

        START_SET("S3: dataset walker + visitor protocol")
            testDatasetWalker();
        END_SET("S3: dataset walker + visitor protocol")

        START_SET("#140: dataset walker emits v0.11 prelude + genomic AUs")
            testDatasetWalkerV011Prelude();
        END_SET("#140: dataset walker emits v0.11 prelude + genomic AUs")

        START_SET("M67.5: transport ingest (incremental parser)")
            testTransportIngest();
        END_SET("M67.5: transport ingest (incremental parser)")

        START_SET("M68: transport client (WebSocket, v0.10)")
            testTransportClient();
        END_SET("M68: transport client (WebSocket, v0.10)")

        START_SET("M69: acquisition simulator (v0.10)")
            testAcquisitionSimulator();
        END_SET("M69: acquisition simulator (v0.10)")

        START_SET("M68.5: transport server (v0.10)")
            testTransportServer();
        END_SET("M68.5: transport server (v0.10)")

        START_SET("M70: bidirectional conformance (v0.10)")
            testTransportConformance();
        END_SET("M70: bidirectional conformance (v0.10)")

        START_SET("M71: selective access + protection metadata (v0.10)")
            testSelectiveAccess();
        END_SET("M71: selective access + protection metadata (v0.10)")

        START_SET("v1.0: per-AU encryption primitives")
            testPerAUEncryption();
        END_SET("v1.0: per-AU encryption primitives")

        START_SET("v1.0: per-AU encryption file round-trip")
            testPerAUFile();
            testPerAUFileDecryptInPlace();
        END_SET("v1.0: per-AU encryption file round-trip")

        START_SET("v1.0: encrypted transport emission")
            testEncryptedTransport();
        END_SET("v1.0: encrypted transport emission")

        START_SET("M73: Raman / IR spectra + imaging (v0.11)")
            testMilestone73();
        END_SET("M73: Raman / IR spectra + imaging (v0.11)")

        START_SET("M73.1: JCAMP-DX compression + UV-Vis + 2D-COS (v0.11.1)")
            testMilestone73_1();
        END_SET("M73.1: JCAMP-DX compression + UV-Vis + 2D-COS (v0.11.1)")

        START_SET("vibrational .tio round-trip (IR/Raman/UV-Vis materialize)")
            testJcampVibrationalRoundTrip();
        END_SET("vibrational .tio round-trip (IR/Raman/UV-Vis materialize)")

        START_SET("FD-1 Phase A: multi-recipient ProtectionMetadata carriage")
            testFD1MultiRecipient();
        END_SET("FD-1 Phase A: multi-recipient ProtectionMetadata carriage")

        START_SET("FD-1 Phase A-4: multi-recipient cross-language conformance")
            testMultiRecipientXLang();
        END_SET("FD-1 Phase A-4: multi-recipient cross-language conformance")

        START_SET("dek_wrapped envelope-encryption cross-language conformance")
            testDekWrappedXLang();
        END_SET("dek_wrapped envelope-encryption cross-language conformance")

        START_SET("FD-1 Phase C-0: standalone key-wrap primitive")
            testC0StandaloneKeyWrap();
        END_SET("FD-1 Phase C-0: standalone key-wrap primitive")

        START_SET("FD-1 Phase C-2a: server_kek_id in ProtectionMetadata")
            testC2aServerKekId();
        END_SET("FD-1 Phase C-2a: server_kek_id in ProtectionMetadata")

        START_SET("M76: JCAMP-DX compressed writer cross-language byte-parity")
            testM76JcampConformance();
        END_SET("M76: JCAMP-DX compressed writer cross-language byte-parity")

        START_SET("M77: 2D-COS compute primitives (v0.12.0)")
            testMilestone77();
        END_SET("M77: 2D-COS compute primitives (v0.12.0)")

        START_SET("M78: mzTab PEH/PEP + SFH/SMF + SEH/SME (v0.12.0)")
            testMilestone78();
        END_SET("M78: mzTab PEH/PEP + SFH/SMF + SEH/SME (v0.12.0)")

        START_SET("v1.1 parity: encrypt -> close -> reopen -> decrypt -> read")
            testV11EncryptionParity();
        END_SET("v1.1 parity: encrypt -> close -> reopen -> decrypt -> read")

        START_SET("v1.1.1 parity: decryptInPlaceAtPath:withKey:error:")
            testV111DecryptInPlace();
        END_SET("v1.1.1 parity: decryptInPlaceAtPath:withKey:error:")

        START_SET("M79: modality + genomic enums (v0.11)")
            testM79GenomicEnums();
        END_SET("M79: modality + genomic enums (v0.11)")

        START_SET("M82: GenomicRun + AlignedRead + signal channels (v0.11)")
            testM82GenomicRun();
        END_SET("M82: GenomicRun + AlignedRead + signal channels (v0.11)")

        START_SET("M83: rANS codec")
            testM83Rans();
        END_SET("M83: rANS codec")

        START_SET("M84: BASE_PACK codec")
            testM84BasePack();
        END_SET("M84: BASE_PACK codec")

        START_SET("Packed reference storage")
            testPackedReference();
        END_SET("Packed reference storage")

        START_SET("FLOAT_DELTA_ZSTD codec")
            testFloatDeltaZstd();
            testFloatDeltaZstdPlainTransform();
            testQualitiesV5();
        END_SET("FLOAT_DELTA_ZSTD codec")

        START_SET("M85: QUALITY_BINNED codec")
            testM85Quality();
        END_SET("M85: QUALITY_BINNED codec")

        // M85B NAME_TOKENIZED v1, M93 REF_DIFF v1
        // unit + pipeline, M94.Z V1 unit, and M94.Z V2 dispatch test
        // suites removed alongside the v1 codec impls.

        START_SET("M86: codec wiring")
            testM86GenomicCodecWiring();
        END_SET("M86: codec wiring")

        START_SET("M94.Z: CRAM-mimic FQZCOMP_NX16 throughput")
            testM94ZFqzcompPerf();
        END_SET("M94.Z: CRAM-mimic FQZCOMP_NX16 throughput")

        START_SET("Task 17: libttio_rans native backend introspection")
            testTtioRansNative();
        END_SET("Task 17: libttio_rans native backend introspection")

        START_SET("M94.Z V4 CRAM-mimic dispatch (Stage 3 Task 8)")
            testM94ZV4Dispatch();
        END_SET("M94.Z V4 CRAM-mimic dispatch (Stage 3 Task 8)")

        START_SET("M94.Z V4 cross-corpus byte-exact vs Python (Stage 3 Task 9)")
            testM94ZV4ByteExact();
        END_SET("M94.Z V4 cross-corpus byte-exact vs Python (Stage 3 Task 9)")

        START_SET("mate_info v2 codec round-trip + invalid-input")
            testMateInfoV2();
        END_SET("mate_info v2 codec round-trip + invalid-input")

        START_SET("v1.7 #11 Task 14: mate_info v2 ObjC writer/reader dispatch")
            testMateInfoV2Dispatch();
        END_SET("v1.7 #11 Task 14: mate_info v2 ObjC writer/reader dispatch")

        START_SET("ref_diff v2 codec round-trip + invalid-input")
            testRefDiffV2();
        END_SET("ref_diff v2 codec round-trip + invalid-input")

        START_SET("v1.8 #11 Task 14: ref_diff v2 ObjC writer/reader dispatch")
            testRefDiffV2Dispatch();
        END_SET("v1.8 #11 Task 14: ref_diff v2 ObjC writer/reader dispatch")

        START_SET("M97: read_role + QUALITY_BINNED guard + slice_bytes")
            testM97LongReadProfile();
        END_SET("M97: read_role + QUALITY_BINNED guard + slice_bytes")

        START_SET("M98: assembly graph GFA + storage round-trip")
            testM98AssemblyGraph();
        END_SET("M98: assembly graph GFA + storage round-trip")

        START_SET("M99: blocks_v1 per-AU streaming walkers")
            testM99BlocksV1PerAU();
        END_SET("M99: blocks_v1 per-AU streaming walkers")

        START_SET("M101: SAM optional tags")
            testM101SamTags();
        END_SET("M101: SAM optional tags")

        START_SET("name_tok v2 codec round-trip + invalid-input")
            testNameTokenizerV2();
        END_SET("name_tok v2 codec round-trip + invalid-input")

        START_SET("v1.8 #11 ch3 Task 14: name_tok v2 ObjC writer/reader dispatch")
            testNameTokenizedV2Dispatch();
        END_SET("v1.8 #11 ch3 Task 14: name_tok v2 ObjC writer/reader dispatch")

        START_SET("v1.10 #10: TTIOOffsetsFromLengths cumsum + uint32->uint64 overflow")
            testOffsetsCumsum();
        END_SET("v1.10 #10: TTIOOffsetsFromLengths cumsum + uint32->uint64 overflow")

        START_SET("PO1: TTIOGenomicIndex region/flag query equivalence (disk interned-id + in-memory fallback)")
            testGenomicIndexRegionQuery();
        END_SET("PO1: TTIOGenomicIndex region/flag query equivalence (disk interned-id + in-memory fallback)")

        START_SET("M95: DELTA_RANS_ORDER0 codec unit")
            testM95DeltaRansUnit();
        END_SET("M95: DELTA_RANS_ORDER0 codec unit")

        START_SET("M87: BAM importer")
            testM87BamImporter();
        END_SET("M87: BAM importer")

        START_SET("M88: CRAM/BAM round-trip")
            testM88CramBamRoundTrip();
        END_SET("M88: CRAM/BAM round-trip")

        START_SET("M89: GenomicRead AU + filter + transport (v0.11)")
            testM89GenomicTransport();
        END_SET("M89: GenomicRead AU + filter + transport (v0.11)")

        START_SET("M90: genomic protection (per-AU + signatures + anon + region)")
            testM90GenomicProtection();
        END_SET("M90: genomic protection (per-AU + signatures + anon + region)")

        START_SET("M90.8/M90.9/M90.10: enc-transport genomic + AU compound + wire codec")
            testM90Parity();
        END_SET("M90.8/M90.9/M90.10: enc-transport genomic + AU compound + wire codec")

        START_SET("M90.11/12/13/14/15: encrypted headers + uint8 MPAD + SAM-overlap + seeded RNG + chromosomes sign")
            testM90Final();
        END_SET("M90.11/12/13/14/15: encrypted headers + uint8 MPAD + SAM-overlap + seeded RNG + chromosomes sign")

        START_SET("V4: edge case hardening")
            testV4EdgeCases();
        END_SET("V4: edge case hardening")

        START_SET("V8: HDF5 corruption recovery")
            testV8Hdf5Corruption();
        END_SET("V8: HDF5 corruption recovery")

        START_SET("C1: CLI tools coverage")
            testC1ToolsCli();
        END_SET("C1: CLI tools coverage")

        START_SET("C2: HDF5 error-path coverage")
            testC2HDF5ErrorPaths();
        END_SET("C2: HDF5 error-path coverage")

        START_SET("C3: providers error-path coverage")
            testC3ProvidersErrorPaths();
        END_SET("C3: providers error-path coverage")

        START_SET("C5: protection package gap")
            testC5ProtectionGap();
        END_SET("C5: protection package gap")

        START_SET("C3b: providers write-path coverage")
            testC3bProvidersWritePaths();
        END_SET("C3b: providers write-path coverage")

        START_SET("Providers: extendable datasets and UInt64")
            testProviderExtendable();
        END_SET("Providers: extendable datasets and UInt64")

        START_SET("Genomics: storage-path writer and blocks")
            testGenomicBlocks();
        END_SET("Genomics: storage-path writer and blocks")

        START_SET("Genomics: stream writer, blocks_v1 layout, lazy reference")
            testGenomicStreamWriter();
        END_SET("Genomics: stream writer, blocks_v1 layout, lazy reference")

        START_SET("Genomics: blocks_v1 golden fixture")
            testBlocksV1Golden();
        END_SET("Genomics: blocks_v1 golden fixture")

        START_SET("Genomics: blocks_v1 reader, signatures, transport")
            testGenomicBlocksReader();
        END_SET("Genomics: blocks_v1 reader, signatures, transport")

        START_SET("Threads: knob, pool, block-parallel writer")
            testThreads();
            testGenomicStreamWriterThreads();
            testSpectralStreamWriterThreads();
            testFastqRecordScanner();
        END_SET("Threads: knob, pool, block-parallel writer")

        START_SET("spectral iterBlocks")
            testSpectralIterBlocks();
        END_SET("spectral iterBlocks")

        START_SET("Codecs: reference resolver, external FASTA")
            testReferenceResolverExternal();
        END_SET("Codecs: reference resolver, external FASTA")

        START_SET("Spectral: FDZ1 block API")
            testFloatDeltaZstdBlocks();
        END_SET("Spectral: FDZ1 block API")

        START_SET("Spectral: stream writer and range reads")
            testSpectralStreamWriter();
        END_SET("Spectral: stream writer and range reads")

        START_SET("Importers: streaming BAM, FASTQ, mzML")
            testStreamingImporters();
        END_SET("Importers: streaming BAM, FASTQ, mzML")

        START_SET("Exporters: streaming BAM, FASTQ, mzML")
            testStreamingExporters();
        END_SET("Exporters: streaming BAM, FASTQ, mzML")

        START_SET("C2b: HDF5CompoundType coverage")
            testC2bHDF5CompoundType();
        END_SET("C2b: HDF5CompoundType coverage")

        START_SET("C3c: CanonicalBytes + SqliteProvider")
            testC3cCanonicalBytes();
        END_SET("C3c: CanonicalBytes + SqliteProvider")

        START_SET("Phase 1+2: TTIORun protocol + mixed-dict write API")
            testPhase12RunProtocol();
        END_SET("Phase 1+2: TTIORun protocol + mixed-dict write API")

        START_SET("Task 30: MS runs via memory/sqlite/zarr provider URL")
            testTask30MSProviderURL();
        END_SET("Task 30: MS runs via memory/sqlite/zarr provider URL")

        START_SET("Task 31: instance writer (-writeToFilePath:) via memory/sqlite/zarr URL")
            testTask31InstanceWriterParity();
        END_SET("Task 31: instance writer (-writeToFilePath:) via memory/sqlite/zarr URL")

        START_SET("FASTA / FASTQ I/O")
            testFastaFastqIo();
        END_SET("FASTA / FASTQ I/O")

        START_SET("TTIOProgressSink + Stage B reader hooks")
            testProgressSink();
        END_SET("TTIOProgressSink + Stage B reader hooks")

        START_SET("TTIOProgressSink + Stage C reader hooks")
            testProgressSinkStageC();
        END_SET("TTIOProgressSink + Stage C reader hooks")

        START_SET("TTIOProgressSink + Stage D writer hooks")
            testProgressSinkStageD();
        END_SET("TTIOProgressSink + Stage D writer hooks")

        START_SET("FASTA / FASTQ -> .tio -> FASTA / FASTQ round-trip")
            testFastaFastqTioRoundTrip();
        END_SET("FASTA / FASTQ -> .tio -> FASTA / FASTQ round-trip")

        START_SET("1.1.0: SpectralDataset.references read-back (Phase 0 tio-browser)")
            testReferencesAccessor();
        END_SET("1.1.0: SpectralDataset.references read-back (Phase 0 tio-browser)")

        START_SET("1.1.0: ReferenceImport.writeToDataset round-trip (Phase 0 tio-browser)")
            testReferenceImportWriteToDataset();
        END_SET("1.1.0: ReferenceImport.writeToDataset round-trip (Phase 0 tio-browser)")

        START_SET("M100: SpectralDataset writable open (readFromFilePath:writable:)")
            testSpectralDatasetWritableOpen();
        END_SET("M100: SpectralDataset writable open (readFromFilePath:writable:)")

        START_SET("1.2.0: TTIOMSImage.mzAxis property + persistence + pixelSpectra")
            testMSImageMzAxis();
        END_SET("1.2.0: TTIOMSImage.mzAxis property + persistence + pixelSpectra")

        START_SET("1.2.0: TTIOMSImage legacy (no mzAxis) round-trip")
            testMSImageLegacyMzAxisAbsent();
        END_SET("1.2.0: TTIOMSImage legacy (no mzAxis) round-trip")

        START_SET("1.2.0: TTIOSpectralDataset.msImage accessor (category)")
            testSpectralDatasetMsImageAccessor();
        END_SET("1.2.0: TTIOSpectralDataset.msImage accessor (category)")

        START_SET("1.2.0: TTIOMSImage pixelSpectra raises on missing mzAxis")
        {
            testPixelSpectraRaisesWhenMzAxisAbsent();
        }
        END_SET("1.2.0: TTIOMSImage pixelSpectra raises on missing mzAxis")

        START_SET("transport v0.11: packet-type constants (0x10-0x1B) + feature flag")
            testTransportPacketTypeV011();
        END_SET("transport v0.11: packet-type constants (0x10-0x1B) + feature flag")

        START_SET("transport v0.11: reader skips unknown packet types (forward compat)")
            testTransportReaderSkipUnknown();
        END_SET("transport v0.11: reader skips unknown packet types (forward compat)")

        START_SET("transport v0.11 Task 3.2: writeReferenceGroup: emits 0x10/0x11/0x12")
            testTransportWriterReference();
        END_SET("transport v0.11 Task 3.2: writeReferenceGroup: emits 0x10/0x11/0x12")

        START_SET("transport v0.11 Task 3.3: reader decodes 0x10/0x11/0x12 into SpectralDataset")
            testTransportReaderReference();
        END_SET("transport v0.11 Task 3.3: reader decodes 0x10/0x11/0x12 into SpectralDataset")

        START_SET("transport v0.11: TTIOArrowIpcCodec round-trip (IDENTIFICATIONS_TABLE 0x16 + QUANTIFICATIONS_TABLE 0x17)")
            testArrowIpcCodec();
        END_SET("transport v0.11: TTIOArrowIpcCodec round-trip (IDENTIFICATIONS_TABLE 0x16 + QUANTIFICATIONS_TABLE 0x17)")

        START_SET("transport v0.11 Task 3.4: ENCRYPTION_ALGORITHM (0x1B) writer + reader")
            testTransportEncryptionAlgorithm();
        END_SET("transport v0.11 Task 3.4: ENCRYPTION_ALGORITHM (0x1B) writer + reader")

        START_SET("transport v0.11 Task 3.5: DATASET_PROVENANCE (0x18) writer + reader")
            testTransportDatasetProvenance();
        END_SET("transport v0.11 Task 3.5: DATASET_PROVENANCE (0x18) writer + reader")

        START_SET("transport v0.11 Task 3.6: IMAGE_HEADER + IMAGE_PIXEL + END_OF_IMAGE (0x13-0x15) writer + reader")
            testTransportImage();
        END_SET("transport v0.11 Task 3.6: IMAGE_HEADER + IMAGE_PIXEL + END_OF_IMAGE (0x13-0x15) writer + reader")

        START_SET("transport v0.11 Task 5.1: IMAGE_PIXEL processed-mode (sparse) writer + reader")
            testTransportImageProcessed();
        END_SET("transport v0.11 Task 5.1: IMAGE_PIXEL processed-mode (sparse) writer + reader")

        START_SET("transport v0.11 Task 5.2: TTIOSpectralDataset -irImage / -ramanImage accessors")
            testSpectralDatasetIRImage();
        END_SET("transport v0.11 Task 5.2: TTIOSpectralDataset -irImage / -ramanImage accessors")

        START_SET("P2.5 OIT2: TTIOSpectralDataset imageForKind: / images uniform accessors")
            testSpectralDatasetImages();
        END_SET("P2.5 OIT2: TTIOSpectralDataset imageForKind: / images uniform accessors")

        START_SET("transport v0.11 Task 5.6 follow-up: TTIOIRImage typed HDF5 attrs (ir_mode i64) for cross-lang parity")
            testIRImageWriteReadTypedAttrs();
        END_SET("transport v0.11 Task 5.6 follow-up: TTIOIRImage typed HDF5 attrs (ir_mode i64) for cross-lang parity")

        START_SET("transport v0.11 Task 5.3: IMAGE modality dispatch — Raman (modality=1)")
            testTransportRamanImage();
        END_SET("transport v0.11 Task 5.3: IMAGE modality dispatch — Raman (modality=1)")

        START_SET("transport v0.11 Task 5.3: IMAGE modality dispatch — IR (modality=2)")
            testTransportIRImage();
        END_SET("transport v0.11 Task 5.3: IMAGE modality dispatch — IR (modality=2)")

        START_SET("transport v0.11 Task 3.7: IDENTIFICATIONS_TABLE (0x16) + QUANTIFICATIONS_TABLE (0x17) writer + reader")
            testTransportIdentificationsQuantifications();
        END_SET("transport v0.11 Task 3.7: IDENTIFICATIONS_TABLE (0x16) + QUANTIFICATIONS_TABLE (0x17) writer + reader")

        START_SET("transport v0.11 Task 3.9: -writeDataset: emits full v0.11 prelude per §5.4")
            testTransportWriteDatasetPrelude();
        END_SET("transport v0.11 Task 3.9: -writeDataset: emits full v0.11 prelude per §5.4")

        START_SET("Stage 6.4: TTIOSubject + TTIOSample validation + attributesJson byte parity")
            testSubjectSample();
        END_SET("Stage 6.4: TTIOSubject + TTIOSample validation + attributesJson byte parity")

        START_SET("Stage 6.4: TTIOSpectralDataset subjects + samples HDF5 round-trip + soft-FK warning")
            testSpectralDatasetSubjectsSamples();
        END_SET("Stage 6.4: TTIOSpectralDataset subjects + samples HDF5 round-trip + soft-FK warning")

        START_SET("Stage 6.4: TTIOArrowIpcCodec subjects + samples encode/decode (0x19 + 0x1A null convention)")
            testArrowIpcCodecSubjectsSamples();
        END_SET("Stage 6.4: TTIOArrowIpcCodec subjects + samples encode/decode (0x19 + 0x1A null convention)")

        START_SET("Stage 6.4: SUBJECT_METADATA (0x19) + SAMPLE_METADATA (0x1A) writer + reader wire round-trip")
            testTransportSubjectsSamples();
        END_SET("Stage 6.4: SUBJECT_METADATA (0x19) + SAMPLE_METADATA (0x1A) writer + reader wire round-trip")

        START_SET("transport v0.11 Task 3.10: AccessorMatrixConformance (8 first-class accessors)")
            testAccessorMatrixConformance();
        END_SET("transport v0.11 Task 3.10: AccessorMatrixConformance (8 first-class accessors)")

        START_SET("transport v0.11 Task 3.10: CoverageGapWatchdog (.tis >= 1% of .tio + every accessor preserved)")
            testCoverageGapWatchdog();
        END_SET("transport v0.11 Task 3.10: CoverageGapWatchdog (.tis >= 1% of .tio + every accessor preserved)")

        START_SET("codec registry")
            testCodecRegistry();
        END_SET("codec registry")

        START_SET("OT1: TTIOImportedDataset draft + write-through delegate")
            testImportedDataset();
        END_SET("OT1: TTIOImportedDataset draft + write-through delegate")

        START_SET("OT3: TTIORunSelection analytical/NMR/genomic run selection")
            testRunSelection();
        END_SET("OT3: TTIORunSelection analytical/NMR/genomic run selection")

        START_SET("OT4: TTIOBrukerTDFReader readDataset write-through draft")
            testBrukerReaderDataset();
        END_SET("OT4: TTIOBrukerTDFReader readDataset write-through draft")

        START_SET("OT5: per-format TTIOReader adapters")
            testFormatReaders();
        END_SET("OT5: per-format TTIOReader adapters")

        START_SET("OT6: per-format TTIOWriter adapters")
            testFormatWriters();
        END_SET("OT6: per-format TTIOWriter adapters")

        START_SET("OT7: importer + exporter registries")
            testImporterExporterRegistry();
        END_SET("OT7: importer + exporter registries")

        START_SET("OT7: registry cross-language parity fence")
            testRegistryParity();
        END_SET("OT7: registry cross-language parity fence")

        START_SET("OT8: TtioEncode / TtioExport CLI tools")
            testEncodeExportCli();
        END_SET("OT8: TtioEncode / TtioExport CLI tools")

        START_SET("OIT1: TTIOImage base extraction (MS/Raman/IR round-trip + kind/spectralAxis)")
            testImageBase();
        END_SET("OIT1: TTIOImage base extraction (MS/Raman/IR round-trip + kind/spectralAxis)")
    }
    return 0;
}
