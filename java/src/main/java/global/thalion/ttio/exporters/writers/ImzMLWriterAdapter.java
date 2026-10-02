/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.exporters.writers;

import global.thalion.ttio.AcquisitionRun;
import global.thalion.ttio.Enums;
import global.thalion.ttio.MSImage;
import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.exporters.ImzMLWriter;
import global.thalion.ttio.exporters.Writer;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;

/**
 * {@link Writer} adapter for imzML export. Mirrors Python
 * {@code ttio.exporters.writers.ImzMLWriter} (null-guard message + sibling
 * {@code .ibd} via {@code with_suffix}) and the tio-browser GUI
 * {@code ExportTask.exportImzML} call shape (mode + grid + pixel-size +
 * scan-pattern arguments from the {@link MSImage}).
 *
 * <p>Mode default {@code "continuous"} matches the GUI {@code ExportConfig}
 * fallback; override via the {@code "mode"} opt.</p>
 *
 * <p>(M102) A dataset without an MS image is exported from its first
 * mass-spectrum run that carries pixel coordinates (an imported
 * processed-mode file's {@code imzml_pixels} run, or an older file's
 * legacy coordinate parameter) via {@link ImzMLWriter#writeFromRun}; the
 * mode then defaults to the recorded {@code imzml_mode}, else to what the
 * m/z arrays allow.</p>
 *
 * @since 1.7.0
 */
public final class ImzMLWriterAdapter implements Writer {

    @Override
    public void write(SpectralDataset ds, String layer, Path output,
                      Map<String, Object> opts) throws IOException {
        MSImage img = (MSImage) ds.imageForKind(Enums.ImageKind.MS);
        // Python: ibd = Path(output).with_suffix(".ibd").
        Path ibd = withSuffixIbd(output);
        Object modeOpt = opts != null ? opts.get("mode") : null;
        if (img == null) {
            AcquisitionRun pixelRun = ImzMLWriter.findPixelRun(ds);
            if (pixelRun == null) {
                throw new IllegalArgumentException(
                    "dataset has no MS image to export as imzML");
            }
            ImzMLWriter.writeFromRun(pixelRun, ds.provenanceRecords(),
                output, ibd, modeOpt != null ? modeOpt.toString() : null, null);
            return;
        }
        String mode = modeOpt != null ? modeOpt.toString() : "continuous";
        // GUI: ImzMLWriter.write(img.toPixelSpectra(), targetPath, ibd, mode,
        //          width, height, 1, pixelSizeX, pixelSizeY, scanPattern,
        //          uuidHex=null, sink).
        ImzMLWriter.write(
            img.toPixelSpectra(),
            output,
            ibd,
            mode,
            img.width(), img.height(), 1,
            img.pixelSizeX(), img.pixelSizeY(),
            img.scanPattern() != null ? img.scanPattern() : "flyback",
            null,
            null);
    }

    /** Replace the output file's extension with {@code .ibd} (Python's
     *  {@code Path.with_suffix(".ibd")}). */
    private static Path withSuffixIbd(Path output) {
        String name = output.getFileName().toString();
        int dot = name.lastIndexOf('.');
        String base = dot >= 0 ? name.substring(0, dot) : name;
        Path parent = output.getParent();
        Path ibdName = Path.of(base + ".ibd");
        return parent != null ? parent.resolve(ibdName) : ibdName;
    }
}
