/*
 * TTI-O Java Implementation
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.tools;

import global.thalion.ttio.importers.ImporterRegistry;
import global.thalion.ttio.importers.UnknownFormatError;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Encode one or more input files of a given source format into a
 * {@code .tio} container, dispatching through {@link ImporterRegistry}.
 *
 * <p>Parallel to Python {@code ttio.tools.workbench_cli} {@code cmd_encode}.
 *
 * <p>Usage:
 * <pre>
 *   EncodeCli --input &lt;path&gt; [--input &lt;path&gt; ...] \
 *             --format &lt;fmt&gt; --output &lt;out.tio&gt; [--extra k=v ...]
 *   EncodeCli --list-formats
 * </pre>
 *
 * <p>Exit codes (mirroring Python {@code cmd_encode}): {@code 0} success,
 * {@code 2} importer failure or bad/missing args, {@code 3} unsupported or
 * delegated ({@code fasta}/{@code fastq}) format.
 */
public final class EncodeCli {

    private static final String PROG = "EncodeCli";

    private EncodeCli() {
    }

    public static void main(String[] args) {
        System.exit(run(args));
    }

    /** Testable entry point: returns the process exit code. */
    public static int run(String[] args) {
        List<String> inputs = new ArrayList<>();
        String rawFormat = null;
        String output = null;
        Map<String, Object> opts = new LinkedHashMap<>();

        for (int i = 0; i < args.length; i++) {
            String a = args[i];
            switch (a) {
                case "--list-formats" -> {
                    for (String f : ImporterRegistry.supportedEncodeFormats()) {
                        System.out.println(f);
                    }
                    return 0;
                }
                case "--threads" -> System.setProperty("ttio.threads", args[++i]);
                case "--input" -> {
                    if (i + 1 >= args.length) {
                        return usage();
                    }
                    inputs.add(args[++i]);
                }
                case "--format" -> {
                    if (i + 1 >= args.length) {
                        return usage();
                    }
                    rawFormat = args[++i];
                }
                case "--output" -> {
                    if (i + 1 >= args.length) {
                        return usage();
                    }
                    output = args[++i];
                }
                case "--extra" -> {
                    if (i + 1 >= args.length) {
                        return usage();
                    }
                    String kv = args[++i];
                    int eq = kv.indexOf('=');
                    if (eq < 0) {
                        System.err.println(
                            PROG + ": --extra expects key=value, got: " + kv);
                        return usage();
                    }
                    opts.put(kv.substring(0, eq), kv.substring(eq + 1));
                }
                default -> {
                    System.err.println(PROG + ": unknown argument: " + a);
                    return usage();
                }
            }
        }

        if (inputs.isEmpty() || rawFormat == null || output == null) {
            return usage();
        }

        String fmt = ImporterRegistry.normalize(rawFormat);

        // fasta / fastq keep their richer dedicated CLIs; not yet wired
        // through this command (Python delegates them — here we report).
        if (ImporterRegistry.CLI_DELEGATED.contains(fmt)) {
            System.err.println(PROG + ": format '" + fmt
                + "' is delegated to the dedicated fasta/fastq tools");
            return 3;
        }

        if (!ImporterRegistry.isRegistryFormat(fmt)) {
            return unsupported(rawFormat);
        }

        try {
            ImporterRegistry.encode(fmt, inputs, Path.of(output), opts);
        } catch (UnknownFormatError e) {
            return unsupported(rawFormat);
        } catch (Exception e) {
            System.err.println(PROG + ": encode failed (" + rawFormat + "): "
                + e.getMessage());
            return 2;
        }
        System.out.println("encoded " + output);
        return 0;
    }

    private static int unsupported(String rawFormat) {
        System.err.println(PROG + ": unsupported --format: " + rawFormat
            + ". Supported: "
            + String.join(", ", ImporterRegistry.supportedEncodeFormats())
            + ".");
        return 3;
    }

    private static int usage() {
        System.err.println(
            "usage: EncodeCli --input <path> [--input <path> ...] "
            + "--format <fmt> --output <out.tio> [--extra k=v ...]\n"
            + "       EncodeCli --list-formats\n"
            + "  streaming extras (bam/sam/cram, mzml): block_reads=N block_bytes=N "
            + "legacy_whole_channel=1 reference=<fasta> embed_reference=1 "
            + "batch_reads=N batch_spectra=N group_reads=1");
        return 2;
    }
}
