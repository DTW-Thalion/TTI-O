/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.importers;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * imzML arrays declared through referenceableParamGroups.
 *
 * <p>Many writers (the HR2MSI PXD001283 files among them) state each
 * binaryDataArray's kind, precision and compression in a
 * {@code <referenceableParamGroup>} and reference it with
 * {@code <referenceableParamGroupRef>}, keeping only the external
 * offset/length cvParams inline. The reader must resolve the reference
 * to know which array an offset belongs to. The fixture declares m/z
 * 64-bit and intensity 32-bit only in the groups, so the precision is
 * checked as well as the array kind.
 *
 * <p>Cross-language counterparts:
 *   python/tests/test_imzml_param_groups.py
 *   objc/Tests/TestImzMLParamGroups.m
 */
final class ImzMLParamGroupsTest {

    private static final String UUID_HEX = "0123456789abcdef0123456789abcdef";
    private static final int N_PIXELS = 3;

    private static double mzAt(int pixel, int j) { return 100.0 + 50.0 * j + pixel; }
    private static float intensityAt(int pixel, int j) { return 10.0f * (pixel + 1) + j; }

    private static Path writeFixture(Path dir) throws IOException {
        ByteArrayOutputStream ibd = new ByteArrayOutputStream();
        ibd.write(java.util.HexFormat.of().parseHex(UUID_HEX));

        StringBuilder spectra = new StringBuilder();
        for (int p = 0; p < N_PIXELS; p++) {
            int n = 3 + p;
            int mzOffset = ibd.size();
            ByteBuffer mz = ByteBuffer.allocate(n * 8).order(ByteOrder.LITTLE_ENDIAN);
            for (int j = 0; j < n; j++) mz.putDouble(mzAt(p, j));
            ibd.write(mz.array());
            int intOffset = ibd.size();
            ByteBuffer in = ByteBuffer.allocate(n * 4).order(ByteOrder.LITTLE_ENDIAN);
            for (int j = 0; j < n; j++) in.putFloat(intensityAt(p, j));
            ibd.write(in.array());

            spectra.append("<spectrum id=\"s").append(p).append("\" index=\"").append(p)
                   .append("\" defaultArrayLength=\"0\">")
                   .append("<scanList count=\"1\"><scan>")
                   .append("<cvParam cvRef=\"IMS\" accession=\"IMS:1000050\" name=\"position x\" value=\"")
                   .append(p + 1).append("\"/>")
                   .append("<cvParam cvRef=\"IMS\" accession=\"IMS:1000051\" name=\"position y\" value=\"1\"/>")
                   .append("</scan></scanList>")
                   .append("<binaryDataArrayList count=\"2\">")
                   .append(array("mzArray", mzOffset, n, n * 8))
                   .append(array("intensityArray", intOffset, n, n * 4))
                   .append("</binaryDataArrayList></spectrum>");
        }

        String xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
            + "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1\">"
            + "<fileDescription><fileContent>"
            + "<cvParam cvRef=\"IMS\" accession=\"IMS:1000031\" name=\"processed\"/>"
            + "<cvParam cvRef=\"IMS\" accession=\"IMS:1000080\" name=\"universally unique identifier\" value=\"{"
            + UUID_HEX + "}\"/>"
            + "</fileContent></fileDescription>"
            + "<referenceableParamGroupList count=\"2\">"
            + "<referenceableParamGroup id=\"mzArray\">"
            + "<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"/>"
            + "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>"
            + "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>"
            + "</referenceableParamGroup>"
            + "<referenceableParamGroup id=\"intensityArray\">"
            + "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>"
            + "<cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\"/>"
            + "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>"
            + "</referenceableParamGroup>"
            + "</referenceableParamGroupList>"
            + "<scanSettingsList count=\"1\"><scanSettings id=\"ss\">"
            + "<cvParam cvRef=\"IMS\" accession=\"IMS:1000042\" name=\"max count of pixels x\" value=\""
            + N_PIXELS + "\"/>"
            + "<cvParam cvRef=\"IMS\" accession=\"IMS:1000043\" name=\"max count of pixels y\" value=\"1\"/>"
            + "</scanSettings></scanSettingsList>"
            + "<run id=\"r\"><spectrumList count=\"" + N_PIXELS + "\">"
            + spectra
            + "</spectrumList></run></mzML>";

        Path imzml = dir.resolve("grouped.imzML");
        Files.writeString(imzml, xml, StandardCharsets.UTF_8);
        Files.write(dir.resolve("grouped.ibd"), ibd.toByteArray());
        return imzml;
    }

    private static String array(String ref, int offset, int length, int encoded) {
        return "<binaryDataArray encodedLength=\"0\">"
            + "<referenceableParamGroupRef ref=\"" + ref + "\"/>"
            + "<cvParam cvRef=\"IMS\" accession=\"IMS:1000102\" name=\"external offset\" value=\""
            + offset + "\"/>"
            + "<cvParam cvRef=\"IMS\" accession=\"IMS:1000103\" name=\"external array length\" value=\""
            + length + "\"/>"
            + "<cvParam cvRef=\"IMS\" accession=\"IMS:1000104\" name=\"external encoded length\" value=\""
            + encoded + "\"/>"
            + "<binary/></binaryDataArray>";
    }

    @Test
    void arraysDeclaredByParamGroupRef(@TempDir Path tmp) throws IOException {
        ImzMLReader.ImzMLImport imp = ImzMLReader.read(writeFixture(tmp));
        assertEquals("processed", imp.mode());
        assertEquals(UUID_HEX, imp.uuidHex());
        List<ImzMLReader.PixelSpectrum> spectra = imp.spectra();
        assertEquals(N_PIXELS, spectra.size());
        for (int p = 0; p < N_PIXELS; p++) {
            ImzMLReader.PixelSpectrum s = spectra.get(p);
            assertEquals(p + 1, s.x());
            assertEquals(1, s.y());
            int n = 3 + p;
            assertEquals(n, s.mz().length);
            assertEquals(n, s.intensity().length);
            for (int j = 0; j < n; j++) {
                assertEquals(mzAt(p, j), s.mz()[j], 0.0);
                assertEquals(intensityAt(p, j), s.intensity()[j], 0.0);
            }
        }
    }
}
