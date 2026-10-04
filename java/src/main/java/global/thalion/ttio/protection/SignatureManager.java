/* TTI-O Java Implementation / Copyright (c) 2026 The Thalion Initiative / SPDX-License-Identifier: LGPL-3.0-or-later */
package global.thalion.ttio.protection;

import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.security.GeneralSecurityException;
import java.security.MessageDigest;
import java.util.Base64;
import java.util.LinkedHashMap;
import java.util.Map;
import javax.crypto.Mac;
import javax.crypto.spec.SecretKeySpec;

/**
 * HMAC-SHA256 signatures with v2 canonical little-endian format for TTI-O datasets.
 *
 * <p><b>API status:</b> Stable.</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOSignatureManager}, Python {@code ttio.signatures}.</p>
 *
 *
 */
public final class SignatureManager {

    private static final String HMAC_ALGORITHM = "HmacSHA256";
    private static final int KEY_BYTES = 32;
    private static final String V2_PREFIX = "v2:";
    /** Post-quantum signature prefix (ML-DSA-87). */
    public static final String V3_PREFIX = "v3:";

    private SignatureManager() {}

    // -------------------------------------------------------------- sign / verify

    /**
     * Sign a byte array.
     *
     * @param data bytes to sign
     * @param key  32-byte HMAC key
     * @return signature string in the form {@code "v2:" + base64(hmac)}
     */
    public static String sign(byte[] data, byte[] key) {
        byte[] mac = hmac(data, key);
        return V2_PREFIX + Base64.getEncoder().encodeToString(mac);
    }

    /**
     * Verify a signature string against data and key.
     * Handles both {@code "v2:"}-prefixed signatures and unprefixed v1 signatures.
     *
     * @param data      original bytes
     * @param signature signature string to verify
     * @param key       32-byte HMAC key
     * @return {@code true} if the signature is valid
     */
    public static boolean verify(byte[] data, String signature, byte[] key) {
        byte[] expected = hmac(data, key);
        byte[] actual;
        if (signature.startsWith(V2_PREFIX)) {
            actual = Base64.getDecoder().decode(signature.substring(V2_PREFIX.length()));
        } else {
            // v1 compatibility: treat entire string as base64-encoded HMAC
            actual = Base64.getDecoder().decode(signature);
        }
        return MessageDigest.isEqual(expected, actual);
    }

    // -------------------------------------------------- algorithm-dispatched

    /**
     * Sign {@code data} with the named signature algorithm..
     *
     * <ul>
     *   <li>{@code "hmac-sha256"} — {@code key} is a 32-byte HMAC
     *       secret; output is {@code "v2:" + base64(hmac)}.</li>
     *   <li>{@code "ml-dsa-87"} — {@code key} is the 4896-byte
     *       ML-DSA-87 signing key; output is
     *       {@code "v3:" + base64(signature)}.</li>
     * </ul>
     *
     * @throws CipherSuite.UnsupportedAlgorithmException for unknown
     *         algorithm names.
     *
     */
    public static String sign(byte[] data, byte[] key, String algorithm) {
        if ("hmac-sha256".equals(algorithm)) {
            CipherSuite.validateKey(algorithm, key);
            return sign(data, key);
        }
        if ("ml-dsa-87".equals(algorithm)) {
            CipherSuite.validatePrivateKey(algorithm, key);
            byte[] sig = PostQuantumCrypto.sigSign(key, data);
            return V3_PREFIX + Base64.getEncoder().encodeToString(sig);
        }
        throw new CipherSuite.UnsupportedAlgorithmException(
            algorithm + ": signature path not implemented");
    }

    /**
     * Verify {@code signature} against {@code data}. Accepts v2: HMAC
     * and v3: ML-DSA-87. The algorithm parameter must match the
     * stored prefix — mismatches raise
     * {@link CipherSuite.UnsupportedAlgorithmException} so callers do
     * not silently pass verification of a file signed with a different
     * primitive.
     *
     * @param data       original bytes
     * @param signature  stored signature string (with prefix)
     * @param key        HMAC key (32 bytes) or ML-DSA-87 public key
     *                   (2592 bytes), depending on {@code algorithm}.
     *
     */
    public static boolean verify(byte[] data, String signature, byte[] key,
                                  String algorithm) {
        if (signature != null && signature.startsWith(V3_PREFIX)) {
            if (!"ml-dsa-87".equals(algorithm)) {
                throw new CipherSuite.UnsupportedAlgorithmException(
                    "stored signature is v3 (ml-dsa-87) but caller "
                    + "passed algorithm=" + algorithm);
            }
            CipherSuite.validatePublicKey(algorithm, key);
            byte[] sig = Base64.getDecoder().decode(
                signature.substring(V3_PREFIX.length()));
            return PostQuantumCrypto.sigVerify(key, data, sig);
        }
        if ("ml-dsa-87".equals(algorithm)) {
            throw new CipherSuite.UnsupportedAlgorithmException(
                "stored signature is not v3 (ml-dsa-87) — pass "
                + "algorithm=\"hmac-sha256\" to verify legacy signatures");
        }
        if ("hmac-sha256".equals(algorithm)) {
            CipherSuite.validateKey(algorithm, key);
            return verify(data, signature, key);
        }
        throw new CipherSuite.UnsupportedAlgorithmException(
            algorithm + ": verification path not implemented");
    }

    // ----------------------------------------------------------------- hmac

    /**
     * Compute raw HMAC-SHA256 bytes.
     *
     * @param data bytes to authenticate
     * @param key  32-byte HMAC key
     * @return 32-byte HMAC digest
     */
    public static byte[] hmac(byte[] data, byte[] key) {
        try {
            Mac mac = Mac.getInstance(HMAC_ALGORITHM);
            mac.init(new SecretKeySpec(key, HMAC_ALGORITHM));
            return mac.doFinal(data);
        } catch (GeneralSecurityException e) {
            throw new RuntimeException("HMAC-SHA256 computation failed", e);
        }
    }

    // -------------------------------------------------------- canonicalisation

    /**
     * Canonicalise a {@code double[]} to little-endian bytes for signing.
     *
     * @param data signal samples
     * @return little-endian byte representation
     */
    public static byte[] canonicalBytes(double[] data) {
        ByteBuffer buf = ByteBuffer.allocate(data.length * Double.BYTES)
                .order(ByteOrder.LITTLE_ENDIAN);
        for (double v : data) {
            buf.putDouble(v);
        }
        return buf.array();
    }

    /**
     * Canonicalise an {@code int[]} to little-endian bytes for signing.
     *
     * @param data integer samples
     * @return little-endian byte representation
     */
    public static byte[] canonicalBytes(int[] data) {
        ByteBuffer buf = ByteBuffer.allocate(data.length * Integer.BYTES)
                .order(ByteOrder.LITTLE_ENDIAN);
        for (int v : data) {
            buf.putInt(v);
        }
        return buf.array();
    }

    /**
     * Canonicalise a {@code long[]} to little-endian bytes for signing.
     *
     * @param data long samples
     * @return little-endian byte representation
     */
    public static byte[] canonicalBytes(long[] data) {
        ByteBuffer buf = ByteBuffer.allocate(data.length * Long.BYTES)
                .order(ByteOrder.LITTLE_ENDIAN);
        for (long v : data) {
            buf.putLong(v);
        }
        return buf.array();
    }

    /**
     * Canonicalise a variable-length string for compound signing.
     * Format: {@code uint32_le(len) || utf8_bytes}. A {@code null} string
     * produces a 4-byte zero-length prefix with no payload.
     *
     * @param s string to canonicalise (may be {@code null})
     * @return canonical byte representation
     */
    public static byte[] canonicalStringBytes(String s) {
        byte[] utf8 = (s != null) ? s.getBytes(StandardCharsets.UTF_8) : new byte[0];
        ByteBuffer buf = ByteBuffer.allocate(Integer.BYTES + utf8.length)
                .order(ByteOrder.LITTLE_ENDIAN);
        buf.putInt(utf8.length);
        buf.put(utf8);
        return buf.array();
    }

    // ------------------------------------------------------------ test key

    /**
     * Canonical signature test key: {@code key[i] = (byte)((0x5A ^ (i * 7)) & 0xFF)}
     * for i in [0,32).
     *
     * @return 32-byte deterministic test key
     */
    public static byte[] testKey() {
        byte[] key = new byte[KEY_BYTES];
        for (int i = 0; i < KEY_BYTES; i++) {
            key[i] = (byte) ((0x5A ^ (i * 7)) & 0xFF);
        }
        return key;
    }

    // ────────────────────────────────────────────── M90.2 genomic runs

    /** Channels signed by {@link #signGenomicRun}. M101: the SAM tags
     *  channel (a flat uint8 SAM_TAGS stream) is signed with the bases
     *  it describes, when present. */
    private static final String[] GENOMIC_SIGNAL_CHANNELS = {
        "sequences", "qualities", "tags"
    };
    /** Index columns signed by {@link #signGenomicRun}. L1
     *  (Task #82 Phase B.1, 2026-05-01): the M82-era VL-string
     *  {@code chromosomes} compound was replaced with
     *  {@code chromosome_ids} (uint16) +
     *  {@code chromosome_names} (compound) — both new columns
     *  are signed in place of the old single column. */
    private static final String[] GENOMIC_INDEX_COLUMNS = {
        "offsets", "lengths", "positions", "mapping_qualities", "flags",
        "chromosome_ids", "chromosome_names",
        // blocks_v1_grouped only (M103): the permutation back to input order.
        "input_index"
    };

    /** sign every signal channel and every genomic_index
     *  column under one {@code /study/genomic_runs/<name>/} group in
     *  one call, storing each signature on the dataset's
     *  {@code @ttio_signature} attribute.
     *
     *  <p>Returns a map from {@code "<sub>/<dataset>"} (e.g.
     *  {@code "signal_channels/sequences"},
     *  {@code "genomic_index/positions"}) to the prefixed signature
     *  string. Datasets that don't exist on disk are silently skipped
     *  (e.g. encrypted files have segments instead of plaintext signal
     *  channels).
     *
     *  <p>{@code algorithm} dispatches identically to
     *  {@link #sign(byte[], byte[], String)} —
     *  {@code "hmac-sha256"} (default) or {@code "ml-dsa-87"} (PQC).
     *
     *  <p><b>Cross-language equivalents:</b> Python
     *  {@code ttio.signatures.sign_genomic_run}, Objective-C
     *  {@code TTIOSignatureManager#signGenomicRun:}.
     *
     * M90.2
     */
    public static Map<String, String> signGenomicRun(
            StorageGroup runGroup, byte[] key, String algorithm) {
        Map<String, String> out = new LinkedHashMap<>();
        if (runGroup.hasChild("signal_channels")) {
            try (StorageGroup sig = runGroup.openGroup("signal_channels")) {
                for (String cname : GENOMIC_SIGNAL_CHANNELS) {
                    if (!sig.hasChild(cname)) continue;
                    for (String[] pd : channelDatasets(sig, cname)) {
                        StorageGroup parent = pd[0].isEmpty() ? sig : sig.openGroup(pd[0]);
                        try (StorageDataset ds = parent.openDataset(pd[1])) {
                            byte[] canonical = ds.readCanonicalBytes();
                            String s = sign(canonical, key, algorithm);
                            ds.setAttribute("ttio_signature", s);
                            out.put("signal_channels/" + (pd[0].isEmpty() ? "" : pd[0] + "/") + pd[1], s);
                        } finally {
                            if (parent != sig) parent.close();
                        }
                    }
                }
            }
        }
        if (runGroup.hasChild("genomic_index")) {
            try (StorageGroup idx = runGroup.openGroup("genomic_index")) {
                for (String cname : GENOMIC_INDEX_COLUMNS) {
                    if (!idx.hasChild(cname)) continue;
                    try (StorageDataset ds = idx.openDataset(cname)) {
                        byte[] canonical = ds.readCanonicalBytes();
                        String s = sign(canonical, key, algorithm);
                        ds.setAttribute("ttio_signature", s);
                        out.put("genomic_index/" + cname, s);
                    }
                }
            }
        }
        if (runGroup.hasChild("blocks")) {
            try (StorageGroup blocks = runGroup.openGroup("blocks")) {
                if (blocks.hasChild("index")) {
                    try (StorageDataset ds = blocks.openDataset("index")) {
                        String s = sign(ds.readCanonicalBytes(), key, algorithm);
                        ds.setAttribute("ttio_signature", s);
                        out.put("blocks/index", s);
                    }
                }
            }
        }
        return out;
    }

    /** The dataset(s) behind a genomic signal channel as
     *  {@code {parentGroupName, datasetName}} pairs ({@code ""} parent =
     *  the flat dataset; a group channel such as {@code sequences/data}
     *  under blocks_v1 or {@code sequences/refdiff_v2} lists its child
     *  datasets). */
    private static java.util.List<String[]> channelDatasets(StorageGroup sig, String cname) {
        java.util.List<String[]> out = new java.util.ArrayList<>();
        StorageGroup g = null;
        try {
            g = sig.openGroup(cname);
        } catch (RuntimeException e) {
            g = null;
        }
        if (g == null) {
            out.add(new String[]{ "", cname });
            return out;
        }
        try (StorageGroup grp = g) {
            for (String child : grp.childNames()) {
                boolean isGroup;
                try (StorageGroup sub = grp.openGroup(child)) { isGroup = true; }
                catch (RuntimeException e) { isGroup = false; }
                if (!isGroup) out.add(new String[]{ cname, child });
            }
        }
        return out;
    }

    /** {@code signGenomicRun} convenience overload defaulting to
     *  {@code "hmac-sha256"}. */
    public static Map<String, String> signGenomicRun(
            StorageGroup runGroup, byte[] key) {
        return signGenomicRun(runGroup, key, "hmac-sha256");
    }

    /** verify every signal channel and every genomic_index
     *  column under one genomic run. Returns {@code true} iff every
     *  present, signed dataset verifies under {@code key}.
     *
     *  <p>A dataset that was signed but is now tampered returns
     *  {@code false}. A dataset that exists but has no
     *  {@code @ttio_signature} attribute also returns {@code false}
     *  — that's intentional, since a partial-signature run is not a
     *  fully-signed run. Datasets that don't exist on disk are skipped.
     *
     * M90.2
     */
    public static boolean verifyGenomicRun(
            StorageGroup runGroup, byte[] key, String algorithm) {
        if (runGroup.hasChild("signal_channels")) {
            try (StorageGroup sig = runGroup.openGroup("signal_channels")) {
                for (String cname : GENOMIC_SIGNAL_CHANNELS) {
                    if (!sig.hasChild(cname)) continue;
                    for (String[] pd : channelDatasets(sig, cname)) {
                        StorageGroup parent = pd[0].isEmpty() ? sig : sig.openGroup(pd[0]);
                        try {
                            if (!verifyOneDataset(parent, pd[1], key, algorithm)) return false;
                        } finally {
                            if (parent != sig) parent.close();
                        }
                    }
                }
            }
        }
        if (runGroup.hasChild("genomic_index")) {
            try (StorageGroup idx = runGroup.openGroup("genomic_index")) {
                for (String cname : GENOMIC_INDEX_COLUMNS) {
                    if (!idx.hasChild(cname)) continue;
                    if (!verifyOneDataset(idx, cname, key, algorithm)) {
                        return false;
                    }
                }
            }
        }
        if (runGroup.hasChild("blocks")) {
            try (StorageGroup blocks = runGroup.openGroup("blocks")) {
                if (blocks.hasChild("index")
                        && !verifyOneDataset(blocks, "index", key, algorithm)) {
                    return false;
                }
            }
        }
        return true;
    }

    /** {@code verifyGenomicRun} convenience overload defaulting to
     *  {@code "hmac-sha256"}. */
    public static boolean verifyGenomicRun(
            StorageGroup runGroup, byte[] key) {
        return verifyGenomicRun(runGroup, key, "hmac-sha256");
    }

    /** Sign the byte channels of one assembly graph
     *  ({@code /study/assembly_graphs/<name>/}) — M98. The signed
     *  set is the {@code segments/sequences} channel; signatures
     *  cover it the way they cover genomic-run channels. The
     *  compound tables (records, links, paths, extras, line_index)
     *  are outside the signed set.
     *
     *  <p>Returns a map from {@code "segments/sequences"} to the
     *  prefixed signature stored on the dataset's
     *  {@code @ttio_signature} attribute. Absent datasets (a graph
     *  whose segments all carry {@code *}, or an encrypted file
     *  where the channel is a {@code sequences_segments} compound)
     *  are silently skipped.
     *
     *  <p><b>Cross-language equivalents:</b> Python
     *  {@code ttio.signatures.sign_assembly_graph}, Objective-C
     *  {@code TTIOSignatureManager#signAssemblyGraph:}.</p> */
    public static Map<String, String> signAssemblyGraph(
            StorageGroup graphGroup, byte[] key, String algorithm) {
        Map<String, String> out = new LinkedHashMap<>();
        if (graphGroup.hasChild("segments")) {
            try (StorageGroup seg = graphGroup.openGroup("segments")) {
                if (seg.hasChild("sequences")) {
                    try (StorageDataset ds = seg.openDataset("sequences")) {
                        byte[] canonical = ds.readCanonicalBytes();
                        String s = sign(canonical, key, algorithm);
                        ds.setAttribute("ttio_signature", s);
                        out.put("segments/sequences", s);
                    }
                }
            }
        }
        return out;
    }

    /** {@code signAssemblyGraph} convenience overload defaulting to
     *  {@code "hmac-sha256"}. */
    public static Map<String, String> signAssemblyGraph(
            StorageGroup graphGroup, byte[] key) {
        return signAssemblyGraph(graphGroup, key, "hmac-sha256");
    }

    /** Verify every dataset {@link #signAssemblyGraph} would sign.
     *  Returns {@code true} iff every present, signed dataset
     *  verifies under {@code key}; a present unsigned dataset fails,
     *  matching {@link #verifyGenomicRun}. Absent datasets are
     *  skipped. */
    public static boolean verifyAssemblyGraph(
            StorageGroup graphGroup, byte[] key, String algorithm) {
        if (graphGroup.hasChild("segments")) {
            try (StorageGroup seg = graphGroup.openGroup("segments")) {
                if (seg.hasChild("sequences")
                        && !verifyOneDataset(seg, "sequences", key,
                                             algorithm)) {
                    return false;
                }
            }
        }
        return true;
    }

    /** {@code verifyAssemblyGraph} convenience overload defaulting to
     *  {@code "hmac-sha256"}. */
    public static boolean verifyAssemblyGraph(
            StorageGroup graphGroup, byte[] key) {
        return verifyAssemblyGraph(graphGroup, key, "hmac-sha256");
    }

    private static boolean verifyOneDataset(StorageGroup parent, String name,
                                              byte[] key, String algorithm) {
        try (StorageDataset ds = parent.openDataset(name)) {
            if (!ds.hasAttribute("ttio_signature")) return false;
            Object sigObj = ds.getAttribute("ttio_signature");
            if (sigObj == null) return false;
            String stored;
            if (sigObj instanceof byte[] b) {
                stored = new String(b, StandardCharsets.UTF_8);
            } else {
                stored = sigObj.toString();
            }
            byte[] canonical = ds.readCanonicalBytes();
            return verify(canonical, stored, key, algorithm);
        }
    }
}
