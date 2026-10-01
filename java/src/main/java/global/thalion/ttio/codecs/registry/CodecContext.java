package global.thalion.ttio.codecs.registry;

import global.thalion.ttio.codecs.ReferenceResolver;
import java.util.function.Supplier;

/** Run-derived context for codecs. All fields nullable; plain codecs ignore it.
 *  Built once per GenomicRun (decode) or per channel (encode). */
public record CodecContext(
        int[] readLengths,
        int[] revcompFlags,
        Integer elementSize,
        Integer readCount,
        long[] positions,
        Supplier<String[]> cigarsProvider,
        Long totalBases,
        String[] chromosomes,
        short[] ownChromIds,
        long[] ownPositions,
        Integer nRecords,
        ReferenceResolver referenceResolver,
        long[] offsets,
        byte[] reference,
        byte[] referenceMd5,
        String referenceUri,
        Integer readsPerSlice,
        Long sliceBytes,
        byte[] sequences,
        Supplier<byte[]> sequencesProvider,
        Integer qualStrategyHint,
        /** SAM_TAGS (M101) encode: the reference bases per
         *  {@code ownChromIds} value, for MD/NM derivation; null or
         *  all-null = no derivation. */
        java.util.List<byte[]> tagReferences,
        /** SAM_TAGS (M101) decode: the same list, built lazily and
         *  called once per decode. */
        Supplier<java.util.List<byte[]>> tagReferencesProvider) {

    public static CodecContext empty() { return builder().build(); }

    public static Builder builder() { return new Builder(); }

    public static final class Builder {
        private int[] readLengths;
        private int[] revcompFlags;
        private Integer elementSize;
        private Integer readCount;
        private long[] positions;
        private Supplier<String[]> cigarsProvider;
        private Long totalBases;
        private String[] chromosomes;
        private short[] ownChromIds;
        private long[] ownPositions;
        private Integer nRecords;
        private ReferenceResolver referenceResolver;
        private long[] offsets;
        private byte[] reference;
        private byte[] referenceMd5;
        private String referenceUri;
        private Integer readsPerSlice;
        private Long sliceBytes;
        private byte[] sequences;
        private Supplier<byte[]> sequencesProvider;
        private Integer qualStrategyHint;
        private java.util.List<byte[]> tagReferences;
        private Supplier<java.util.List<byte[]>> tagReferencesProvider;

        public Builder readLengths(int[] v) { this.readLengths = v; return this; }
        public Builder revcompFlags(int[] v) { this.revcompFlags = v; return this; }
        public Builder elementSize(Integer v) { this.elementSize = v; return this; }
        public Builder readCount(Integer v) { this.readCount = v; return this; }
        public Builder positions(long[] v) { this.positions = v; return this; }
        public Builder cigarsProvider(Supplier<String[]> v) { this.cigarsProvider = v; return this; }
        public Builder totalBases(Long v) { this.totalBases = v; return this; }
        public Builder chromosomes(String[] v) { this.chromosomes = v; return this; }
        public Builder ownChromIds(short[] v) { this.ownChromIds = v; return this; }
        public Builder ownPositions(long[] v) { this.ownPositions = v; return this; }
        public Builder nRecords(Integer v) { this.nRecords = v; return this; }
        public Builder referenceResolver(ReferenceResolver v) { this.referenceResolver = v; return this; }
        public Builder offsets(long[] v) { this.offsets = v; return this; }
        public Builder reference(byte[] v) { this.reference = v; return this; }
        public Builder referenceMd5(byte[] v) { this.referenceMd5 = v; return this; }
        public Builder referenceUri(String v) { this.referenceUri = v; return this; }
        public Builder readsPerSlice(Integer v) { this.readsPerSlice = v; return this; }
        /** M97: REF_DIFF_V2 slice byte budget; null/0 = the fixed-count rule. */
        public Builder sliceBytes(Long v) { this.sliceBytes = v; return this; }
        public Builder sequences(byte[] v) { this.sequences = v; return this; }
        public Builder sequencesProvider(Supplier<byte[]> v) { this.sequencesProvider = v; return this; }
        public Builder qualStrategyHint(Integer v) { this.qualStrategyHint = v; return this; }
        public Builder tagReferences(java.util.List<byte[]> v) { this.tagReferences = v; return this; }
        public Builder tagReferencesProvider(Supplier<java.util.List<byte[]>> v) {
            this.tagReferencesProvider = v; return this;
        }

        public CodecContext build() {
            return new CodecContext(readLengths, revcompFlags, elementSize, readCount,
                positions, cigarsProvider, totalBases, chromosomes, ownChromIds,
                ownPositions, nRecords, referenceResolver, offsets, reference,
                referenceMd5, referenceUri, readsPerSlice, sliceBytes,
                sequences, sequencesProvider, qualStrategyHint,
                tagReferences, tagReferencesProvider);
        }
    }
}
