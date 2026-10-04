cmake_minimum_required(VERSION 3.16)
project(ttio_rans C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(TTIO_RANS_ENABLE_TSAN
    "Build a separate test_thread_safety_tsan target with -fsanitize=thread" OFF)
option(TTIO_RANS_BUILD_JNI
    "Build the JNI wrapper library (libttio_rans_jni.so) for Java integration" OFF)
# BUILD_TESTING gates the test executables. Defaults ON so local dev sees
# the full suite. The release matrix and JNI-only builds set OFF because
# some tests (e.g. test_fqzcomp_qual_threaded) use POSIX setenv/unsetenv
# unavailable on MinGW UCRT.
option(BUILD_TESTING "Build native test executables" ON)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# Optional gcov/llvm-cov instrumentation (R8 coverage visibility).
# Default OFF so production / wheel / ObjC / Java builds are unaffected.
# Placed before any target so it instruments the library AND every test
# executable defined below in one shot. MSVC has no --coverage; gcov is a
# GCC/Clang feature, so this is a no-op there. Use a fresh build dir
# (native/_covbuild) — never the shared native/_build.
option(TTIO_COVERAGE "Instrument with --coverage for gcov/llvm-cov reporting" OFF)
if(TTIO_COVERAGE AND NOT MSVC)
    add_compile_options(--coverage -O0 -g)
    add_link_options(--coverage)
    message(STATUS "ttio_rans: gcov coverage instrumentation ENABLED")
elseif(TTIO_COVERAGE AND MSVC)
    message(WARNING "TTIO_COVERAGE=ON has no effect on MSVC (no --coverage support); ignoring")
endif()
if(TTIO_COVERAGE AND TTIO_RANS_ENABLE_TSAN)
    message(WARNING "TTIO_COVERAGE and TTIO_RANS_ENABLE_TSAN are both ON; gcov + ThreadSanitizer instrumentation together produce unreliable coverage and race data")
endif()

add_library(ttio_rans SHARED
    src/rans_core.c
    src/rans_encode_scalar.c
    src/rans_decode_scalar.c
    src/rans_decode_streaming.c
    src/rans_decode_m94z.c
    src/rans_encode_adaptive.c
    src/rans_decode_adaptive.c
    src/rc_cram.c
    src/fqzcomp_qual.c
    src/fqzcomp_seqctx.c
    src/m94z_v6.c
    src/ttio_engine.c
    src/m94z_qual.c
    src/m94z_v4_wire.c
    src/mate_info_v2.c
    src/ref_diff_v2.c
    src/name_tok_v2.c
    src/sam_tags.c
    src/seq_cm.c
    src/seq_group.c
    src/dispatch.c
    src/rans_o0.c
    src/threadpool.c
    src/wire_format.c
)
target_include_directories(ttio_rans
    PUBLIC  include
    PRIVATE src
)
# Compiler warnings: GCC/Clang spellings on Unix; MSVC's /W3 on Windows.
# MSVC chokes on -Wextra (parses 'e' as a numeric option) and rejects
# -Wpedantic outright; both flags are GCC/Clang-only.
if(MSVC)
    target_compile_options(ttio_rans PRIVATE /W3)
else()
    target_compile_options(ttio_rans PRIVATE -Wall -Wextra -Wpedantic)
endif()
# libm is implicit on MSVC (math.h symbols live in the C runtime).
if(MSVC)
    target_link_libraries(ttio_rans PRIVATE Threads::Threads ZLIB::ZLIB)
else()
    target_link_libraries(ttio_rans PRIVATE Threads::Threads m ZLIB::ZLIB)
endif()

# x86_64 detection. CMAKE_SYSTEM_PROCESSOR is empty under some
# toolchains (MSYS2 ucrt64 among them), which silently dropped the
# SIMD sources and left dispatch.c referencing symbols that were never
# built. Ask the compiler instead, and keep the old match as well so
# cross builds that set the processor explicitly still work.
include(CheckCSourceCompiles)
check_c_source_compiles("
#if !defined(__x86_64__) && !defined(_M_X64)
#error not x86_64
#endif
int main(void) { return 0; }
" TTIO_TARGET_IS_X86_64)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64" OR TTIO_TARGET_IS_X86_64)
    set(TTIO_X86_64 TRUE)
else()
    set(TTIO_X86_64 FALSE)
endif()

# x86_64: build SSE4.1 / AVX2 kernels.
#
# GCC/Clang need explicit per-file -msse4.1 / -mavx2 flags to enable
# the wider intrinsics in those compilation units (without them, the
# intrinsics-bearing headers refuse to expose SSE4.1+ symbols).
#
# MSVC has no per-file SSE4.1 flag (intrinsics in <smmintrin.h> are
# always available for x64 targets) and uses /arch:AVX2 at the
# whole-compilation-unit level for AVX2. For per-file AVX2 enablement
# under MSVC we'd need to split the AVX2 sources into a separate
# library target — accepting the simpler "no per-file flag, intrinsics
# work via <immintrin.h>" path here gets MSVC compiling; an MSVC-side
# AVX2 dispatch optimisation is a separate enhancement.
if(TTIO_X86_64)
    target_sources(ttio_rans PRIVATE
        src/rans_encode_sse41.c
        src/rans_decode_sse41.c
        src/rans_encode_avx2.c
        src/rans_decode_avx2.c
    )
    if(NOT MSVC)
        set_source_files_properties(
            src/rans_encode_sse41.c src/rans_decode_sse41.c
            PROPERTIES COMPILE_FLAGS "-msse4.1"
        )
        set_source_files_properties(
            src/rans_encode_avx2.c src/rans_decode_avx2.c
            PROPERTIES COMPILE_FLAGS "-mavx2"
        )
    endif()
endif()

# Install. Defaults ON for standalone / ObjC / Java consumers. The Python wheel
# build sets TTIO_RANS_INSTALL=OFF (it installs the library into ttio/.libs via
# python/CMakeLists.txt) so these rules don't drop a stray lib/ + include/ into
# the wheel root.
option(TTIO_RANS_INSTALL "Install libttio_rans and its public header" ON)
if(TTIO_RANS_INSTALL)
    install(TARGETS ttio_rans LIBRARY DESTINATION lib)
    install(FILES include/ttio_rans.h DESTINATION include)
endif()

# Optional JNI wrapper.  Enable with `-DTTIO_RANS_BUILD_JNI=ON`.
# Java loads the library by base name "ttio_rans_jni" via System.loadLibrary,
# which expands to libttio_rans_jni.so / libttio_rans_jni.dylib /
# ttio_rans_jni.dll on Linux/macOS/Windows respectively.
if(TTIO_RANS_BUILD_JNI)
    # We only need the jni.h headers — not AWT, which is missing on
    # headless JDK installs. Use Java::JNI from FindJavaInclude when
    # available; otherwise fall back to manual JAVA_HOME probing.
    find_package(JNI COMPONENTS JVM)
    if(NOT JNI_INCLUDE_DIRS)
        # Manual fallback: probe $JAVA_HOME/include and the standard
        # Debian/Ubuntu JDK location.
        if(DEFINED ENV{JAVA_HOME})
            set(_jdk_root "$ENV{JAVA_HOME}")
        else()
            execute_process(
                COMMAND bash -c "readlink -f $(which javac) | sed 's|/bin/javac||'"
                OUTPUT_VARIABLE _jdk_root
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        endif()
        if(_jdk_root AND EXISTS "${_jdk_root}/include/jni.h")
            set(JNI_INCLUDE_DIRS "${_jdk_root}/include" "${_jdk_root}/include/linux")
            message(STATUS "JNI: using ${_jdk_root}/include (manual probe)")
        else()
            message(FATAL_ERROR "JNI headers not found; set JAVA_HOME to a JDK")
        endif()
    endif()
    add_library(ttio_rans_jni SHARED src/ttio_rans_jni.c)
    target_link_libraries(ttio_rans_jni PRIVATE ttio_rans)
    target_include_directories(ttio_rans_jni PRIVATE ${JNI_INCLUDE_DIRS})
    if(MSVC)
        target_compile_options(ttio_rans_jni PRIVATE /W3)
    else()
        target_compile_options(ttio_rans_jni PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    # Library naming follows each platform's native convention:
    #   Linux:   libttio_rans_jni.so    (lib prefix kept; CMake default)
    #   macOS:   libttio_rans_jni.dylib (lib prefix kept; CMake default)
    #   Windows: ttio_rans_jni.dll      (no lib prefix — Windows convention)
    #
    # MinGW-w64 (used on Windows) defaults to "lib" prefix even though
    # Windows convention is no prefix; explicitly set PREFIX="" on
    # WIN32 to override. tio-browser's NativeLibraryLoader uses
    # System.load on an absolute extracted path so either form would
    # work there; this preserves the convention the existing comment
    # block describes.
    if(WIN32)
        set_target_properties(ttio_rans_jni PROPERTIES PREFIX "")
    else()
        set_target_properties(ttio_rans_jni PROPERTIES PREFIX "lib")
    endif()
    install(TARGETS ttio_rans_jni LIBRARY DESTINATION lib)
endif()

# Tests
if(BUILD_TESTING)
enable_testing()
add_executable(test_roundtrip tests/test_roundtrip.c)
target_link_libraries(test_roundtrip ttio_rans)
target_compile_options(test_roundtrip PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME roundtrip COMMAND test_roundtrip)

add_executable(test_fqzcomp_seqctx tests/test_fqzcomp_seqctx.c)
target_link_libraries(test_fqzcomp_seqctx ttio_rans)
target_compile_options(test_fqzcomp_seqctx PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME fqzcomp_seqctx COMMAND test_fqzcomp_seqctx)

add_executable(test_engine_selection tests/test_engine_selection.c)
target_link_libraries(test_engine_selection PRIVATE ttio_rans Threads::Threads)
target_compile_options(test_engine_selection PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME engine_selection COMMAND test_engine_selection)

add_executable(test_engine_api tests/test_engine_api.c)
target_link_libraries(test_engine_api PRIVATE ttio_rans Threads::Threads)
target_compile_options(test_engine_api PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME engine_api COMMAND test_engine_api)

add_executable(test_v6_model tests/test_v6_model.c)
target_link_libraries(test_v6_model ttio_rans)
target_compile_options(test_v6_model PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME v6_model COMMAND test_v6_model)

add_executable(test_m94z_v6_chain tests/test_m94z_v6_chain.c)
target_link_libraries(test_m94z_v6_chain ttio_rans)
target_compile_options(test_m94z_v6_chain PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME m94z_v6_chain COMMAND test_m94z_v6_chain)

add_executable(test_m94z_v6_block tests/test_m94z_v6_block.c)
target_link_libraries(test_m94z_v6_block PRIVATE ttio_rans Threads::Threads)
target_compile_options(test_m94z_v6_block PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME m94z_v6_block COMMAND test_m94z_v6_block)

add_executable(test_m94z_v6_seqctx tests/test_m94z_v6_seqctx.c)
target_link_libraries(test_m94z_v6_seqctx PRIVATE ttio_rans Threads::Threads)
target_compile_options(test_m94z_v6_seqctx PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME m94z_v6_seqctx COMMAND test_m94z_v6_seqctx)

add_executable(test_m94z_qual_umbrella tests/test_m94z_qual_umbrella.c)
target_link_libraries(test_m94z_qual_umbrella ttio_rans)
target_compile_options(test_m94z_qual_umbrella PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME m94z_qual_umbrella COMMAND test_m94z_qual_umbrella)

add_executable(test_thread_safety tests/test_thread_safety.c)
target_link_libraries(test_thread_safety PRIVATE ttio_rans Threads::Threads)
target_compile_options(test_thread_safety PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME thread_safety COMMAND test_thread_safety)

add_executable(test_v2_format tests/test_v2_format.c)
target_link_libraries(test_v2_format PRIVATE ttio_rans Threads::Threads)
target_compile_options(test_v2_format PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME v2_format COMMAND test_v2_format)

add_executable(test_streaming tests/test_streaming.c)
target_link_libraries(test_streaming PRIVATE ttio_rans)
target_compile_options(test_streaming PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME streaming COMMAND test_streaming)

add_executable(test_m94z_decode tests/test_m94z_decode.c)
target_link_libraries(test_m94z_decode PRIVATE ttio_rans)
target_compile_options(test_m94z_decode PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME m94z_decode COMMAND test_m94z_decode)

# L2: adaptive M94.Z (CRAM-mimic) tests.
add_executable(test_adaptive_roundtrip tests/test_adaptive_roundtrip.c)
target_link_libraries(test_adaptive_roundtrip PRIVATE ttio_rans)
target_compile_options(test_adaptive_roundtrip PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME adaptive_roundtrip COMMAND test_adaptive_roundtrip)

add_executable(test_adaptive_byte_pairing tests/test_adaptive_byte_pairing.c)
target_link_libraries(test_adaptive_byte_pairing PRIVATE ttio_rans)
target_compile_options(test_adaptive_byte_pairing PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME adaptive_byte_pairing COMMAND test_adaptive_byte_pairing)

add_executable(test_adaptive_halve_boundary tests/test_adaptive_halve_boundary.c)
target_link_libraries(test_adaptive_halve_boundary PRIVATE ttio_rans)
target_compile_options(test_adaptive_halve_boundary PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME adaptive_halve_boundary COMMAND test_adaptive_halve_boundary)

add_executable(test_adaptive_max_sym_bounds tests/test_adaptive_max_sym_bounds.c)
target_link_libraries(test_adaptive_max_sym_bounds PRIVATE ttio_rans)
target_compile_options(test_adaptive_max_sym_bounds PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME adaptive_max_sym_bounds COMMAND test_adaptive_max_sym_bounds)
# Phase 1 byte-equality gate: rc_cram vs htscodecs reference.
# Run manually: ./test_rc_cram_byte_equal /tmp/rc_cram_ours.bin
# then cmp /tmp/rc_cram_ours.bin /tmp/rc_cram_htscodecs.bin
add_executable(test_rc_cram_byte_equal tests/test_rc_cram_byte_equal.c)
target_include_directories(test_rc_cram_byte_equal PRIVATE src)
target_link_libraries(test_rc_cram_byte_equal ttio_rans)
target_compile_options(test_rc_cram_byte_equal PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME rc_cram_byte_equal COMMAND test_rc_cram_byte_equal)

# Phase 2 byte-equality + round-trip gate: ttio_fqzcomp_qual_compress
# strategy 1 vs htscodecs --strategy=1 on chr22.
# Inputs: extract via tools/perf/m94z_v4_prototype/extract_chr22_inputs.py
# Reference: tools/perf/m94z_v4_prototype/fqzcomp_htscodecs_ref
# Run manually: ./test_fqzcomp_qual_strategy1 /tmp/our_chr22_strategy1.fqz
# then cmp /tmp/our_chr22_strategy1.fqz /tmp/htscodecs_chr22_strategy1.fqz
add_executable(test_fqzcomp_qual_strategy1 tests/test_fqzcomp_qual_strategy1.c)
target_include_directories(test_fqzcomp_qual_strategy1 PRIVATE src)
target_link_libraries(test_fqzcomp_qual_strategy1 ttio_rans)
target_compile_options(test_fqzcomp_qual_strategy1 PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME fqzcomp_qual_strategy1 COMMAND test_fqzcomp_qual_strategy1)

# Phase 3 byte-equality gate: auto-tune (strategy_hint = -1) vs htscodecs
# fqz_compress(gp=NULL) across multiple corpora.
# Argv-parametrized; driven by tools/perf/htscodecs_compare.sh which
# extracts each corpus's inputs and compares bytes against
# tools/perf/m94z_v4_prototype/fqzcomp_htscodecs_ref_autotune.
# Not registered as a ctest (requires per-invocation argv + corpus on disk).
add_executable(test_fqzcomp_qual_autotune tests/test_fqzcomp_qual_autotune.c)
target_include_directories(test_fqzcomp_qual_autotune PRIVATE src)
target_link_libraries(test_fqzcomp_qual_autotune ttio_rans)
target_compile_options(test_fqzcomp_qual_autotune PRIVATE -Wall -Wextra -Wpedantic)

# Stage 4 fast unit gate: threaded auto-tune histogram pass produces
# byte-identical output to the sequential path on synthetic input.
# Cross-corpus integration is gated by python/tests/integration/
# test_m94z_v4_byte_exact.py (and Java/ObjC mirrors); this is the
# fast-running ctest that catches partition / merge regressions.
add_executable(test_fqzcomp_qual_threaded tests/test_fqzcomp_qual_threaded.c)
target_include_directories(test_fqzcomp_qual_threaded PRIVATE src)
target_link_libraries(test_fqzcomp_qual_threaded PRIVATE ttio_rans)
target_compile_options(test_fqzcomp_qual_threaded PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME fqzcomp_qual_threaded COMMAND test_fqzcomp_qual_threaded)

# #11 mate_info v2: varint + zigzag helper unit tests.
add_executable(test_mate_info_v2_helpers tests/test_mate_info_v2_helpers.c)
target_link_libraries(test_mate_info_v2_helpers PRIVATE ttio_rans)
target_include_directories(test_mate_info_v2_helpers PRIVATE src)
target_compile_options(test_mate_info_v2_helpers PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME mate_info_v2_helpers COMMAND test_mate_info_v2_helpers)


add_executable(test_rans_o0_byte_exact tests/test_rans_o0_byte_exact.c)
target_link_libraries(test_rans_o0_byte_exact PRIVATE ttio_rans)
target_compile_options(test_rans_o0_byte_exact PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME rans_o0_byte_exact COMMAND test_rans_o0_byte_exact)

# #11 mate_info v2 task 3: invariants test (MF exhaustiveness, round-trip, error paths).
add_executable(test_mate_info_v2_invariants tests/test_mate_info_v2_invariants.c)
target_link_libraries(test_mate_info_v2_invariants PRIVATE ttio_rans)
target_include_directories(test_mate_info_v2_invariants PRIVATE src)
add_test(NAME mate_info_v2_invariants COMMAND test_mate_info_v2_invariants)

# #11 mate_info v2 task 4: parameterized stress test + I4 + I5 invariants.
add_executable(test_mate_info_v2_stress tests/test_mate_info_v2_stress.c)
target_link_libraries(test_mate_info_v2_stress PRIVATE ttio_rans)
target_include_directories(test_mate_info_v2_stress PRIVATE src)
add_test(NAME mate_info_v2_stress COMMAND test_mate_info_v2_stress)

# #11 ref_diff v2 task 2: ACGT 2-bit + cigar count helpers + tests.
add_executable(test_ref_diff_v2_helpers tests/test_ref_diff_v2_helpers.c)
target_link_libraries(test_ref_diff_v2_helpers PRIVATE ttio_rans)
target_include_directories(test_ref_diff_v2_helpers PRIVATE src)
add_test(NAME ref_diff_v2_helpers COMMAND test_ref_diff_v2_helpers)

# #11 ref_diff v2 task 3: full encoder + decoder + invariants I3/I5/I6.
add_executable(test_ref_diff_v2_invariants tests/test_ref_diff_v2_invariants.c)
target_link_libraries(test_ref_diff_v2_invariants PRIVATE ttio_rans)
target_include_directories(test_ref_diff_v2_invariants PRIVATE src)
add_test(NAME ref_diff_v2_invariants COMMAND test_ref_diff_v2_invariants)

# #11 ref_diff v2 task 4: parameterized stress (3 sizes x 6 patterns) + I2 invariant.
add_executable(test_ref_diff_v2_stress tests/test_ref_diff_v2_stress.c)
target_link_libraries(test_ref_diff_v2_stress PRIVATE ttio_rans)
target_include_directories(test_ref_diff_v2_stress PRIVATE src)
add_test(NAME ref_diff_v2_stress COMMAND test_ref_diff_v2_stress)

# #11 ch3 task 2: name_tok_v2 helpers + tests.
add_executable(test_name_tok_v2_helpers tests/test_name_tok_v2_helpers.c)
target_link_libraries(test_name_tok_v2_helpers PRIVATE ttio_rans)
target_include_directories(test_name_tok_v2_helpers PRIVATE src)
add_test(NAME name_tok_v2_helpers COMMAND test_name_tok_v2_helpers)

# #11 ch3 task 3: name_tok_v2 encoder + decoder + invariants.
add_executable(test_name_tok_v2_invariants tests/test_name_tok_v2_invariants.c)
target_link_libraries(test_name_tok_v2_invariants PRIVATE ttio_rans)
target_include_directories(test_name_tok_v2_invariants PRIVATE src)
add_test(NAME name_tok_v2_invariants COMMAND test_name_tok_v2_invariants)

# #11 ch3 task 4: stress + malformed-input tests.
add_executable(test_name_tok_v2_stress tests/test_name_tok_v2_stress.c)
target_link_libraries(test_name_tok_v2_stress PRIVATE ttio_rans)
target_include_directories(test_name_tok_v2_stress PRIVATE src)
add_test(NAME name_tok_v2_stress COMMAND test_name_tok_v2_stress)

# M101: SAM_TAGS codec (codec id 18).
add_executable(test_sam_tags tests/test_sam_tags.c)
target_link_libraries(test_sam_tags PRIVATE ttio_rans)
target_include_directories(test_sam_tags PRIVATE src)
add_test(NAME sam_tags COMMAND test_sam_tags)

add_executable(test_seq_cm tests/test_seq_cm.c)
target_link_libraries(test_seq_cm PRIVATE ttio_rans)
target_include_directories(test_seq_cm PRIVATE src)
add_test(NAME seq_cm COMMAND test_seq_cm)

add_executable(test_seq_group tests/test_seq_group.c)
target_link_libraries(test_seq_group PRIVATE ttio_rans)
add_test(NAME seq_group COMMAND test_seq_group)
endif()  # BUILD_TESTING


# Optional ThreadSanitizer build of the same test, gated by an option so
# the default ctest run remains fast and portable.  Enable with:
#   cmake -DTTIO_RANS_ENABLE_TSAN=ON ..
if(TTIO_RANS_ENABLE_TSAN)
    add_library(ttio_rans_tsan SHARED
        src/rans_core.c
        src/rans_encode_scalar.c
        src/rans_decode_scalar.c
        src/rans_decode_streaming.c
        src/dispatch.c
        src/threadpool.c
        src/wire_format.c
    )
    target_include_directories(ttio_rans_tsan
        PUBLIC  include
        PRIVATE src
    )
    target_compile_options(ttio_rans_tsan PRIVATE
        -Wall -Wextra -Wpedantic -g -O1 -fsanitize=thread)
    target_link_options(ttio_rans_tsan PRIVATE -fsanitize=thread)
    target_link_libraries(ttio_rans_tsan PRIVATE Threads::Threads)
    if(TTIO_X86_64)
        target_sources(ttio_rans_tsan PRIVATE
            src/rans_encode_sse41.c
            src/rans_decode_sse41.c
            src/rans_encode_avx2.c
            src/rans_decode_avx2.c
        )
        set_source_files_properties(
            src/rans_encode_sse41.c src/rans_decode_sse41.c
            PROPERTIES COMPILE_FLAGS "-msse4.1"
        )
        set_source_files_properties(
            src/rans_encode_avx2.c src/rans_decode_avx2.c
            PROPERTIES COMPILE_FLAGS "-mavx2"
        )
    endif()

    add_executable(test_thread_safety_tsan tests/test_thread_safety.c)
    target_compile_options(test_thread_safety_tsan PRIVATE
        -Wall -Wextra -Wpedantic -g -O1 -fsanitize=thread)
    target_link_options(test_thread_safety_tsan PRIVATE -fsanitize=thread)
    target_link_libraries(test_thread_safety_tsan
        PRIVATE ttio_rans_tsan Threads::Threads)
    add_test(NAME thread_safety_tsan COMMAND test_thread_safety_tsan)
endif()

# Benchmark (not registered as a ctest; run manually). Skipped under
# coverage: it is never exercised by ctest, and its -O3 would clash with
# the directory-level --coverage -O0 instrumentation.
if(NOT TTIO_COVERAGE)
    add_executable(bench_throughput bench/bench_throughput.c)
    target_link_libraries(bench_throughput PRIVATE ttio_rans)
    target_compile_options(bench_throughput PRIVATE -O3)

    add_executable(v6_acceptance tools/v6_acceptance.c)
    target_link_libraries(v6_acceptance PRIVATE ttio_rans)
    target_compile_options(v6_acceptance PRIVATE -O2)

    add_executable(bench_m94z_qual bench/bench_m94z_qual.c)
    target_link_libraries(bench_m94z_qual PRIVATE ttio_rans Threads::Threads)
    target_compile_options(bench_m94z_qual PRIVATE -O3)

    add_executable(bench_v6_ratio bench/bench_v6_ratio.c)
    target_link_libraries(bench_v6_ratio PRIVATE ttio_rans Threads::Threads)
    target_compile_options(bench_v6_ratio PRIVATE -O3)
endif()
