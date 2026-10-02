/* SPDX-License-Identifier: Apache-2.0 */
/*
 * Builds libobjc2's weak-reference table before main() so that it is
 * torn down after GNUstep's exit clean-up, not before it.
 *
 * libobjc2 (v2.3 and master) keeps the table in a function-local C++
 * static, weakRefs()::w, so its destructor is registered with
 * __cxa_atexit the first time anything stores a weak reference. In a
 * TTI-O tool that first store happens inside main() (NSNotificationCenter
 * holds its observers weakly), after __libc_start_main has registered
 * _dl_fini. gnustep-base registers its own clean-up, handleExit, while
 * the libraries initialise, before _dl_fini. exit() runs handlers in
 * reverse order, so it frees the table first and only then reaches
 * _dl_fini, which runs handleExit from gnustep-base's __cxa_finalize.
 * +[NSUserDefaults atExit] then releases the shared defaults, an
 * observer with weak references, and objc_delete_weak_refs probes the
 * freed buckets. Valgrind reports that read on every exit; it is a
 * SIGSEGV only when the freed memory has been reused or returned to the
 * kernel, which is the intermittent CI crash in TtioTransportEncode.
 *
 * Storing one weak reference from a library constructor builds the
 * table before _dl_fini is registered. Its destructor then runs from
 * libobjc's own __cxa_finalize, and _dl_fini finalises libobjc after
 * gnustep-base, which depends on it, so handleExit sees a live table.
 * Every tool and the test runner link libTTIO, so all of them get it.
 */
#import <Foundation/Foundation.h>

/* A file-scope weak variable, so the store is an objc_storeWeak call
 * the ARC optimiser cannot drop the way it drops an unused local. */
static __weak id TTIOWeakTableProbe;

__attribute__((constructor))
static void TTIOConstructWeakRefTable(void)
{
    NSObject *probe = [NSObject new];
    TTIOWeakTableProbe = probe;
}
