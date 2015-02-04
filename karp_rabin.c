#include "karp_rabin.h"
#include <gmp.h>
#include <stdio.h>
#include <assert.h>

int main(void) {
    int n = 100, m = 20;
    fingerprinter printer = fingerprinter_build(n, 0);
    gmp_printf("p = %Zd\n", printer->p);
    gmp_printf("r = %Zd\n", printer->r);

    fingerprint print = init_fingerprint();
    set_fingerprint(printer, "aaaaabbbbbcccccaaaaa", m, print);

    gmp_printf("uv finger = %Zd\n", print->finger);
    gmp_printf("uv r_k = %Zd\n", print->r_k);
    gmp_printf("uv r_mk = %Zd\n", print->r_mk);

    fingerprint prefix = init_fingerprint();
    set_fingerprint(printer, "aaaaa", 5, prefix);

    fingerprint v = init_fingerprint();
    fingerprint_suffix(printer, print, prefix, v);

    fingerprint suffix = init_fingerprint();
    set_fingerprint(printer, "bbbbbcccccaaaaa", 15, suffix);
    assert(fingerprint_equals(v, suffix));

    fingerprint u = init_fingerprint();
    fingerprint_prefix(printer, print, suffix, u);
    assert(fingerprint_equals(u, prefix));

    fingerprint uv = init_fingerprint();
    fingerprint_concat(printer, prefix, suffix, uv);
    assert(fingerprint_equals(uv, print));

    fingerprint empty = init_fingerprint();
    fingerprint_suffix(printer, print, empty, v);
    assert(fingerprint_equals(v, print));

    fingerprint_prefix(printer, print, empty, v);
    assert(fingerprint_equals(v, print));

    fingerprint_concat(printer, print, empty, v);
    assert(fingerprint_equals(v, print));

    fingerprint_concat(printer, empty, print, v);
    assert(fingerprint_equals(v, print));

    fingerprint_free(print);
    fingerprint_free(prefix);
    fingerprint_free(suffix);
    fingerprint_free(u);
    fingerprint_free(v);
    fingerprint_free(uv);
    fingerprint_free(empty);

    return 0;
}