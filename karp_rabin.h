/*
    karp_rabin.h
    Library for Karp-Rabin fingerprints.
    Utilises the GNU Multile Precision Arithmetic library (https://gmplib.org/) and dev/urandom.
*/

#ifndef KARP_RABIN
#define KARP_RABIN

#include <gmp.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>

/*
    mpz_equals
    Small function to check if two MP-Integers are equal.
    Parameters:
        mpz_t x - First number
        mpz_t y - Second number
    Returns int:
        1 if x == y
        0 otherwise
*/
inline int mpz_equals(mpz_t x, mpz_t y) {
    return ((x->_mp_size == y->_mp_size) && mpn_cmp(x->_mp_d, y->_mp_d, x->_mp_size) == 0);
}

/*
    compare
    Small function to compare two MP-Integers.
    Parameters:
        mpz_t x - First number
        mpz_t y - Second number
    Returns int:
        1 if x > y
        -1 if x < y
        0 otherwise
*/
inline int compare(mpz_t x, mpz_t y) {
    if (x->_mp_size > y->_mp_size) return 1;
    else if (y->_mp_size > x->_mp_size) return -1;
    else return mpn_cmp(x->_mp_d, y->_mp_d, x->_mp_size);
}

/*
    typedef struct fingerprinter_t *fingerprinter
    Structure to hold numbers for computing fingerprints.
    Components:
        mpz_t p - Prime number
        mpz_t r - Random number such that 0 <= r < p
*/
typedef struct fingerprinter_t {
    mpz_t p, r;
} *fingerprinter;

/*
    fingerprinter_build
    Constructs a fingerprint for a problem size and accuracy.
    Parameters:
        unsigned int n     - Size of the text
        unsigned int alpha - Desired accuracy
    Returns fingerprinter:
        The constructed fingerprint
    Notes:
        Primality is tested using a probabilistic algorithm. For practical purposes it is adequate.
        Chances of a collision are at most 1/n^(1+alpha).
*/
fingerprinter fingerprinter_build(unsigned int n, unsigned int alpha) {
    fingerprinter printer = malloc(sizeof(struct fingerprinter_t));

    mpz_init_set_ui(printer->p, n);
    mpz_pow_ui(printer->p, printer->p, 2 + alpha);
    mpz_nextprime(printer->p, printer->p);

    gmp_randstate_t state;
    unsigned long seed;
    size_t seed_len = 0;
    int f = open("/dev/urandom", O_RDONLY);
    while (seed_len < sizeof seed) {
        size_t result = read(f, ((char*)&seed) + seed_len, (sizeof seed) - seed_len);
        seed_len += result;
    }
    close(f);
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, seed);

    mpz_init(printer->r);
    mpz_sub_ui(printer->r, printer->p, 1);
    mpz_urandomm(printer->r, state, printer->r);
    mpz_add_ui(printer->r, printer->r, 1);

    return printer;
}

/*
    fingerprinter_free
    Frees a fingerprinter from memory.
    Parameters:
        fingerprinter printer - The fingerprinter to free
*/
void fingerprinter_free(fingerprinter printer) {
    mpz_clear(printer->p);
    mpz_clear(printer->r);
    free(printer);
}

/*
    typedef struct fingerprint_t *fingerprint
    Structure to hold fingerprints.
    Components:
        mpz_t finger - The fingerprint itself
        mpz_t r_k    - r^k, where k = ceiling(log_p(finger))
        mpz_t r_mk   - r^-k
*/
typedef struct fingerprint_t {
    mpz_t finger, r_k, r_mk;
} *fingerprint;

/*
    init_fingerprint
    Constructs an empty fingerprint.
    Returns fingerprint:
        finger = 0
        r^k = r^mk = 1
*/
fingerprint init_fingerprint() {
    fingerprint finger = malloc(sizeof(struct fingerprint_t));
    mpz_init(finger->finger);
    mpz_init_set_ui(finger->r_k, 1);
    mpz_init_set_ui(finger->r_mk, 1);
    return finger;
}

/*
    set_fingerprint
    Sets a fingerprint to a given string.
    Parameters:
        fingerprinter printer - The printer to use
        char          *T      - The text string
        unsigned      int l   - The length of the string
        fingerprint   print   - The fingerprint to change
    Returns void:
        Parameter print modified by reference to new fingerprint.
*/
void set_fingerprint(fingerprinter printer, char *T, unsigned int l, fingerprint print) {
    mpz_set_ui(print->r_k, 1);
    int i;

    mpz_set_ui(print->finger, T[0]);
    mpz_mod(print->finger, print->finger, printer->p);

    for (i = 1; i < l; i++) {
        mpz_mul(print->r_k, print->r_k, printer->r);
        mpz_addmul_ui(print->finger, print->r_k, T[i]);
        mpz_mod(print->finger, print->finger, printer->p);
    }
    mpz_mul(print->r_k, print->r_k, printer->r);
    mpz_mod(print->r_k, print->r_k, printer->p);

    mpz_invert(print->r_mk, print->r_k, printer->p);
}

/*
    fingerprint_assign
    Copies a value between fingerprints.
    Parameters:
        fingerprint from - The fingerprint to copy from
        fingerprint to   - The fingerprint to copy to
    Returns void:
        Parameter to modified by reference to copied fingerprint.
*/
void fingerprint_assign(fingerprint from, fingerprint to) {
    mpz_set(to->finger, from->finger);
    mpz_set(to->r_k, from->r_k);
    mpz_set(to->r_mk, from->r_mk);
}

/*
    fingerprint_suffix
    Removes the prefix from a fingerprint.
    Parameters:
        fingerprinter printer - The printer to use
        fingerprint uv        - The total fingerprint
        fingerprint u         - The fingerprint prefix
        fingerprint v         - The fingerprint suffix
    Returns void:
        Parameter v modified by reference to suffix.
*/
void fingerprint_suffix(fingerprinter printer, fingerprint uv, fingerprint u, fingerprint v) {
    mpz_mul(v->r_k, uv->r_k, u->r_mk);
    mpz_mod(v->r_k, v->r_k, printer->p);
    mpz_invert(v->r_mk, v->r_k, printer->p);

    mpz_sub(v->finger, uv->finger, u->finger);
    if (mpz_cmp_si(v->finger, 0) < 0) mpz_add(v->finger, v->finger, printer->p);
    mpz_mul(v->finger, v->finger, u->r_mk);
    mpz_mod(v->finger, v->finger, printer->p);
}

/*
    fingerprint_prefix
    Removes the suffix from a fingerprint.
    Parameters:
        fingerprinter printer - The printer to use
        fingerprint uv        - The total fingerprint
        fingerprint v         - The fingerprint suffix
        fingerprint u         - The fingerprint prefix
    Returns void:
        Parameter u modified by reference to prefix.
*/
void fingerprint_prefix(fingerprinter printer, fingerprint uv, fingerprint v, fingerprint u) {
    mpz_mul(u->r_k, uv->r_k, v->r_mk);
    mpz_mod(u->r_k, u->r_k, printer->p);
    mpz_invert(u->r_mk, u->r_k, printer->p);

    mpz_mul(u->finger, v->finger, u->r_k);
    mpz_mod(u->finger, u->finger, printer->p);
    mpz_sub(u->finger, uv->finger, u->finger);
    if (mpz_cmp_si(u->finger, 0) < 0) mpz_add(u->finger, u->finger, printer->p);
}

/*
    fingerprint_concat
    Concatenates two fingerprints together.
    Parameters:
        fingerprinter printer - The printer to use
        fingerprint u         - The fingerprint prefix
        fingerprint v         - The fingerprint suffix
        fingerprint uv        - The total fingerprint
    Returns void:
        Parameter uv modified by reference to concatenation.
*/
void fingerprint_concat(fingerprinter printer, fingerprint u, fingerprint v, fingerprint uv) {
    mpz_mul(uv->r_k, u->r_k, v->r_k);
    mpz_mod(uv->r_k, uv->r_k, printer->p);
    mpz_invert(uv->r_mk, uv->r_k, printer->p);

    mpz_mul(uv->finger, v->finger, u->r_k);
    mpz_mod(uv->finger, uv->finger, printer->p);
    mpz_add(uv->finger, u->finger, uv->finger);
    if (compare(uv->finger, printer->p) > 0) mpz_sub(uv->finger, uv->finger, printer->p);
}

/*
    fingerprint_equals
    Checks if two fingerprints are equal.
    Parameters:
        fingerprint T_f - The first fingerprint
        fingerprint P_f - The second fingerprint
    Returns int:
        1 if T_f = P_f
        0 otherwise
*/
int fingerprint_equals(fingerprint T_f, fingerprint P_f) {
    return mpz_equals(T_f->finger, P_f->finger);
}

/*
    fingerprint_cmp
    Compares two fingerprints
    Parameters:
        fingerprint T_f - The first fingerprint
        fingerprint P_f - The second fingerprint
    Returns int:
        Positive if T_f > P_f
        0 if T_f = P_f
        Negative if T_f < P_f
*/
int fingerprint_cmp(fingerprint T_f, fingerprint P_f) {
    return mpz_cmp(T_f->finger, P_f->finger);
}

/*
    fingerprint_free
    Frees a fingerprint from memory.
    Parameters:
        fingerprint finger - The fingerprint to free
*/
void fingerprint_free(fingerprint finger) {
    mpz_clear(finger->finger);
    mpz_clear(finger->r_k);
    mpz_clear(finger->r_mk);
    free(finger);
}

#endif
