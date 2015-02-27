#include "karp_rabin.h"
#include "hash_lookup.h"

typedef struct short_dict_matcher_t {
    hash_lookup lookup;
    fingerprint *t_prev;
    int k, start;
} *short_dict_matcher;

int suffix(char *t, int n, char **p, int *m, int num_patterns) {
    int i, j, diff;
    for (i = 0; i < num_patterns; i++) {
        diff = n - m[i];
        for (j = 0; j < m[i]; j++) {
            if (t[diff + j] != p[i][j]) {
                break;
            }
            if (j == m[i]) return 1;
        }
    }
    return 0;
}

short_dict_matcher short_dict_matching_build(int k, fingerprinter printer, char **p, int *m) {
    short_dict_matcher state = malloc(sizeof(struct short_dict_matcher_t));
    int i = 1, j, k_p = 0;
    for (i = 1; i != 0; i <<= 1) {
        if (k & i) {
            k_p = i;
        }
    }
    if (k_p ^ k) k_p <<= 1;
    state->k = k_p;
    state->t_prev = malloc(sizeof(fingerprint) * k_p);
    for (i = 0; i < k_p; i++) {
        state->t_prev[i] = init_fingerprint();
    }
    fingerprint *pattern_prints = malloc(sizeof(fingerprint) * k_p * k_p);
    int *suffix_match = malloc(sizeof(int) * k_p * k_p);
    for (i = 0; i < k_p * k_p; i++) {
        pattern_prints[i] = init_fingerprint();
    }

    int num_prints = 0, x, y;
    for (i = 0; i < k; i++) {
        x = k_p;
        y = m[i];
        while (x != 0) {
            if (y >= (x >> 1)) {
                set_fingerprint(printer, &p[i][m[i] - (x >> 1)], (x >> 1), pattern_prints[num_prints]);
                suffix_match[num_prints] = suffix(&p[i][m[i] - (x >> 1)], (x >> 1), p, m, k);
                for (j = 0; j < num_prints; j++) {
                    if (fingerprint_equals(pattern_prints[num_prints], pattern_prints[j])) {
                        break;
                    }
                }
                if (j == num_prints) {
                    num_prints++;
                }
                y -= x >> 1;
            }
            x >>= 1;
        }
    }

    state->lookup = hashlookup_build(pattern_prints, suffix_match, num_prints, printer);

    for (i = 0; i < k_p * k_p; i++) {
        fingerprint_free(pattern_prints[i]);
    }
    free(pattern_prints);
    free(suffix_match);
    return state;
}

void short_dict_matching_free(short_dict_matcher state) {
    int i;
    for (i = 0; i < state->k; i++) {
        fingerprint_free(state->t_prev[i]);
    }
    free(state->t_prev);
    free(state);
}