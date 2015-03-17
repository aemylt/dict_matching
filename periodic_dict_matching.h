#ifndef __PERIODIC_DICT_MATCHING__
#define __PERIODIC_DICT_MATCHING__

#include "karp_rabin.h"
#include "hash_lookup.h"

typedef struct periodic_dict_matching_state_t {
    int num_heads, num_tails, k, *period, *count, *first_location, *last_location;
    hash_lookup head, tail;
    fingerprint *first_print;
    fingerprint *last_print;
    fingerprint *period_f;
    fingerprint *t_prev;
} *periodic_dict_matching_state;

periodic_dict_matching_state periodic_dict_matching_build(char **P, int *m, int *period, int k, fingerprinter printer) {
    periodic_dict_matching_state state = malloc(sizeof(struct periodic_dict_matching_state_t));
    state->k = k;
    state->t_prev = malloc(sizeof(fingerprint) * k);
    int i, j, num_heads = 0, num_tails = 0;
    fingerprint *prints = malloc(sizeof(fingerprint) * k);
    for (i = 0; i < k; i++) {
        state->t_prev[i] = init_fingerprint();
    }
    for (i = 0; i < k; i++) {
        prints[i] = init_fingerprint();
        set_fingerprint(printer, P[i], k, prints[num_heads]);
        for (j = 0; j < num_heads; j++) {
            if (fingerprint_equals(prints[j], prints[num_heads])) break;
        }
        if (j == num_heads) num_heads++;
    }
    state->head = hashlookup_build(prints, NULL, num_heads, printer);
    state->num_heads = num_heads;
    state->period = malloc(sizeof(int) * num_heads);
    state->count = malloc(sizeof(int) * num_heads);
    state->first_location = malloc(sizeof(int) * num_heads);
    state->last_location = malloc(sizeof(int) * num_heads);
    state->first_print = malloc(sizeof(fingerprint) * num_heads);
    state->last_print = malloc(sizeof(fingerprint) * num_heads);
    state->period_f = malloc(sizeof(fingerprint) * num_heads);
    for (i = 0; i < num_heads; i++) {
        state->period[i] = 0;
        state->count[i] = 0;
        state->first_location[i] = 0;
        state->last_location[i] = 0;
        state->first_print[i] = init_fingerprint();
        state->last_print[i] = init_fingerprint();
        state->period_f[i] = init_fingerprint();
    }
    for (i = 0; i < k; i++) {
        set_fingerprint(printer, P[i], k, prints[num_tails]);
        for (j = 0; j < num_tails; j++) {
            if (fingerprint_equals(prints[j], prints[num_tails])) break;
        }
        if (j == num_tails) num_tails++;
    }
    state->tail = hashlookup_build(prints, NULL, num_tails, printer);
    state->num_tails = num_tails;
    for (i = 0; i < k; i++) {
        fingerprint_free(prints[i]);
    }
    free(prints);
    return state;
}

periodic_dict_matching_free(periodic_dict_matching_state state) {
    int i;
    for (i = 0; i < state->k; i++) {
        fingerprint_free(state->t_prev[i]);
    }
    free(state->t_prev);
    for (i = 0; i < state->num_heads; i++) {
        fingerprint_free(state->first_print[i]);
        fingerprint_free(state->last_print[i]);
        fingerprint_free(state->period_f[i]);
    }
    free(state->first_print);
    free(state->last_print);
    free(state->period_f);
    free(state->period);
    free(state->count);
    free(state->first_location);
    free(state->last_location);
    hashlookup_free(state->head);
    hashlookup_free(state->tail);
    free(state);
}

#endif