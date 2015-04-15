#ifndef __PERIODIC_DICT_MATCHING__
#define __PERIODIC_DICT_MATCHING__

#include "karp_rabin.h"
#include "hash_lookup.h"

typedef struct periodic_dict_matcher_t {
    int num_heads, num_tails, k, *period, *count, *first_location, *last_location, *m;
    hash_lookup head, tail;
    fingerprint *first_print;
    fingerprint *last_print;
    fingerprint *period_f;
} *periodic_dict_matcher;

periodic_dict_matcher periodic_dict_matching_build(char **P, int *m, int *period, int k, fingerprinter printer) {
    periodic_dict_matcher state = malloc(sizeof(struct periodic_dict_matcher_t));
    state->k = k;
    int i, j, num_heads = 0, num_tails = 0, location;
    fingerprint *prints = malloc(sizeof(fingerprint) * k);
    fingerprint *period_prints = malloc(sizeof(fingerprint) * k);
    int *periods = malloc(sizeof(int) * k);
    for (i = 0; i < k; i++) {
        prints[i] = init_fingerprint();
        period_prints[i] = init_fingerprint();
        if ((period[i] <= k) && (m[i] > k)) {
            set_fingerprint(printer, P[i], k, prints[num_heads]);
            set_fingerprint(printer, &P[i][k - period[i]], period[i], period_prints[num_heads]);
            periods[num_heads] = period[i];
            for (j = 0; j < num_heads; j++) {
                if (fingerprint_equals(prints[j], prints[num_heads])) break;
            }
            if (j == num_heads) num_heads++;
        }
    }
    state->num_heads = num_heads;
    state->head = hashlookup_build(prints, NULL, num_heads, printer);
    if (num_heads) {
        state->period = malloc(sizeof(int) * num_heads);
        state->count = malloc(sizeof(int) * num_heads);
        state->first_location = malloc(sizeof(int) * num_heads);
        state->last_location = malloc(sizeof(int) * num_heads);
        state->first_print = malloc(sizeof(fingerprint) * num_heads);
        state->last_print = malloc(sizeof(fingerprint) * num_heads);
        state->period_f = malloc(sizeof(fingerprint) * num_heads);
        for (i = 0; i < num_heads; i++) {
            location = hashlookup_search(state->head, prints[i], NULL);
            state->period[location] = periods[i];
            state->count[location] = 0;
            state->first_location[location] = 0;
            state->last_location[location] = 0;
            state->first_print[location] = init_fingerprint();
            state->last_print[location] = init_fingerprint();
            state->period_f[location] = init_fingerprint();
            fingerprint_assign(period_prints[i], state->period_f[location]);
        }
    }
    int *locations = malloc(sizeof(int) * k);
    int *lengths = malloc(sizeof(int) * k);
    for (i = 0; i < k; i++) {
        if ((period[i] <= k) && (m[i] > (k << 1))) {
            set_fingerprint(printer, &P[i][m[i] - k], k, prints[num_tails]);
            set_fingerprint(printer, P[i], k, period_prints[0]);
            locations[num_tails] = hashlookup_search(state->head, period_prints[0], NULL);
            lengths[num_tails] = m[i];
            for (j = 0; j < num_tails; j++) {
                if (fingerprint_equals(prints[j], prints[num_tails])) {
                    if (lengths[num_tails] < lengths[j]) {
                        locations[j] = locations[num_tails];
                        lengths[j] = lengths[num_tails];
                    }
                    break;
                }
            }
            if (j == num_tails) num_tails++;
        }
    }
    state->num_tails = num_tails;
    state->tail = hashlookup_build(prints, locations, num_tails, printer);
    if (num_tails) {
        state->m = malloc(sizeof(int) * num_tails);
        for (i = 0; i < num_tails; i++) {
            location = hashlookup_search(state->tail, prints[i], NULL);
            state->m[location] = lengths[i];
        }
    }
    free(locations);
    free(lengths);
    for (i = 0; i < k; i++) {
        fingerprint_free(prints[i]);
        fingerprint_free(period_prints[i]);
    }
    free(prints);
    free(period_prints);
    free(periods);
    return state;
}

int periodic_dict_matching_stream(periodic_dict_matcher state, fingerprinter printer, fingerprint t_f, fingerprint *t_prev, fingerprint tmp, int j) {
    fingerprint_suffix(printer, t_f, t_prev[(j + state->k) % (state->k << 1)], tmp);
    int head_location = hashlookup_search(state->head, tmp, NULL);
    int head_pointer = 0;
    int tail_location = hashlookup_search(state->tail, tmp, &head_pointer);
    if (head_location != -1) {
        int period = j - state->last_location[head_location];
        fingerprint_suffix(printer, t_f, state->last_print[head_location], tmp);
        if ((period == state->period[head_location]) && (fingerprint_equals(tmp, state->period_f[head_location]))) {
            state->count[head_location]++;
            state->last_location[head_location] = j;
            fingerprint_assign(t_f, state->last_print[head_location]);
        } else {
            state->count[head_location] = 0;
        }
        if (!state->count[head_location]) {
            state->first_location[head_location] = j;
            fingerprint_assign(t_f, state->first_print[head_location]);
            state->last_location[head_location] = j;
            fingerprint_assign(t_f, state->last_print[head_location]);
            state->count[head_location] = state->k / state->period[head_pointer];
        }
    }
    int result = -1;
    if (tail_location != -1) {
        int num_occurances = (state->m[tail_location] - (state->k % state->period[head_pointer])) / state->period[head_pointer];
        int last_occurance = j - (state->m[tail_location] - state->k) % state->period[head_pointer];
        if ((num_occurances <= state->count[head_pointer]) && (last_occurance == state->last_location[head_pointer])) {
            result = j;
        }
    }
    return result;
}

void periodic_dict_matching_free(periodic_dict_matcher state) {
    int i;
    for (i = 0; i < state->num_heads; i++) {
        fingerprint_free(state->first_print[i]);
        fingerprint_free(state->last_print[i]);
        fingerprint_free(state->period_f[i]);
    }
    if (state->num_heads) {
        free(state->first_print);
        free(state->last_print);
        free(state->period_f);
        free(state->period);
        free(state->count);
        free(state->first_location);
        free(state->last_location);
        hashlookup_free(&state->head);
    }
    if (state->num_tails) {
        free(state->m);
        hashlookup_free(&state->tail);
    }
    free(state);
}

int periodic_dict_matching_size(periodic_dict_matcher state) {
    int size = sizeof(struct periodic_dict_matcher_t) + sizeof(int) * (4 * state->num_heads + state->num_tails) + sizeof(fingerprint) * (3 * state->num_heads) + hashlookup_size(state->head) + hashlookup_size(state->tail);
    int i;
    for (i = 0; i < state->num_heads; i++) {
        size += fingerprint_size(state->first_print[i]) + fingerprint_size(state->last_print[i]) + fingerprint_size(state->period_f[i]);
    }
    return size;
}

#endif