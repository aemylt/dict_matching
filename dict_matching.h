#ifndef __DICT_MATCHING__
#define __DICT_MATCHING__

#include <stdio.h>
#include <stdlib.h>
#include "karp_rabin.h"
#include "hash_lookup.h"
#include "first_lookup.h"
#include "short_dict_matching.h"
#include "periodic_dict_matching.h"

typedef struct {
    int row_size, num_patterns, *period, *count, *first_location, *last_location;
    hash_lookup lookup;
    fingerprint *first_print, *last_print;
    fingerprint *period_f;
} pattern_row;

int shift_row(fingerprinter printer, pattern_row *row, fingerprint finger, int j, fingerprint tmp) {
    int i;
    for (i = 0; i < row->num_patterns; i++) {
        if ((row->count[i]) && (row->first_location[i] + row->row_size == j)) {
            fingerprint_assign(row->first_print[i], finger);
            if (row->count[i] > 1) {
                fingerprint_concat(printer, row->first_print[i], row->period_f[i], tmp);
                fingerprint_assign(tmp, row->first_print[i]);
                row->first_location[i] += row->period[i];
            }
            row->count[i]--;
            return 1;
        }
    }
    return 0;
}

void add_occurance(fingerprinter printer, pattern_row *row, fingerprint finger, int location, int i, fingerprint tmp) {
    if (row->count[i]) {
        if (row->count[i] == 1) {
            row->period[i] = location - row->first_location[i];
            fingerprint_suffix(printer, finger, row->first_print[i], row->period_f[i]);
            row->last_location[i] = location;
            fingerprint_assign(finger, row->last_print[i]);
            row->count[i] = 2;
        } else {
            fingerprint_suffix(printer, finger, row->last_print[i], tmp);
            int period = location - row->last_location[i];
            if ((period == row->period[i]) && (fingerprint_equals(tmp, row->period_f[i]))) {
                row->last_location[i] = location;
                fingerprint_assign(finger, row->last_print[i]);
                row->count[i]++;
            } else {
                fprintf(stderr, "Warning: Non-Periodic occurance spotted. Occurance ignored.");
            }
        }
    } else {
        fingerprint_assign(finger, row->first_print[i]);
        row->first_location[i] = location;
        row->count[i] = 1;
    }
}

typedef struct dict_matcher_t {
    fingerprinter printer;
    fingerprint T_j;
    fingerprint tmp;

    short_dict_matcher short_matcher;
    periodic_dict_matcher periodic_matcher;

    pattern_row *rows;
    int num_rows, *end_pattern;
    first_lookup first_round;
    fingerprint T_f;
    fingerprint T_prev;
    fingerprint current;
} *dict_matcher;

void get_periods(char **P, int *m, int num_patterns, int *period) {
    int i, j, k;
    int *failure;
    failure = malloc(sizeof(int));
    int failure_size = 1;
    for (i = 0; i < num_patterns; i++) {
        if (m[i] > failure_size) {
            failure = realloc(failure, sizeof(int) * m[i]);
            failure_size = m[i];
        }
        failure[0] = -1;
        k = -1;
        for (j = 1; j < m[i]; j++) {
            while (k > -1 && P[k + 1] != P[j]) k = failure[k];
            if (P[k + 1] == P[j]) k++;
            failure[j] = k;
        }
        period[i] = failure[m[i] - 1];
    }
    free(failure);
}

dict_matcher dict_matching_build(char **P, int *m, int num_patterns, int n, int alpha) {
    dict_matcher matcher = malloc(sizeof(struct dict_matcher_t));
    matcher->printer = fingerprinter_build(n, alpha);
    matcher->T_j = init_fingerprint();
    matcher->tmp = init_fingerprint();
    matcher->T_f = init_fingerprint();

    int *periods = malloc(sizeof(int) * num_patterns);
    get_periods(P, m, num_patterns, periods);

    int m_max = 0;
    int i;
    for (i = 0; i < num_patterns; i++) {
        if (((periods[i] == -1) && (m[i] > num_patterns)) || (periods[i] > num_patterns)) {
            if (m[i] > m_max) {
                m_max = m[i];
            }
        }
    }
    matcher->num_rows = 0;
    if (m_max > num_patterns) {
        while ((1 << matcher->num_rows) < m_max) {
            matcher->num_rows++;
        }
        matcher->rows = malloc(sizeof(pattern_row) * matcher->num_rows);

        int j, k, lookup_size, old_lookup_size;
        char *first_letters = malloc(sizeof(char) * num_patterns);
        int *end_pattern = malloc(sizeof(int) * num_patterns);
        lookup_size = 0;
        for (j = 0; j < num_patterns; j++) {
            if (((periods[j] == -1) && (m[j] > num_patterns)) || (periods[j] > num_patterns)) {
                first_letters[lookup_size] = P[j][0];
                for (k = 0; k < lookup_size; k++) {
                    if (first_letters[k] == first_letters[lookup_size]) break;
                }
                if (k == lookup_size) lookup_size++;
                if (m[j] == 1) end_pattern[k] = 1;
                else end_pattern[k] = 0;
            }
        }
        matcher->first_round = firstlookup_build(first_letters, end_pattern, lookup_size);
        free(first_letters);
        fingerprint *patterns = malloc(sizeof(fingerprint) * num_patterns);
        matcher->T_prev = init_fingerprint();
        matcher->current = init_fingerprint();
        for (i = 0; i < num_patterns; i++) {
            patterns[i] = init_fingerprint();
        }
        i = 0;
        old_lookup_size = lookup_size;
        lookup_size = 0;
        while ((1 << (i + 1)) <= m_max) {
            matcher->rows[i].row_size = 1 << i;
            lookup_size = 0;
            for (j = 0; j < num_patterns; j++) {
                if ((((periods[j] == -1) && (m[j] > num_patterns)) || (periods[j] > num_patterns)) && (m[j] >= matcher->rows[i].row_size << 1)) {
                    set_fingerprint(matcher->printer, P[j], matcher->rows[i].row_size << 1, patterns[lookup_size]);
                    if (m[j] == matcher->rows[i].row_size << 1) end_pattern[lookup_size] = 1;
                    else end_pattern[lookup_size] = 0;
                    for (k = 0; k < lookup_size; k++) {
                        if (fingerprint_equals(patterns[k], patterns[lookup_size])) {
                            end_pattern[k] |= end_pattern[lookup_size];
                            break;
                        }
                    }
                    if (k == lookup_size) lookup_size++;
                }
            }
            matcher->rows[i].lookup = hashlookup_build(patterns, end_pattern, lookup_size, matcher->printer);
            matcher->rows[i].num_patterns = old_lookup_size;
            matcher->rows[i].first_print = malloc(sizeof(fingerprint) * old_lookup_size);
            matcher->rows[i].first_location = malloc(sizeof(int) * old_lookup_size);
            matcher->rows[i].last_print = malloc(sizeof(fingerprint) * old_lookup_size);
            matcher->rows[i].last_location = malloc(sizeof(int) * old_lookup_size);
            matcher->rows[i].period = malloc(sizeof(int) * old_lookup_size);
            matcher->rows[i].count = malloc(sizeof(int) * old_lookup_size);
            matcher->rows[i].period_f = malloc(sizeof(fingerprint) * old_lookup_size);
            for (j = 0; j < old_lookup_size; j++) {
                matcher->rows[i].first_print[j] = init_fingerprint();
                matcher->rows[i].first_location[j] = 0;
                matcher->rows[i].last_print[j] = init_fingerprint();
                matcher->rows[i].last_location[j] = 0;
                matcher->rows[i].period[j] = 0;
                matcher->rows[i].count[j] = 0;
                matcher->rows[i].period_f[j] = init_fingerprint();
            }
            old_lookup_size = lookup_size;
            i++;
        }

        for (i = 0; i < num_patterns; i++) {
            fingerprint_free(patterns[i]);
        }
        free(patterns);
        free(end_pattern);
    }

    matcher->short_matcher = short_dict_matching_build(num_patterns, matcher->printer, P, m);
    matcher->periodic_matcher = periodic_dict_matching_build(P, m, periods, num_patterns, matcher->printer);

    return matcher;
}

int dict_matching_stream(dict_matcher matcher, char T_j, int j) {
    if (matcher->num_rows) fingerprint_assign(matcher->T_f, matcher->T_prev);
    set_fingerprint(matcher->printer, &T_j, 1, matcher->T_j);
    fingerprint_concat(matcher->printer, matcher->T_f, matcher->T_j, matcher->tmp);
    fingerprint_assign(matcher->tmp, matcher->T_f);

    int long_result = 0, short_result, periodic_result;

    if (matcher->num_rows) {
        int i, occurance, match;
        for (i = matcher->num_rows - 1; i >= 0; i--) {
            occurance = shift_row(matcher->printer, &matcher->rows[i], matcher->current, j, matcher->tmp);
            if (occurance) {
                fingerprint_suffix(matcher->printer, matcher->T_f, matcher->current, matcher->tmp);
                match = hashlookup_search(matcher->rows[i].lookup, matcher->tmp, &long_result);
                if (match != -1) {
                    if (i != matcher->num_rows - 1) {
                        add_occurance(matcher->printer, &matcher->rows[i + 1], matcher->current, j, match, matcher->tmp);
                    }
                    else {
                        long_result = 1;
                    }
                }
            }
        }

        int first_round = firstlookup_search(matcher->first_round, T_j, &long_result);
        if (first_round != -1) {
            add_occurance(matcher->printer, &matcher->rows[0], matcher->T_prev, j, first_round, matcher->tmp);
        }
    }

    short_result = short_dict_matching_stream(matcher->short_matcher, matcher->printer, matcher->T_f, matcher->tmp, j);
    periodic_result = periodic_dict_matching_stream(matcher->periodic_matcher, matcher->printer, matcher->T_f, matcher->tmp, j);

    return (long_result || (short_result != -1) || (periodic_result != -1)) ? j : -1;
}

void dict_matching_free(dict_matcher matcher) {
    fingerprint_free(matcher->T_j);
    fingerprint_free(matcher->tmp);
    fingerprinter_free(matcher->printer);
    fingerprint_free(matcher->T_f);

    if (matcher->num_rows){
        fingerprint_free(matcher->T_prev);
        fingerprint_free(matcher->current);
        firstlookup_free(&matcher->first_round);
        int i, j;
        for (i = 0; i < matcher->num_rows; i++) {
            for (j = 0; j < matcher->rows[i].num_patterns; j++) {
                fingerprint_free(matcher->rows[i].first_print[j]);
                fingerprint_free(matcher->rows[i].last_print[j]);
                fingerprint_free(matcher->rows[i].period_f[j]);
            }
            free(matcher->rows[i].first_print);
            free(matcher->rows[i].first_location);
            free(matcher->rows[i].last_print);
            free(matcher->rows[i].last_location);
            free(matcher->rows[i].period_f);
            free(matcher->rows[i].period);
            free(matcher->rows[i].count);

            hashlookup_free(&matcher->rows[i].lookup);
        }
        free(matcher->rows);
    }

    short_dict_matching_free(matcher->short_matcher);
    periodic_dict_matching_free(matcher->periodic_matcher);

    free(matcher);
}

#endif