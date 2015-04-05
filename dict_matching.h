#ifndef __DICT_MATCHING__
#define __DICT_MATCHING__

#include <stdio.h>
#include <stdlib.h>
#include "karp_rabin.h"
#include "hash_lookup.h"
#include "first_lookup.h"
#include "short_dict_matching.h"
#include "periodic_dict_matching.h"
#include "rbtree.c"

typedef struct {
    int row_size;
    int num_patterns;
    int num_complete;
    int cur_prefix;
    int *period;
    int *count;
    int *first_location;
    int *last_location;
    int *end_pattern;
    int *progression_location;
    int *prefix_length;
    hash_lookup lookup;
    hash_lookup *suffixes;
    fingerprint *first_print;
    fingerprint *last_print;
    fingerprint *period_f;
    fingerprint *complete_f;
    rbtree next_progression;
} pattern_row;

int compare_int(void *leftp, void *rightp) {
    int left = (int)leftp;
    int right = (int)rightp;
    if (left < right) return -1;
    else if (left > right) return 1;
    else return 0;
}

int shift_row(fingerprinter printer, pattern_row *row, fingerprint finger, int j, fingerprint tmp) {
    int i = (int)rbtree_lookup(row->next_progression, (void*)j - row->row_size, (void*)-1, compare_int);
    if (i != -1) {
        fingerprint_assign(row->first_print[i], finger);
        if (!row->end_pattern[i]) {
            if (row->count[i] > 1) {
                fingerprint_concat(printer, row->first_print[i], row->period_f[i], tmp);
                fingerprint_assign(tmp, row->first_print[i]);
                row->first_location[i] += row->period[i];
                rbtree_insert(row->next_progression, (void*)j - row->row_size + row->period[i], (void*)i, compare_int);
            }
            rbtree_delete(row->next_progression, (void*)j - row->row_size, compare_int);
            row->count[i]--;
        }
    }
    if (row->num_complete) {
        int del = (int)rbtree_lookup(row->next_progression, (void*)j - row->row_size - row->num_complete, (void*)-1, compare_int);
        if (del != -1) {
            if (row->count[del] > 1) {
                fingerprint_concat(printer, row->first_print[del], row->period_f[del], tmp);
                fingerprint_assign(tmp, row->first_print[del]);
                row->first_location[del] += row->period[del];
                rbtree_insert(row->next_progression, (void*)j - row->row_size - row->num_complete + row->period[del], (void*)del, compare_int);
            }
            rbtree_delete(row->next_progression, (void*)j - row->row_size - row->num_complete, compare_int);
            row->count[del]--;
        }
    }
    return (i != -1) ? 1 : 0;
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
                fprintf(stderr, "Warning: Non-Periodic occurance at %d. Occurance ignored.\n", location);
            }
        }
    } else {
        fingerprint_assign(finger, row->first_print[i]);
        row->first_location[i] = location;
        row->count[i] = 1;
        rbtree_insert(row->next_progression, (void*)location, (void*)i, compare_int);
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
    int num_patterns;
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
            while (k > -1 && P[i][k + 1] != P[i][j]) k = failure[k];
            if (P[i][k + 1] == P[i][j]) k++;
            failure[j] = k;
        }
        period[i] = m[i] - failure[m[i] - 1] - 1;
    }
    free(failure);
}

dict_matcher dict_matching_build(char **P, int *m, int num_patterns, int n, int alpha) {
    dict_matcher matcher = malloc(sizeof(struct dict_matcher_t));
    matcher->printer = fingerprinter_build(n, alpha);
    matcher->T_j = init_fingerprint();
    matcher->tmp = init_fingerprint();
    matcher->T_f = init_fingerprint();
    matcher->num_patterns = num_patterns;

    int *periods = malloc(sizeof(int) * num_patterns);
    get_periods(P, m, num_patterns, periods);

    int m_max = 0;
    int i;
    for (i = 0; i < num_patterns; i++) {
        if (periods[i] > num_patterns) {
            if (m[i] > m_max) {
                m_max = m[i] - num_patterns;
            }
        }
    }
    matcher->num_rows = 0;
    if (m_max > num_patterns) {
        while ((1 << matcher->num_rows) <= m_max) {
            matcher->num_rows++;
        }
        matcher->rows = malloc(sizeof(pattern_row) * matcher->num_rows);

        int j, k, l, lookup_size, old_lookup_size;
        char *first_letters = malloc(sizeof(char) * num_patterns);
        lookup_size = 0;
        for (j = 0; j < num_patterns; j++) {
            if (periods[j] > num_patterns) {
                first_letters[lookup_size] = P[j][0];
                for (k = 0; k < lookup_size; k++) {
                    if (first_letters[k] == first_letters[lookup_size]) break;
                }
                if (k == lookup_size) lookup_size++;
            }
        }
        matcher->first_round = firstlookup_build(first_letters, lookup_size);
        free(first_letters);
        fingerprint *patterns = malloc(sizeof(fingerprint) * num_patterns);
        matcher->T_prev = init_fingerprint();
        matcher->current = init_fingerprint();
        for (i = 0; i < num_patterns; i++) {
            patterns[i] = init_fingerprint();
        }
        old_lookup_size = lookup_size;
        lookup_size = 0;
        int *end_pattern_tmp = malloc(sizeof(int) * num_patterns);
        int *end_pattern = malloc(sizeof(int) * old_lookup_size);
        for (i = 0; i < old_lookup_size; i++) {
            end_pattern[i] = 0;
        }
        int *pattern_hash_id = malloc(sizeof(int) * num_patterns);
        int *progression_location_tmp = malloc(sizeof(int) * num_patterns);
        int *progression_location = NULL;
        int *period_location = NULL;
        int *prefix_length_tmp = malloc(sizeof(int) * num_patterns);
        int *prefix_length = NULL;
        fingerprint *complete_patterns = malloc(sizeof(fingerprint) * num_patterns);
        fingerprint **suffixes = malloc(sizeof(fingerprint*) * num_patterns);
        int *num_suffixes = malloc(sizeof(int) * num_patterns);
        for (i = 0; i < num_patterns; i++) {
            complete_patterns[i] = init_fingerprint();
            suffixes[i] = malloc(sizeof(fingerprint) * num_patterns);
            for (j = 0; j < num_patterns; j++) {
                suffixes[i][j] = init_fingerprint();
            }
        }
        int complete_size = 0;
        i = 0;
        while ((1 << (i + 1)) <= m_max) {
            matcher->rows[i].num_complete = complete_size;
            if (complete_size) {
                matcher->rows[i].complete_f = malloc(sizeof(fingerprint) * complete_size);
                matcher->rows[i].suffixes = malloc(sizeof(hash_lookup) * complete_size);
                for (j = 0; j < complete_size; j++) {
                    matcher->rows[i].complete_f[j] = init_fingerprint();
                    fingerprint_assign(complete_patterns[j], matcher->rows[i].complete_f[j]);
                    matcher->rows[i].suffixes[j] = hashlookup_build(suffixes[j], NULL, num_suffixes[j], matcher->printer);
                }
                matcher->rows[i].prefix_length = prefix_length;
                matcher->rows[i].progression_location = progression_location;
            }
            matcher->rows[i].row_size = 1 << i;
            matcher->rows[i].end_pattern = end_pattern;
            lookup_size = 0;
            complete_size = 0;
            for (j = 0; j < num_patterns; j++) {
                if ((periods[j] > num_patterns) && (m[j] - num_patterns >= matcher->rows[i].row_size << 1)) {
                    set_fingerprint(matcher->printer, P[j], matcher->rows[i].row_size << 1, patterns[lookup_size]);
                    if (m[j] - num_patterns < (matcher->rows[i].row_size << 2)) {
                        end_pattern_tmp[lookup_size] = 1;
                        set_fingerprint(matcher->printer, P[j], m[j] - num_patterns, complete_patterns[complete_size]);
                        progression_location_tmp[complete_size] = lookup_size;
                        prefix_length_tmp[complete_size] = m[j] - num_patterns;
                        num_suffixes[complete_size] = 1;
                        set_fingerprint(matcher->printer, &P[j][m[j] - num_patterns], num_patterns, suffixes[complete_size][0]);
                        for (k = 0; k < complete_size; k++) {
                            if (fingerprint_equals(complete_patterns[k], complete_patterns[complete_size])) {
                                end_pattern_tmp[lookup_size] = 0;
                                for (l = 0; l < num_suffixes[k]; l++) {
                                    if (fingerprint_equals(suffixes[k][l], suffixes[complete_size][0])) {
                                        break;
                                    }
                                    if (l == num_suffixes[k]) {
                                        fingerprint_assign(suffixes[complete_size][0], suffixes[k][l]);
                                        num_suffixes[k]++;
                                    }
                                }
                                break;
                            }
                        }
                        if (k == complete_size) complete_size++;
                    } else {
                        end_pattern_tmp[lookup_size] = 0;
                    }
                    for (k = 0; k < lookup_size; k++) {
                        if (fingerprint_equals(patterns[k], patterns[lookup_size])) {
                            end_pattern_tmp[k] |= end_pattern_tmp[lookup_size];
                            if (end_pattern_tmp[lookup_size]) progression_location_tmp[complete_size - 1] = k;
                            break;
                        }
                    }
                    if (k == lookup_size) lookup_size++;
                }
            }
            matcher->rows[i].lookup = hashlookup_build(patterns, NULL, lookup_size, matcher->printer);
            matcher->rows[i].num_patterns = old_lookup_size;
            matcher->rows[i].first_print = malloc(sizeof(fingerprint) * old_lookup_size);
            matcher->rows[i].first_location = malloc(sizeof(int) * old_lookup_size);
            matcher->rows[i].last_print = malloc(sizeof(fingerprint) * old_lookup_size);
            matcher->rows[i].last_location = malloc(sizeof(int) * old_lookup_size);
            matcher->rows[i].period = malloc(sizeof(int) * old_lookup_size);
            matcher->rows[i].count = malloc(sizeof(int) * old_lookup_size);
            matcher->rows[i].period_f = malloc(sizeof(fingerprint) * old_lookup_size);
            matcher->rows[i].next_progression = rbtree_create();
            matcher->rows[i].cur_prefix = 0;
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
            end_pattern = malloc(sizeof(int) * lookup_size);
            for (j = 0; j < lookup_size; j++) {
                int location = hashlookup_search(matcher->rows[i].lookup, patterns[j], NULL);
                pattern_hash_id[j] = location;
                end_pattern[location] = end_pattern_tmp[j];
            }
            if (complete_size) {
                period_location = malloc(sizeof(int) * complete_size);
                prefix_length = malloc(sizeof(int) * complete_size);
                progression_location = malloc(sizeof(int) * complete_size);
                for (j = 0; j < complete_size; j++) {
                    period_location[j] = pattern_hash_id[progression_location_tmp[j]];
                    progression_location[j] = progression_location_tmp[j];
                    prefix_length[j] = prefix_length_tmp[j];
                }
            }
            i++;
        }
        matcher->rows[i].num_complete = complete_size;
        if (complete_size) {
            matcher->rows[i].complete_f = malloc(sizeof(fingerprint) * complete_size);
            matcher->rows[i].suffixes = malloc(sizeof(hash_lookup) * complete_size);
            for (j = 0; j < complete_size; j++) {
                matcher->rows[i].complete_f[j] = init_fingerprint();
                fingerprint_assign(complete_patterns[j], matcher->rows[i].complete_f[j]);
                matcher->rows[i].suffixes[j] = hashlookup_build(suffixes[j], NULL, num_suffixes[j], matcher->printer);
            }
            matcher->rows[i].prefix_length = prefix_length;
            matcher->rows[i].progression_location = progression_location;
        }
        matcher->rows[i].row_size = 1 << i;
        matcher->rows[i].end_pattern = end_pattern;
        matcher->rows[i].num_patterns = old_lookup_size;
        matcher->rows[i].first_print = malloc(sizeof(fingerprint) * old_lookup_size);
        matcher->rows[i].first_location = malloc(sizeof(int) * old_lookup_size);
        matcher->rows[i].last_print = malloc(sizeof(fingerprint) * old_lookup_size);
        matcher->rows[i].last_location = malloc(sizeof(int) * old_lookup_size);
        matcher->rows[i].period = malloc(sizeof(int) * old_lookup_size);
        matcher->rows[i].count = malloc(sizeof(int) * old_lookup_size);
        matcher->rows[i].period_f = malloc(sizeof(fingerprint) * old_lookup_size);
        matcher->rows[i].next_progression = rbtree_create();
        matcher->rows[i].cur_prefix = 0;
        for (j = 0; j < old_lookup_size; j++) {
            matcher->rows[i].first_print[j] = init_fingerprint();
            matcher->rows[i].first_location[j] = 0;
            matcher->rows[i].last_print[j] = init_fingerprint();
            matcher->rows[i].last_location[j] = 0;
            matcher->rows[i].period[j] = 0;
            matcher->rows[i].count[j] = 0;
            matcher->rows[i].period_f[j] = init_fingerprint();
        }

        for (i = 0; i < num_patterns; i++) {
            fingerprint_free(patterns[i]);
        }
        free(patterns);
        free(end_pattern_tmp);
        free(progression_location_tmp);
        free(prefix_length_tmp);
        for (i = 0; i < num_patterns; i++) {
            fingerprint_free(complete_patterns[i]);
            for (j = 0; j < num_patterns; j++) {
                fingerprint_free(suffixes[i][j]);
            }
            free(suffixes[i]);
        }
        free(suffixes);
        free(num_suffixes);
        free(complete_patterns);
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

    int short_result, periodic_result;

    if (matcher->num_rows) {
        int i, occurance, match;
        for (i = matcher->num_rows - 1; i >= 0; i--) {
            if (matcher->rows[i].num_complete) {
                int cur_prefix = matcher->rows[i].cur_prefix;
                int cur_progression = matcher->rows[i].progression_location[cur_prefix];
                if ((matcher->rows[i].count[cur_progression]) && (matcher->rows[i].first_location[cur_progression] <= (j - matcher->rows[i].prefix_length[cur_prefix] + (matcher->rows[i].row_size << 1)))) {
                    printf("%d!\n", j);
                }

                if (matcher->rows[i].num_complete > 1) {
                    if (++cur_prefix == matcher->rows[i].num_complete) cur_prefix = 0;
                    cur_progression = matcher->rows[i].progression_location[cur_prefix];
                    if ((matcher->rows[i].count[cur_progression]) && (matcher->rows[i].first_location[cur_progression] <= (j - matcher->rows[i].prefix_length[cur_prefix] + (matcher->rows[i].row_size << 1)))) {
                        printf("%d!\n", j);
                    }
                    matcher->rows[i].cur_prefix = (++cur_prefix == matcher->rows[i].num_complete) ? 0 : cur_prefix;
                }
            }

            occurance = shift_row(matcher->printer, &matcher->rows[i], matcher->current, j, matcher->tmp);
            if ((i < matcher->num_rows - 1) && (occurance)) {
                fingerprint_suffix(matcher->printer, matcher->T_f, matcher->current, matcher->tmp);
                match = hashlookup_search(matcher->rows[i].lookup, matcher->tmp, NULL);
                if (match != -1) {
                    add_occurance(matcher->printer, &matcher->rows[i + 1], matcher->current, j, match, matcher->tmp);
                }
            }
        }

        int first_round = firstlookup_search(matcher->first_round, T_j);
        if (first_round != -1) {
            add_occurance(matcher->printer, &matcher->rows[0], matcher->T_prev, j, first_round, matcher->tmp);
        }
    }

    short_result = short_dict_matching_stream(matcher->short_matcher, matcher->printer, matcher->T_f, matcher->tmp, j);
    periodic_result = periodic_dict_matching_stream(matcher->periodic_matcher, matcher->printer, matcher->T_f, matcher->tmp, j);

    return ((short_result != -1) || (periodic_result != -1)) ? j : -1;
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
            rbtree_destroy(matcher->rows[i].next_progression);

            if (i != matcher->num_rows - 1) hashlookup_free(&matcher->rows[i].lookup);

            if (matcher->rows[i].num_complete) {
                for (j = 0; j < matcher->rows[i].num_complete; j++) {
                    fingerprint_free(matcher->rows[i].complete_f[j]);
                    hashlookup_free(&matcher->rows[i].suffixes[j]);
                }
                free(matcher->rows[i].suffixes);
                free(matcher->rows[i].complete_f);
                free(matcher->rows[i].prefix_length);
                free(matcher->rows[i].progression_location);
            }
            free(matcher->rows[i].end_pattern);
        }
        free(matcher->rows);
    }

    short_dict_matching_free(matcher->short_matcher);
    periodic_dict_matching_free(matcher->periodic_matcher);

    free(matcher);
}

#endif
