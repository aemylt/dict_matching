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
    int *period;
    int *count;
    int *first_location;
    int *last_location;
    hash_lookup lookup;
    fingerprint *first_print;
    fingerprint *last_print;
    fingerprint *period_f;
    rbtree next_progression;
} pattern_row;

typedef struct {
    int num_prefixes;
    int num_progressions;
    int cur_prefix;
    int *first_location;
    int *last_location;
    int *period;
    int *count;
    int *progression_index;
    int *prefix_length;
    int *num_suffixes;
    fingerprint **suffixes;
    fingerprint *first_print;
    fingerprint *last_print;
    fingerprint *period_f;
    fingerprint *prefix;
    rbtree *suffix_tree;
} final_row;

int compare_int(void *leftp, void *rightp) {
    int left = (int)leftp;
    int right = (int)rightp;
    if (left < right) return -1;
    else if (left > right) return 1;
    else return 0;
}

int compare_finger(void *leftp, void *rightp) {
    fingerprint left = (fingerprint)leftp;
    fingerprint right = (fingerprint)rightp;
    return fingerprint_cmp(left, right);
}

int shift_row(fingerprinter printer, pattern_row *row, fingerprint finger, int j, fingerprint tmp) {
    int i = (int)rbtree_lookup(row->next_progression, (void*)j - row->row_size, (void*)-1, compare_int);
    if (i != -1) {
        fingerprint_assign(row->first_print[i], finger);
        if (row->count[i] > 1) {
            fingerprint_concat(printer, row->first_print[i], row->period_f[i], tmp);
            fingerprint_assign(tmp, row->first_print[i]);
            row->first_location[i] += row->period[i];
            rbtree_insert(row->next_progression, (void*)j - row->row_size + row->period[i], (void*)i, compare_int);
        }
        rbtree_delete(row->next_progression, (void*)j - row->row_size, compare_int);
        row->count[i]--;
        return 1;
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
                fprintf(stderr, "Warning: Non-Periodic occurance at %d in index %d. Occurance ignored.\n", location, i);
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
    final_row final;
    int num_rows;
    first_lookup first_round;
    fingerprint T_f;
    fingerprint *T_prev;
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
        if ((periods[i] > num_patterns) && (m[i] > num_patterns << 1)) {
            if (m[i] > m_max) {
                m_max = m[i] - num_patterns;
            }
        }
    }
    matcher->num_rows = 0;
    if (m_max > num_patterns) {
        while ((1 << (matcher->num_rows + 1)) < m_max) {
            matcher->num_rows++;
        }
        matcher->rows = malloc(sizeof(pattern_row) * matcher->num_rows);

        int j, k, l, lookup_size, old_lookup_size;
        char *first_letters = malloc(sizeof(char) * num_patterns);
        fingerprint *patterns = malloc(sizeof(fingerprint) * num_patterns);
        fingerprint *old_patterns = malloc(sizeof(fingerprint) * num_patterns);
        int *prev_row = malloc(sizeof(int) * num_patterns);
        lookup_size = 0;
        for (j = 0; j < num_patterns; j++) {
            patterns[j] = init_fingerprint();
            old_patterns[j] = init_fingerprint();
            if ((periods[j] > num_patterns) && (m[j] > num_patterns << 1)) {
                first_letters[lookup_size] = P[j][0];
                set_fingerprint(matcher->printer, P[j], 1, old_patterns[lookup_size]);
                prev_row[j] = lookup_size;
                for (k = 0; k < lookup_size; k++) {
                    if (first_letters[k] == first_letters[lookup_size]) {
                        prev_row[j] = k;
                        break;
                    }
                }
                if (k == lookup_size) lookup_size++;
            }
        }
        matcher->first_round = firstlookup_build(first_letters, lookup_size);
        free(first_letters);
        matcher->current = init_fingerprint();
        old_lookup_size = lookup_size;
        lookup_size = 0;
        int *end_pattern = malloc(sizeof(int) * num_patterns);
        int *progression_index = malloc(sizeof(int) * num_patterns);
        int *prefix_length = malloc(sizeof(int) * num_patterns);
        fingerprint *prefix = malloc(sizeof(fingerprint) * num_patterns);
        fingerprint **suffixes = malloc(sizeof(fingerprint*) * num_patterns);
        int *num_suffixes = malloc(sizeof(int) * num_patterns);
        for (i = 0; i < num_patterns; i++) {
            prefix[i] = init_fingerprint();
            suffixes[i] = malloc(sizeof(fingerprint) * num_patterns);
            for (j = 0; j < num_patterns; j++) {
                suffixes[i][j] = init_fingerprint();
            }
        }
        int num_prefixes = 0;
        int num_progressions = 0;
        i = 0;
        while ((1 << (i + 1)) < m_max) {
            matcher->rows[i].row_size = 1 << i;
            lookup_size = 0;
            for (j = 0; j < num_patterns; j++) {
                if ((periods[j] > num_patterns) && (m[j] > num_patterns << 1) && (m[j] - num_patterns > matcher->rows[i].row_size << 1)) {
                    set_fingerprint(matcher->printer, &P[j][matcher->rows[i].row_size], matcher->rows[i].row_size, matcher->tmp);
                    fingerprint_concat(matcher->printer, old_patterns[prev_row[j]], matcher->tmp, patterns[lookup_size]);
                    prev_row[j] = lookup_size;
                    if (m[j] - num_patterns <= (matcher->rows[i].row_size << 2)) {
                        end_pattern[lookup_size] = num_progressions;
                        set_fingerprint(matcher->printer, &P[j][matcher->rows[i].row_size << 1], m[j] - num_patterns - (matcher->rows[i].row_size << 1), matcher->tmp);
                        fingerprint_concat(matcher->printer, patterns[lookup_size], matcher->tmp, prefix[num_prefixes]);
                        progression_index[num_prefixes] = num_progressions;
                        prefix_length[num_prefixes] = m[j] - num_patterns;
                        num_suffixes[num_prefixes] = 1;
                        set_fingerprint(matcher->printer, &P[j][m[j] - num_patterns], num_patterns, suffixes[num_prefixes][0]);
                        for (k = 0; k < num_prefixes; k++) {
                            if (fingerprint_equals(prefix[k], prefix[num_prefixes])) {
                                end_pattern[lookup_size] = progression_index[k];
                                for (l = 0; l < num_suffixes[k]; l++) {
                                    if (fingerprint_equals(suffixes[k][l], suffixes[num_prefixes][0])) {
                                        break;
                                    }
                                }
                                if (l == num_suffixes[k]) {
                                    fingerprint_assign(suffixes[num_prefixes][0], suffixes[k][l]);
                                    num_suffixes[k]++;
                                }
                                break;
                            }
                        }
                        if (k == num_prefixes) num_prefixes++;
                    } else {
                        end_pattern[lookup_size] = -1;
                    }
                    for (k = 0; k < lookup_size; k++) {
                        if (fingerprint_equals(patterns[k], patterns[lookup_size])) {
                            prev_row[j] = k;
                            if ((end_pattern[k] == -1) && (end_pattern[lookup_size] != -1)) {
                                end_pattern[k] = end_pattern[lookup_size];
                                num_progressions++;
                            }
                            else if (end_pattern[lookup_size] != -1) progression_index[num_prefixes - 1] = end_pattern[k];
                            break;
                        }
                    }
                    if (k == lookup_size) {
                        lookup_size++;
                        if (end_pattern[lookup_size -1] != -1) num_progressions++;
                    }
                }
            }
            for (j = 0; j < lookup_size; j++) {
                fingerprint_assign(patterns[j], old_patterns[j]);
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
            matcher->rows[i].next_progression = rbtree_create();
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

        matcher->final.num_prefixes = num_prefixes;
        matcher->final.num_progressions = num_progressions;
        matcher->final.cur_prefix = 0;
        matcher->final.first_location = malloc(sizeof(int) * num_progressions);
        matcher->final.last_location = malloc(sizeof(int) * num_progressions);
        matcher->final.period = malloc(sizeof(int) * num_progressions);
        matcher->final.count = malloc(sizeof(int) * num_progressions);
        matcher->final.first_print = malloc(sizeof(fingerprint) * num_progressions);
        matcher->final.last_print = malloc(sizeof(fingerprint) * num_progressions);
        matcher->final.period_f = malloc(sizeof(fingerprint) * num_progressions);
        matcher->final.progression_index = malloc(sizeof(int) * num_prefixes);
        matcher->final.prefix_length = malloc(sizeof(int) * num_prefixes);
        matcher->final.num_suffixes = malloc(sizeof(int) * num_prefixes);
        matcher->final.suffixes = malloc(sizeof(fingerprint*) * num_prefixes);
        matcher->final.prefix = malloc(sizeof(fingerprint) * num_prefixes);
        matcher->final.suffix_tree = malloc(sizeof(rbtree) * num_patterns);
        for (i = 0; i < num_patterns; i++) {
            matcher->final.suffix_tree[i] = rbtree_create();
        }
        for (i = 0; i < num_progressions; i++) {
            matcher->final.count[i] = 0;
            matcher->final.first_print[i] = init_fingerprint();
            matcher->final.last_print[i] = init_fingerprint();
            matcher->final.period_f[i] = init_fingerprint();
        }
        for (i = 0; i < num_prefixes; i++) {
            matcher->final.progression_index[i] = progression_index[i];
            matcher->final.prefix_length[i] = prefix_length[i];
            matcher->final.num_suffixes[i] = num_suffixes[i];
            matcher->final.suffixes[i] = malloc(sizeof(fingerprint) * num_suffixes[i]);
            for (j = 0; j < num_suffixes[i]; j++) {
                matcher->final.suffixes[i][j] = init_fingerprint();
                fingerprint_assign(suffixes[i][j], matcher->final.suffixes[i][j]);
            }
            matcher->final.prefix[i] = init_fingerprint();
            fingerprint_assign(prefix[i], matcher->final.prefix[i]);
        }

        for (i = 0; i < num_patterns; i++) {
            fingerprint_free(patterns[i]);
            fingerprint_free(old_patterns[i]);
        }
        free(patterns);
        free(old_patterns);
        free(prev_row);
        free(end_pattern);
        free(progression_index);
        free(prefix_length);
        for (i = 0; i < num_patterns; i++) {
            fingerprint_free(prefix[i]);
            for (j = 0; j < num_patterns; j++) {
                fingerprint_free(suffixes[i][j]);
            }
            free(suffixes[i]);
        }
        free(suffixes);
        free(num_suffixes);
        free(prefix);
    }

    matcher->short_matcher = short_dict_matching_build(num_patterns, matcher->printer, P, m);
    matcher->periodic_matcher = periodic_dict_matching_build(P, m, periods, num_patterns, matcher->printer);

    matcher->T_prev = malloc(sizeof(fingerprint) * (num_patterns << 1));
    for (i = 0; i < (num_patterns << 1); i++) {
        matcher->T_prev[i] = init_fingerprint();
    }

    return matcher;
}

int dict_matching_stream(dict_matcher matcher, char T_j, int j) {
    set_fingerprint(matcher->printer, &T_j, 1, matcher->T_j);
    fingerprint_concat(matcher->printer, matcher->T_f, matcher->T_j, matcher->tmp);
    fingerprint_assign(matcher->tmp, matcher->T_f);

    int result = -1, short_result, periodic_result;

    if (matcher->num_rows) {
        int i;
        int cur_prefix = matcher->final.cur_prefix;
        int cur_progression = matcher->final.progression_index[cur_prefix];
        if (matcher->final.count[cur_progression]) {
            int test_location = matcher->final.first_location[cur_progression] + matcher->final.prefix_length[cur_prefix];
            if (test_location + matcher->final.num_prefixes < j) {
                if (matcher->final.count[cur_progression] > 1) {
                    fingerprint_concat(matcher->printer, matcher->final.first_print[cur_progression], matcher->final.period_f[cur_progression], matcher->tmp);
                    fingerprint_assign(matcher->tmp, matcher->final.first_print[cur_progression]);
                    matcher->final.first_location[cur_progression] += matcher->final.period[cur_progression];
                }
                matcher->final.count[cur_progression]--;
                test_location = matcher->final.first_location[cur_progression] + matcher->final.prefix_length[cur_prefix];
            }
            if ((matcher->final.count[cur_progression]) && (test_location < j)) {
                fingerprint_suffix(matcher->printer, matcher->T_prev[test_location % (matcher->num_patterns << 1)], matcher->final.first_print[cur_progression], matcher->tmp);
                if (fingerprint_equals(matcher->tmp, matcher->final.prefix[cur_prefix])) {
                    for (i = 0; i < matcher->final.num_suffixes[cur_prefix]; i++) {
                        rbtree_insert(matcher->final.suffix_tree[test_location % matcher->num_patterns], matcher->final.suffixes[cur_prefix][i], (void*)1, compare_finger);
                    }
                }
            }
        }
        if (matcher->final.num_prefixes > 1) {
            if (++cur_prefix == matcher->final.num_prefixes) cur_prefix = 0;
            cur_progression = matcher->final.progression_index[cur_prefix];
            if (matcher->final.count[cur_progression]) {
                int test_location = matcher->final.first_location[cur_progression] + matcher->final.prefix_length[cur_prefix];
                if (test_location + matcher->final.num_prefixes < j) {
                    if (matcher->final.count[cur_progression] > 1) {
                        fingerprint_concat(matcher->printer, matcher->final.first_print[cur_progression], matcher->final.period_f[cur_progression], matcher->tmp);
                        fingerprint_assign(matcher->tmp, matcher->final.first_print[cur_progression]);
                        matcher->final.first_location[cur_progression] += matcher->final.period[cur_progression];
                    }
                    matcher->final.count[cur_progression]--;
                    test_location = matcher->final.first_location[cur_progression] + matcher->final.prefix_length[cur_prefix];
                }
                if ((matcher->final.count[cur_progression]) && (test_location < j)) {
                    fingerprint_suffix(matcher->printer, matcher->T_prev[test_location % (matcher->num_patterns << 1)], matcher->final.first_print[cur_progression], matcher->tmp);
                    if (fingerprint_equals(matcher->tmp, matcher->final.prefix[cur_prefix])) {
                        for (i = 0; i < matcher->final.num_suffixes[cur_prefix]; i++) {
                            rbtree_insert(matcher->final.suffix_tree[test_location % matcher->num_patterns], matcher->final.suffixes[cur_prefix][i], (void*)1, compare_finger);
                        }
                    }
                }
            }
            matcher->final.cur_prefix = (++cur_prefix == matcher->final.num_prefixes) ? 0 : cur_prefix;
        }

        fingerprint_suffix(matcher->printer, matcher->T_f, matcher->T_prev[(j + matcher->num_patterns) % (matcher->num_patterns << 1)], matcher->tmp);
        int suffix_match = (int)rbtree_lookup(matcher->final.suffix_tree[j % matcher->num_patterns], matcher->tmp, (void*)0, compare_finger);
        if (suffix_match) {
            result = j;
        }
        rbtree_destroy(matcher->final.suffix_tree[j % matcher->num_patterns]);
        matcher->final.suffix_tree[j % matcher->num_patterns] = rbtree_create();

        int occurance, match;
        for (i = matcher->num_rows - 1; i >= 0; i--) {
            occurance = shift_row(matcher->printer, &matcher->rows[i], matcher->current, j, matcher->tmp);
            if (occurance) {
                fingerprint_suffix(matcher->printer, matcher->T_f, matcher->current, matcher->tmp);
                int final = 0;
                match = hashlookup_search(matcher->rows[i].lookup, matcher->tmp, &final);
                if (match != -1) {
                    if (i < matcher->num_rows - 1) {
                        add_occurance(matcher->printer, &matcher->rows[i + 1], matcher->current, j, match, matcher->tmp);
                    }
                    if (final != -1) {
                        int location = j - (matcher->rows[i].row_size << 1);
                        if (matcher->final.count[final]) {
                            if (matcher->final.count[final] == 1) {
                                matcher->final.period[final] = location - matcher->final.first_location[final];
                                fingerprint_suffix(matcher->printer, matcher->current, matcher->final.first_print[final], matcher->final.period_f[final]);
                                matcher->final.last_location[final] = location;
                                fingerprint_assign(matcher->current, matcher->final.last_print[final]);
                                matcher->final.count[final] = 2;
                            } else {
                                fingerprint_suffix(matcher->printer, matcher->current, matcher->final.last_print[final], matcher->tmp);
                                int period = location - matcher->final.last_location[final];
                                if (period >= matcher->num_patterns) {
                                    if ((period == matcher->final.period[final]) && (fingerprint_equals(matcher->tmp, matcher->final.period_f[final]))) {
                                        matcher->final.last_location[final] = location;
                                        fingerprint_assign(matcher->current, matcher->final.last_print[final]);
                                        matcher->final.count[final]++;
                                    } else {
                                        fprintf(stderr, "Warning: Non-Periodic occurance at %d in index %d. Occurance ignored.\n", location, final);
                                    }
                                } else {
                                    fprintf(stderr, "Warning: Prefix %d has period smaller than number of patterns. Occurance at %d ignored.\n", final, location);
                                }
                            }
                        } else {
                            fingerprint_assign(matcher->current, matcher->final.first_print[final]);
                            matcher->final.first_location[final] = location;
                            matcher->final.count[final] = 1;
                        }
                    }
                }
            }
        }

        int first_round = firstlookup_search(matcher->first_round, T_j);
        int index = j % (matcher->num_patterns << 1);
        if (first_round != -1) {
            int prev_index = (index) ? index - 1 : (matcher->num_patterns << 1) - 1;
            add_occurance(matcher->printer, &matcher->rows[0], matcher->T_prev[prev_index], j, first_round, matcher->tmp);
        }
    }

    short_result = short_dict_matching_stream(matcher->short_matcher, matcher->printer, matcher->T_f, matcher->T_prev, matcher->tmp, j);
    periodic_result = periodic_dict_matching_stream(matcher->periodic_matcher, matcher->printer, matcher->T_f, matcher->T_prev, matcher->tmp, j);

    fingerprint_assign(matcher->T_f, matcher->T_prev[j % (matcher->num_patterns << 1)]);
    return ((result != -1) || (short_result != -1) || (periodic_result != -1)) ? j : -1;
}

void dict_matching_free(dict_matcher matcher) {
    fingerprint_free(matcher->T_j);
    fingerprint_free(matcher->tmp);
    fingerprinter_free(matcher->printer);
    fingerprint_free(matcher->T_f);

    int i;

    if (matcher->num_rows){
        fingerprint_free(matcher->current);
        firstlookup_free(&matcher->first_round);
        int j;
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

            hashlookup_free(&matcher->rows[i].lookup);
        }
        free(matcher->rows);

        free(matcher->final.first_location);
        free(matcher->final.last_location);
        free(matcher->final.period);
        free(matcher->final.progression_index);
        free(matcher->final.prefix_length);
        for (i = 0; i < matcher->final.num_progressions; i++) {
            fingerprint_free(matcher->final.first_print[i]);
            fingerprint_free(matcher->final.last_print[i]);
            fingerprint_free(matcher->final.period_f[i]);
        }
        free(matcher->final.first_print);
        free(matcher->final.last_print);
        free(matcher->final.period_f);
        for (i = 0; i < matcher->final.num_prefixes; i++) {
            for (j = 0; j < matcher->final.num_suffixes[i]; j++) {
                fingerprint_free(matcher->final.suffixes[i][j]);
            }
            free(matcher->final.suffixes[i]);
            fingerprint_free(matcher->final.prefix[i]);
        }
        free(matcher->final.num_suffixes);
        free(matcher->final.suffixes);
        free(matcher->final.prefix);
        for (i = 0; i < matcher->num_patterns; i++) {
            rbtree_destroy(matcher->final.suffix_tree[i]);
        }
        free(matcher->final.suffix_tree);
    }

    short_dict_matching_free(matcher->short_matcher);
    periodic_dict_matching_free(matcher->periodic_matcher);

    for (i = 0; i < (matcher->num_patterns << 1); i++) {
        fingerprint_free(matcher->T_prev[i]);
    }
    free(matcher->T_prev);

    free(matcher);
}

#endif
