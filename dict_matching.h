#ifndef __DICT_MATCHING__
#define __DICT_MATCHING__

#include <stdlib.h>
#include "karp_rabin.h"
#include "hash_lookup.h"
#include "first_lookup.h"

typedef struct {
    int location;
    fingerprint T_f;
} viable_occurance;

typedef struct {
    int row_size, num_patterns, *start, *end;
    hash_lookup lookup;
    viable_occurance **VOs;
} pattern_row;

int shift_row(pattern_row row, fingerprint finger, int j) {
    int i;
    for (i = 0; i < row.num_patterns; i++) {
        if ((row.start[i] != row.end[i]) && (row.VOs[i][row.start[i]].location == j)) {
            fingerprint_assign(row.VOs[i][row.start[i]].T_f, finger);
            if (++row.start[i] == row.row_size) row.start[i] = 0;
            return 1;
        }
    }
    return 0;
}

void add_occurance(pattern_row row, fingerprint finger, int location, int i) {
    fingerprint_assign(finger, row.VOs[i][row.end[i]].T_f);
    row.VOs[i][row.end[i]].location = i + row.row_size;
    if (++row.end[i] == row.row_size) row.end[i] = 0;
}

typedef struct dict_matcher_t {
    fingerprinter printer;
    fingerprint text_f;
    fingerprint T_j;
    fingerprint tmp;
    int num_patterns;
    char *text;
    int m;
    hash_lookup lookup;

    pattern_row *rows;
    int num_rows;
    first_lookup first_round;
} *dict_matcher;

dict_matcher dict_matching_build(char **P, int *m, int num_patterns, int n, int alpha) {
    dict_matcher matcher = malloc(sizeof(struct dict_matcher_t));
    matcher->num_patterns = num_patterns;
    matcher->printer = fingerprinter_build(n, alpha);
    fingerprint *patterns = malloc(sizeof(fingerprint) * num_patterns);
    int i;
    matcher->m = 0;
    for (i = 0; i < num_patterns; i++) {
        patterns[i] = init_fingerprint();
        set_fingerprint(matcher->printer, P[i], m[i], patterns[i]);
    }
    matcher->m = m[0];
    matcher->text = malloc(sizeof(char) * matcher->m);
    matcher->text_f = init_fingerprint();
    matcher->T_j = init_fingerprint();
    matcher->tmp = init_fingerprint();
    matcher->lookup = hashlookup_build(patterns, num_patterns, matcher->printer);

    int m_max = m[0];
    matcher->num_rows = 0;
    while ((1 << matcher->num_rows) < m_max) {
        matcher->num_rows++;
    }
    matcher->rows = malloc(sizeof(pattern_row) * matcher->num_rows);

    int j, k, lookup_size;
    char *first_letters = malloc(sizeof(char) * num_patterns);
    lookup_size = 0;
    for (j = 0; j < num_patterns; j++) {
        first_letters[lookup_size] = P[j][0];
        for (k = 0; k < lookup_size; k++) {
            if (first_letters[k] == first_letters[lookup_size]) break;
        }
        if (k == lookup_size) lookup_size++;
    }
    matcher->first_round = firstlookup_build(first_letters, lookup_size);
    free(first_letters);

    i = 1;
    lookup_size = 0;
    while ((1 << i) <= m_max) {
        matcher->rows[i - 1].row_size = 1 << i;
        lookup_size = 0;
        for (j = 0; j < num_patterns; j++) {
            set_fingerprint(matcher->printer, P[j], matcher->rows[i - 1].row_size << 1, patterns[lookup_size]);
            for (k = 0; k < lookup_size; k++) {
                if (fingerprint_equals(patterns[k], patterns[lookup_size])) break;
            }
            if (k == lookup_size) lookup_size++;
        }
        matcher->rows[i - 1].lookup = hashlookup_build(patterns, lookup_size, matcher->printer);
        matcher->rows[i - 1].num_patterns = lookup_size;
        matcher->rows[i - 1].start = malloc(sizeof(int) * lookup_size);
        matcher->rows[i - 1].end = malloc(sizeof(int) * lookup_size);
        matcher->rows[i - 1].VOs = malloc(sizeof(viable_occurance*) * lookup_size);
        for (j = 0; j < lookup_size; j++) {
            matcher->rows[i - 1].start[j] = 0;
            matcher->rows[i - 1].end[j] = 0;
            matcher->rows[i - 1].VOs[j] = malloc(sizeof(viable_occurance) * matcher->rows[i - 1].row_size);
            for (k = 0; k < matcher->rows[i - 1].row_size; k++) {
                matcher->rows[i - 1].VOs[j][k].T_f = init_fingerprint();
            }
        }
        i++;
    }

    for (i = 0; i < num_patterns; i++) {
        fingerprint_free(patterns[i]);
    }
    free(patterns);

    return matcher;
}

int dict_matching_stream(dict_matcher matcher, char T_j, int j) {
    if (j >= matcher->m) {
        set_fingerprint(matcher->printer, &matcher->text[j % matcher->m], 1, matcher->T_j);
        fingerprint_suffix(matcher->printer, matcher->text_f, matcher->T_j, matcher->tmp);
        fingerprint_assign(matcher->tmp, matcher->text_f);
    }
    matcher->text[j % matcher->m] = T_j;
    set_fingerprint(matcher->printer, &T_j, 1, matcher->T_j);
    fingerprint_concat(matcher->printer, matcher->text_f, matcher->T_j, matcher->tmp);
    fingerprint_assign(matcher->tmp, matcher->text_f);
    if (j >= matcher->m - 1) {
        return (hashlookup_search(matcher->lookup, matcher->text_f) != -1) ? j : -1;
    }
    return -1;
}

void dict_matching_free(dict_matcher matcher) {
    fingerprint_free(matcher->text_f);
    fingerprint_free(matcher->T_j);
    fingerprint_free(matcher->tmp);
    fingerprinter_free(matcher->printer);
    free(matcher->text);
    hashlookup_free(&matcher->lookup);

    firstlookup_free(&matcher->first_round);
    int i, j, k;
    for (i = 0; i < matcher->num_rows; i++) {
        for (j = 0; j < matcher->rows[i].num_patterns; j++) {
            for (k = 0; k < matcher->rows[i].row_size; k++) {
                fingerprint_free(matcher->rows[i].VOs[j][k].T_f);
            }
            free(matcher->rows[i].VOs[j]);
        }
        hashlookup_free(&matcher->rows[i].lookup);
    }
    free(matcher->rows);

    free(matcher);
}

#endif