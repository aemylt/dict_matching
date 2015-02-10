#ifndef __DICT_MATCHING__
#define __DICT_MATCHING__

#include <stdlib.h>
#include "karp_rabin.h"
#include "hash_lookup.h"

typedef struct dict_matcher_t {
    fingerprinter printer;
    fingerprint text_f;
    fingerprint T_j;
    fingerprint tmp;
    int num_patterns;
    char *text;
    int m;
    hash_lookup lookup;
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
    free(matcher);
}

#endif