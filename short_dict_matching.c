#include "short_dict_matching.h"
#include "karp_rabin.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int num_patterns, n, *m, *correct;
char **P;

void build_up() {
    int i;
    P = malloc(sizeof(char*) * num_patterns);
    for (i = 0; i < num_patterns; i++) {
        P[i] = malloc(sizeof(char) * (m[i] + 1));
    }
    correct = malloc(sizeof(int) * n);
}

void stream_test(char *T, int n, char **P, int *m, int num_patterns, int *correct) {
    fingerprinter printer = fingerprinter_build(100, 0);
    short_dict_matcher state = short_dict_matching_build(num_patterns, printer, P, m);
    int i;
    fingerprint t_f = init_fingerprint(), t_j = init_fingerprint(), tmp = init_fingerprint();
    for (i = 0; i < n; i++) {
        set_fingerprint(printer, &T[i], 1, t_j);
        fingerprint_concat(printer, t_f, t_j, tmp);
        fingerprint_assign(tmp, t_f);
        assert(correct[i] == short_dict_matching_stream(state, printer, t_f, tmp, i));
    }
    set_fingerprint(printer, T, n, t_j);
    short_dict_matching_free(state);
}

void tear_down() {
    int i;
    for (i = 0; i < num_patterns; i++) {
        free(P[i]);
    }
    free(P);
    free(m);
    free(correct);
}

void test_seven_patterns() {
    n = 6;
    num_patterns = 7;
    m = malloc(sizeof(int) * num_patterns);
    m[0] = 1; m[1] = 2; m[2] = 3; m[3] = 2; m[4] = 3; m[5] = 1; m[6] = 3;
    build_up();
    strcpy(P[0], "a");
    strcpy(P[1], "ab");
    strcpy(P[2], "bab");
    strcpy(P[3], "bc");
    strcpy(P[4], "bca");
    strcpy(P[5], "c");
    strcpy(P[6], "caa");
    correct[0]  = 0; correct[1]  = 1; correct[2]  = 2; correct[3]  = 3; correct[4]  = 4; correct[5]  = 5;
    stream_test("abccab", n, P, m, num_patterns, correct);
    tear_down();
}

int main(void) {
    test_seven_patterns();
    printf("All tests succeeded!\n");
    return 0;
}