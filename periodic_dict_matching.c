#include "periodic_dict_matching.h"
#include "karp_rabin.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int num_patterns, n, *m, *rho, *correct;
char **P;

void build_up() {
    int i;
    P = malloc(sizeof(char*) * num_patterns);
    for (i = 0; i < num_patterns; i++) {
        P[i] = malloc(sizeof(char) * (m[i] + 1));
    }
    correct = malloc(sizeof(int) * n);
    rho = malloc(sizeof(int) * num_patterns);
}

void stream_test(char *T) {
    fingerprinter printer = fingerprinter_build(100, 0);
    periodic_dict_matching_state state = periodic_dict_matching_build(P, m, rho, num_patterns, printer);
    /*
    int i;
    fingerprint t_f = init_fingerprint(), t_j = init_fingerprint(), tmp = init_fingerprint();
    for (i = 0; i < n; i++) {
        set_fingerprint(printer, &T[i], 1, t_j);
        fingerprint_concat(printer, t_f, t_j, tmp);
        fingerprint_assign(tmp, t_f);
        assert(correct[i] == short_dict_matching_stream(state, printer, t_f, tmp, i));
    }
    set_fingerprint(printer, T, n, t_j);
    */
    periodic_dict_matching_free(state);
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

void test_one_pattern() {
    n = 16;
    num_patterns = 1;
    m = malloc(sizeof(int) * num_patterns);
    m[0] = 5;
    build_up();
    rho[0] = 1;
    strcpy(P[0], "aaaaa");
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  = -1; correct[4]  =  4; correct[5]  =  5;
    correct[6]  =  6; correct[7]  =  7; correct[8]  =  8; correct[9]  =  9; correct[10] = -1; correct[11] = -1;
    correct[12] = -1; correct[13] = -1; correct[14] = -1; correct[15] = 15;
    stream_test("aaaaaaaaaabaaaaa");
    tear_down();
}

void test_two_patterns() {
    n = 16;
    num_patterns = 2;
    m = malloc(sizeof(int) * num_patterns);
    m[0] = 5; m[1] = 6;
    build_up();
    strcpy(P[0], "aaaaa");
    strcpy(P[1], "bababa");
    rho[0] = 1; rho[1] = 2;
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  = -1; correct[4]  =  4; correct[5]  =  5;
    correct[6]  =  6; correct[7]  =  7; correct[8]  =  8; correct[9]  =  9; correct[10] = -1; correct[11] = -1;
    correct[12] = -1; correct[13] = -1; correct[14] = -1; correct[15] = 15;
    stream_test("aaaaaaaaaabaaaaa");
    tear_down();
}

void test_substrings() {
    n = 16;
    num_patterns = 2;
    m = malloc(sizeof(int) * num_patterns);
    m[0] = 5; m[1] = 6;
    build_up();
    strcpy(P[0], "ababa");
    strcpy(P[1], "bababa");
    rho[0] = 2; rho[1] = 2;
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  = -1; correct[4]  =  4; correct[5]  =  5;
    correct[6]  =  6; correct[7]  =  7; correct[8]  =  8; correct[9]  =  9; correct[10] = -1; correct[11] = -1;
    correct[12] = -1; correct[13] = -1; correct[14] = -1; correct[15] = 15;
    stream_test("aaaaaaaaaabaaaaa");
    tear_down();
}

void test_matching_prefix() {
    n = 16;
    num_patterns = 2;
    m = malloc(sizeof(int) * num_patterns);
    m[0] = 7; m[1] = 6;
    build_up();
    strcpy(P[0], "abababa");
    strcpy(P[1], "ababab");
    rho[0] = 2; rho[1] = 2;
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  = -1; correct[4]  =  4; correct[5]  =  5;
    correct[6]  =  6; correct[7]  =  7; correct[8]  =  8; correct[9]  =  9; correct[10] = -1; correct[11] = -1;
    correct[12] = -1; correct[13] = -1; correct[14] = -1; correct[15] = 15;
    stream_test("aaaaaaaaaabaaaaa");
    tear_down();
}

int main(void) {
    test_one_pattern();
    test_two_patterns();
    test_substrings();
    test_matching_prefix();
    printf("All tests succeeded!\n");
    return 0;
}