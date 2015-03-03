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

void test_four_patterns() {
    n = 46;
    num_patterns = 4;
    m = malloc(sizeof(int) * num_patterns);
    m[0] = 4; m[1] = 4; m[2] = 4; m[3] = 2;
    build_up();
    strcpy(P[0], "take");
    strcpy(P[1], "fast");
    strcpy(P[2], "sofa");
    strcpy(P[3], "so");
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  = 3; correct[4]  = -1; correct[5]  = 5;
    correct[6]  = -1; correct[7]  = -1; correct[8]  = -1; correct[9]  = -1; correct[10] = 10; correct[11] = -1;
    correct[12] = 12; correct[13] = -1; correct[14] = 14; correct[15] = -1; correct[16] = -1; correct[17] = -1;
    correct[18] = -1; correct[19] = -1; correct[20] = 20; correct[21] = -1; correct[22] = 22; correct[23] = -1;
    correct[24] = -1; correct[25] = -1; correct[26] = 26; correct[27] = -1; correct[28] = -1; correct[29] = 29;
    correct[30] = -1; correct[31] = -1; correct[32] = 32; correct[33] = -1; correct[34] = -1; correct[35] = 35;
    correct[36] = -1; correct[37] = 37; correct[38] = -1; correct[39] = 39; correct[40] = -1; correct[41] = -1;
    correct[42] = 42; correct[43] = -1; correct[44] = -1; correct[45] = 45;
    stream_test("takeso fasofast fassofatake sosso sofastake so", n, P, m, num_patterns, correct);
    tear_down();
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
    test_four_patterns();
    test_seven_patterns();
    printf("All tests succeeded!\n");
    return 0;
}