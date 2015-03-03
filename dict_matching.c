#include <assert.h>
#include <string.h>
#include "dict_matching.h"

char **P;
int *m, num_patterns, n, alpha, *correct;

void build_up() {
    int i;
    P = malloc(sizeof(char*) * num_patterns);
    for (i = 0; i < num_patterns; i++) {
        P[i] = malloc(sizeof(char) * (m[i] + 1));
    }
    correct = malloc(sizeof(int) * n);
}

void stream_test(char *T) {
    dict_matcher matcher = dict_matching_build(P, m, num_patterns, 100, alpha);
    int i;
    for (i = 0; i < n; i++) {
        assert(dict_matching_stream(matcher, T[i], i) == correct[i]);
    }
    dict_matching_free(matcher);
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

void test_single_pattern() {
    m = malloc(sizeof(int));
    m[0] = 8;
    num_patterns = 1;
    n = 30;
    alpha = 0;
    build_up();
    strcpy(P[0], "ababbabb");
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  = -1; correct[4]  = -1; correct[5]  = -1;
    correct[6]  = -1; correct[7]  = -1; correct[8]  = -1; correct[9]  = -1; correct[10] = -1; correct[11] = -1;
    correct[12] = -1; correct[13] = -1; correct[14] = -1; correct[15] = -1; correct[16] = -1; correct[17] = -1;
    correct[18] = -1; correct[19] = -1; correct[20] = -1; correct[21] = -1; correct[22] = -1; correct[23] = -1;
    correct[24] = -1; correct[25] = -1; correct[26] = -1; correct[27] = -1; correct[28] = -1; correct[29] = 29;
    stream_test("ababbabcababbbbbababbaababbabb");
    tear_down();
}

void test_two_patterns() {
    m = malloc(sizeof(int) * 2);
    m[0] = 8;
    m[1] = 8;
    num_patterns = 2;
    n = 30;
    alpha = 0;
    build_up();
    strcpy(P[0], "ababbabb");
    strcpy(P[1], "ababbabc");
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  = -1; correct[4]  = -1; correct[5]  = -1;
    correct[6]  = -1; correct[7]  =  7; correct[8]  = -1; correct[9]  = -1; correct[10] = -1; correct[11] = -1;
    correct[12] = -1; correct[13] = -1; correct[14] = -1; correct[15] = -1; correct[16] = -1; correct[17] = -1;
    correct[18] = -1; correct[19] = -1; correct[20] = -1; correct[21] = -1; correct[22] = -1; correct[23] = -1;
    correct[24] = -1; correct[25] = -1; correct[26] = -1; correct[27] = -1; correct[28] = -1; correct[29] = 29;
    stream_test("ababbabcababbbbbababbaababbabb");
    tear_down();
}

void test_different_lengths() {
    n = 46;
    num_patterns = 5;
    m = malloc(sizeof(int) * num_patterns);
    m[0] = 4; m[1] = 4; m[2] = 4; m[3] = 2; m[4] = 1;
    build_up();
    strcpy(P[0], "take");
    strcpy(P[1], "fast");
    strcpy(P[2], "sofa");
    strcpy(P[3], "so");
    strcpy(P[4], "s");
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  =  3; correct[4]  =  4; correct[5]  =  5;
    correct[6]  = -1; correct[7]  = -1; correct[8]  = -1; correct[9]  =  9; correct[10] = 10; correct[11] = -1;
    correct[12] = 12; correct[13] = 13; correct[14] = 14; correct[15] = -1; correct[16] = -1; correct[17] = -1;
    correct[18] = 18; correct[19] = 19; correct[20] = 20; correct[21] = -1; correct[22] = 22; correct[23] = -1;
    correct[24] = -1; correct[25] = -1; correct[26] = 26; correct[27] = -1; correct[28] = 28; correct[29] = 29;
    correct[30] = 30; correct[31] = 31; correct[32] = 32; correct[33] = -1; correct[34] = 34; correct[35] = 35;
    correct[36] = -1; correct[37] = 37; correct[38] = 38; correct[39] = 39; correct[40] = -1; correct[41] = -1;
    correct[42] = 42; correct[43] = -1; correct[44] = 44; correct[45] = 45;
    stream_test("takeso fasofast fassofatake sosso sofastake so");
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
    stream_test("abccab");
    tear_down();
}

void test_long_and_short() {
    m = malloc(sizeof(int) * 3);
    m[0] = 8;
    m[1] = 8;
    m[2] = 3;
    num_patterns = 3;
    n = 30;
    alpha = 0;
    build_up();
    strcpy(P[0], "ababbabb");
    strcpy(P[1], "ababbabc");
    strcpy(P[2], "abb");
    correct[0]  = -1; correct[1]  = -1; correct[2]  = -1; correct[3]  = -1; correct[4]  =  4; correct[5]  = -1;
    correct[6]  = -1; correct[7]  =  7; correct[8]  = -1; correct[9]  = -1; correct[10] = -1; correct[11] = -1;
    correct[12] = 12; correct[13] = -1; correct[14] = -1; correct[15] = -1; correct[16] = -1; correct[17] = -1;
    correct[18] = -1; correct[19] = -1; correct[20] = 20; correct[21] = -1; correct[22] = -1; correct[23] = -1;
    correct[24] = -1; correct[25] = -1; correct[26] = 26; correct[27] = -1; correct[28] = -1; correct[29] = 29;
    stream_test("ababbabcababbbbbababbaababbabb");
    tear_down();
}

void test_all_patterns() {
    n = 7;
    num_patterns = 7;
    m = malloc(sizeof(int) * num_patterns);
    m[0] = 1; m[1] = 2; m[2] = 3; m[3] = 4; m[4] = 5; m[5] = 6; m[6] = 7;
    build_up();
    strcpy(P[0], "a");
    strcpy(P[1], "ab");
    strcpy(P[2], "abc");
    strcpy(P[3], "abcd");
    strcpy(P[4], "abcde");
    strcpy(P[5], "abcdef");
    strcpy(P[6], "abcdefg");
    correct[0]  = 0; correct[1]  = 1; correct[2]  = 2; correct[3]  = 3; correct[4]  = 4; correct[5]  = 5; correct[6]  = 6;
    stream_test("abcdefg");
    tear_down();
}

int main(void) {
    test_single_pattern();
    test_two_patterns();
    test_different_lengths();
    test_seven_patterns();
    test_long_and_short();
    test_all_patterns();
    return 0;
}