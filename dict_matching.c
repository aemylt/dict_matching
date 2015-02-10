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
    dict_matcher matcher = dict_matching_build(P, m, num_patterns, n, alpha);
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

int main(void) {
    test_single_pattern();
    test_two_patterns();
    return 0;
}