#include <string.h>
#include <stdio.h>
#include "dict_matching.h"

char **P;
int *m, num_patterns, n, alpha;

void build_up() {
    int i;
    P = malloc(sizeof(char*) * num_patterns);
    for (i = 0; i < num_patterns; i++) {
        P[i] = malloc(sizeof(char) * (m[i] + 1));
    }
}

void stream_test(char *T) {
    dict_matcher matcher = dict_matching_build(P, m, num_patterns, n, alpha);
    int i;
    for (i = 0; i < n; i++) {
        printf("%d,", dict_matching_stream(matcher, T[i], i));
    }
    printf("\n");
    dict_matching_free(matcher);
}

void tear_down() {
    int i;
    for (i = 0; i < num_patterns; i++) {
        free(P[i]);
    }
    free(P);
    free(m);
}

void test_single_pattern() {
    m = malloc(sizeof(int));
    m[0] = 8;
    num_patterns = 1;
    n = 30;
    alpha = 0;
    build_up();
    strcpy(P[0], "ababbabb");
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
    stream_test("ababbabcababbbbbababbaababbabb");
    tear_down();
}

int main(void) {
    test_single_pattern();
    test_two_patterns();
    return 0;
}