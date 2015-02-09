#include "hash_lookup.h"
#include <assert.h>

int main(void) {
    char *keys = malloc(10 * sizeof(char*));
    int i;
    keys[0] = 'a'; keys[1] = 'b'; keys[2] = 'X'; keys[3] = 'Y'; keys[4] = '0';
    keys[5] = '1'; keys[6] = '?'; keys[7] = '$'; keys[8] = '~'; keys[9] = '@';

    int *values = malloc(10 * sizeof(int));
    for (i = 0; i < 10; i++) values[i] = i * i;

    hash_lookup lookup = hashlookup_build(keys, values, 10);
    for (i = 0; i < 10; i++) assert(hashlookup_search(lookup, keys[i]) == values[i]);
    assert(hashlookup_search(lookup, 'Z') == 0);
    hashlookup_free(&lookup);

    lookup = hashlookup_build(keys, values, 1);
    assert(hashlookup_search(lookup, 'a') == 0);
    assert(hashlookup_search(lookup, 'Z') == 0);
    hashlookup_free(&lookup);

    lookup = hashlookup_build(keys, values, 0);
    assert(hashlookup_search(lookup, 'a') == 0);
    assert(hashlookup_search(lookup, 'Z') == 0);
    hashlookup_free(&lookup);

    keys[1] = '1';
    lookup = hashlookup_build(keys, values, 2);
    assert(hashlookup_search(lookup, 'a') == 0);
    assert(hashlookup_search(lookup, '1') == 1);
    assert(hashlookup_search(lookup, 'Z') == 0);
    hashlookup_free(&lookup);

    lookup = hashlookup_build(keys, values, 2);
    keys[1] = 'Z';
    assert(hashlookup_search(lookup, 'a') == 0);
    assert(hashlookup_search(lookup, '1') == 1);
    assert(hashlookup_search(lookup, 'Z') == 0);
    hashlookup_free(&lookup);

    return 0;
}