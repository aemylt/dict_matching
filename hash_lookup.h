/*
    hash_lookup.h
    A dictionary for storing key-value pairs.
    Utilises the C Minimum Perfect Hashing library (http://cmph.sourceforge.net/)
*/

#ifndef HASH_LOOKUP
#define HASH_LOOKUP

#include <cmph.h>
#include <gmp.h>
#include <stdlib.h>
#include "karp_rabin.h"

/*
    typedef struct hash_lookup
    Structure for holding the pairs.
    Components:
        cmph_t *hash   - The hash function
        int    *values - The values
        char   *keys   - The keys
        int    num     - The number of items
*/
typedef struct {
    cmph_t *hash;
    int num;
    char *last_key;
    int key_size;
    fingerprint *keys;
    int *end_pattern;
} hash_lookup;

/*
    hashlookup_build
    Constructs a hash_lookup object.
    Components:
        char **keys - A list of the keys as strings 1 character long
        int *values - A list of the values for each key
        int num - The number of key-value pairs
    Returns hash_lookup:
        The constructed dictionary
*/
hash_lookup hashlookup_build(fingerprint *prints, int *end_pattern, int num, fingerprinter printer) {
    hash_lookup lookup;
    lookup.num = num;

    if (num > 1) {
        int size = printer->p->_mp_size * 8 + 1;
        char **keys = malloc(sizeof(char*) * num);
        int i;
        int *sizes = malloc(sizeof(int) * num);
        for (i = 0; i < num; i++) {
            keys[i] = malloc(sizeof(char) * size);
            sizes[i] = gmp_snprintf(keys[i], size, "%Zx", prints[i]->finger);
        }
        cmph_io_adapter_t *source = cmph_io_vector_adapter(keys, num);
        cmph_config_t *config = cmph_config_new(source);
        cmph_config_set_algo(config, CMPH_CHD);
        lookup.hash = cmph_new(config);
        cmph_config_destroy(config);
        cmph_io_vector_adapter_destroy(source);

        lookup.keys = malloc(sizeof(fingerprint) * num);
        lookup.end_pattern = malloc(sizeof(int));
        int location;
        for (i = 0; i < num; i++) {
            location = cmph_search(lookup.hash, keys[i], sizes[i]);
            lookup.keys[location] = init_fingerprint();
            fingerprint_assign(prints[i], lookup.keys[location]);
            free(keys[i]);
            lookup.end_pattern[location] = end_pattern[i];
        }
        free(keys);
        lookup.key_size = size;
        lookup.last_key = malloc(sizeof(char) * size);
    } else if (num == 1) {
        lookup.keys = malloc(sizeof(fingerprint));
        lookup.keys[0] = init_fingerprint();
        fingerprint_assign(prints[0], lookup.keys[0]);
        lookup.end_pattern = malloc(sizeof(int));
        lookup.end_pattern[0] = end_pattern[0];
    }

    return lookup;
}

/*
    hashlookup_search
    Searches the dictionary for the corresponding value to a key.
    Parameters:
        hash_lookup lookup - The dictionary to search
        char        key    - The key to search for
    Returns int:
        values[lookup.hash(key)] if key \in keys
        0 otherwise
*/
int hashlookup_search(hash_lookup lookup, fingerprint key, int *match) {
    if (lookup.num == 0) return -1;
    if (lookup.num == 1) {
        if (fingerprint_equals(key, lookup.keys[0])) {
            *match |= lookup.end_pattern[0];
            return 0;
        }
        return -1;
    }
    int size = gmp_snprintf(lookup.last_key, lookup.key_size, "%Zx", key->finger);
    int location =  cmph_search(lookup.hash, lookup.last_key, size);
    if ((location < lookup.num) && fingerprint_equals(lookup.keys[location], key)) {
        *match |= lookup.end_pattern[location];
        return location;
    }
    return -1;
}

/*
    hashlookup_free
    Frees the dictionary.
    Parameters:
        hash_lookup *lookup - The dictionary to free
*/
void hashlookup_free(hash_lookup *lookup) {
    if (lookup->num > 1) {
        cmph_destroy(lookup->hash);
        free(lookup->last_key);
        int i;
        for (i = 0; i < lookup->num; i++) {
            fingerprint_free(lookup->keys[i]);
        }
    } else if (lookup->num == 1) {
        fingerprint_free(lookup->keys[0]);
    }
    free(lookup->keys);
    free(lookup->end_pattern);
}

#endif