#ifndef __FIRST_LOOKUP__
#define __FIRST_LOOKUP__

#include <cmph.h>

typedef struct {
    cmph_t *hash;
    int num;
    char *keys;
} first_lookup;

first_lookup firstlookup_build(char *key_string, int num) {
    first_lookup lookup;
    lookup.num = num;

    if (num > 1) {
        char **keys = malloc(sizeof(char*) * num);
        int i;
        for (i = 0; i < num; i++) {
            keys[i] = malloc(sizeof(char) * 2);
            keys[i][0] = key_string[i];
            keys[i][1] = '\0';
        }
        cmph_io_adapter_t *source = cmph_io_vector_adapter(keys, num);
        cmph_config_t *config = cmph_config_new(source);
        cmph_config_set_algo(config, CMPH_CHD);
        lookup.hash = cmph_new(config);
        cmph_config_destroy(config);
        cmph_io_vector_adapter_destroy(source);
        lookup.keys = malloc(num * sizeof(char));

        unsigned int id;
        char *key;
        for (i = 0; i < num; i++) {
            key = keys[i];
            id = cmph_search(lookup.hash, key, 1);
            lookup.keys[id] = key_string[i];
            free(keys[i]);
        }
        free(keys);
    } else if (num == 1) {
        lookup.keys = malloc(sizeof(char));
        lookup.keys[0] = key_string[0];
    }

    return lookup;
}

int firstlookup_search(first_lookup lookup, char key) {
    if (lookup.num == 0) return -1;
    if (lookup.num == 1) {
        if (key == lookup.keys[0]) {
            return 0;
        }
        return -1;
    }
    int id = cmph_search(lookup.hash, &key, 1);
    if ((id < lookup.num) && (key == lookup.keys[id])) {
        return id;
    }
    return -1;
}

void firstlookup_free(first_lookup *lookup) {
    if (lookup->num > 0) {
        if (lookup->num > 1) cmph_destroy(lookup->hash);
        free(lookup->keys);
    }
}

int firstlookup_size(first_lookup lookup) {
    int size = sizeof(char) * lookup.num;
    if (lookup.num > 1) {
        size += cmph_packed_size(lookup.hash);
    }
    return size;
}

#endif