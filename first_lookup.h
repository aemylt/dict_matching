#ifndef __FIRST_LOOKUP__
#define __FIRST_LOOKUP__

typedef struct {
    cmph_t *hash;
    int num;
    char *keys;
    int *end_pattern;
} first_lookup;

first_lookup firstlookup_build(char *key_string, int *end_pattern, int num) {
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
        lookup.end_pattern = malloc(sizeof(int));
        for (i = 0; i < num; i++) {
            key = keys[i];
            id = cmph_search(lookup.hash, key, 1);
            lookup.keys[id] = key_string[i];
            lookup.end_pattern[id] = end_pattern[i];
            free(keys[i]);
        }
        free(keys);
    } else if (num == 1) {
        lookup.keys = malloc(sizeof(char));
        lookup.keys[0] = key_string[0];
        lookup.end_pattern = malloc(sizeof(int));
        lookup.end_pattern[0] = end_pattern[0];
    }

    return lookup;
}

int firstlookup_search(first_lookup lookup, char key, int *match) {
    if (lookup.num == 0) return -1;
    if (lookup.num == 1) {
        if (key == lookup.keys[0]) {
            *match |= lookup.end_pattern[0];
            return 0;
        }
        return -1;
    }
    int id = cmph_search(lookup.hash, &key, 1);
    if ((id < lookup.num) && (key == lookup.keys[id])) {
        *match |= lookup.end_pattern[id];
        return id;
    }
    return -1;
}

void firstlookup_free(first_lookup *lookup) {
    if (lookup->num > 0) {
        if (lookup->num > 1) cmph_destroy(lookup->hash);
        free(lookup->keys);
        free(lookup->end_pattern);
    }
}

#endif