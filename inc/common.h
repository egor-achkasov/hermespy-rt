#define IN
#define OUT

#define FREE_POINTERS(...) \
    do {                            \
        void *ptrs[] = {__VA_ARGS__};\
        for (size_t i = 0; i < sizeof(ptrs) / sizeof(ptrs[0]); i++) \
            free(ptrs[i]);          \
    } while (0)

#define PERROR_CLEANUP_EXIT(msg, rc, ...) \
    do {                            \
        perror(msg);                \
        FREE_POINTERS(__VA_ARGS__); \
        exit(rc);                    \
    } while (0)
