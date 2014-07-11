#include "rsmemory.h"

void *rsMalloc(size_t size)
{
    void *data;
    data = malloc(size);

    if ( data == NULL ) {
        fprintf(stderr, "OutOfMemoryException: failed to allocate %ld bytes\n", size);
        return NULL;
    }

    return data;
}
