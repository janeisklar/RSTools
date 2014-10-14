#include "rsio.h"
#include <stdio.h>
#include <unistd.h>

BOOL rsFileIsReadable(const char *path)
{
    FILE *f;
    f = fopen(path, "r");
    if (f != NULL) {
        fclose(f);
        return TRUE;
    }
    return FALSE;
}

BOOL rsFileIsWritable(const char *path)
{
    FILE *f;
    
    // check if the file can be opened for writing without actually creating it
    f = fopen(path, "r+");
    
    if (f != NULL) {
        // worked, let's close it again
        fclose(f);
        return TRUE;
    }
    
    // maybe it doest't exist yet, so let's try to create it
    f = fopen(path, "a+");
    
    if (f != NULL) {
        // that worked, so remove the created file
        fclose(f);
        unlink(path);
        return TRUE;
    }
    
    // we don't seem to have permissions to write or the path doesn't exist
    return FALSE;
}

BOOL rsCheckInputs(const char **paths)
{
    BOOL readable = TRUE;
    int i = 0;
    while ( paths[i] == NULL || strcmp(paths[i], RSIO_LASTFILE) != 0 ) {
        if ( paths[i] != NULL && ! rsFileIsReadable(paths[i]) ) {
            fprintf(stderr, "Error: File '%s' could not be read!\n", paths[i]);
            readable = FALSE;
        }
        i++;
    }
    return readable;
}

BOOL rsCheckOutputs(const char **paths)
{
    BOOL writable = TRUE;
    int i = 0;
    while ( paths[i] == NULL || strcmp(paths[i], RSIO_LASTFILE) != 0 ) {
        if ( paths[i] != NULL && ! rsFileIsWritable(paths[i]) ) {
            fprintf(stderr, "Error: File '%s' is not writable!\n", paths[i]);
            writable = FALSE;
        }
        i++;
    }
    return writable;
}
