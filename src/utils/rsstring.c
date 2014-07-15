#include "rsstring.h"

/*
 * Takes in a string and breaks it after a specified amount
 * of characters. While doing so it tries not to split up
 * any words.
 * 
 * Result: lineArray - an array of strings (one for each line)
 *         nLines    - the number of lines
 */
void rsStringWordWrap(const char* inputString, char*** lineArray, size_t* nLines, const unsigned int lineLength)
{
    // start with a single line
    size_t wordStart = 0, lastCopied = 0;

    const size_t one         = 1;
    const size_t inputLength = strlen(inputString);
    const size_t lineSize    = sizeof(char)*(lineLength+one);
    
    *lineArray  = (char**)rsMalloc(sizeof(char*)); 
    char*  line = (char*) rsMalloc(lineSize);
    line[0] = '\0';
    (*lineArray)[0] = line;
    *nLines = 1;
    
    for ( size_t i=0; i<inputLength; i++ ) {

        // last character?
        BOOL lastCharacter = (i+one) >= inputLength;

        // check if we completed the word
        BOOL wordComplete = lastCharacter == TRUE || inputString[i+1] == ' ' || inputString[i+1] == '\n' || inputString[i+1] == '\r';
        
        BOOL startedNewLine = FALSE;
        
        // check if we've reached the end of the line or the end of the string
        if ( (i-lastCopied+one) >= (lineLength) || lastCharacter ) {
            
            // if the word doesn't fit on one line we'll have to break it up
            if ( lastCopied >= wordStart ) {
                wordComplete = TRUE;
            }
            
            // copy everything till the beginning of the last word
            const size_t copyUntil = (wordComplete==TRUE) ? i : wordStart-one;
            const size_t copyLength = copyUntil-lastCopied+one;
            memcpy(line, &inputString[lastCopied], sizeof(char)*copyLength);
            line[copyLength] = '\0';
            lastCopied = copyUntil+one;
            
            // allocate a new line
            size_t newSize = sizeof(char*) * ((size_t)((*nLines)+one));
            *lineArray = (char**)realloc(*lineArray, newSize);
            line = (char*)rsMalloc(lineSize);
            line[0] = '\0';
            (*lineArray)[*nLines] = line;
            (*nLines)++;
            startedNewLine = TRUE;
        }
        
        // If this word was completely copied, look for the beginning of the next word
        if ( wordComplete == TRUE ) {
            wordStart = i+1;
            for ( size_t j=wordStart; j<inputLength; j++ ) {
                if ( inputString[j] == ' ' || inputString[j] == '\n' || inputString[i+1] == '\r' ) {
                    continue;
                }
                wordStart = j;
                break;
            }
        }
        
        // if the last word was not yet completely copied, restart processing at its beginning
        if ( startedNewLine == TRUE ) {
            lastCopied = wordStart;
        }
    }
    
    // if the end of the string matches perfectly into the
    // last line we inadvertently created one line too many
    if ( lastCopied >= inputLength ) {
        (*nLines)--;
        rsFree((*lineArray)[*nLines]);
    }
}

/*
 * Concatenates an arbitrary number of strings(char*)
 * and returns the resulting string
 * 
 * Please note that the last element should always be
 * NULL. This is necessary for the method to know when
 * to stop looking for arguments.
 */
char* rsStringConcat(char *first, ...)
{
    // initialize variable arguments
    va_list ap;
    va_start(ap, first);

    size_t count = 0;
    char** strings = NULL;
    char* string = first;
    size_t length  = 1;
    
    // collect strings
    do {
        count++;
        strings = (char**)realloc(strings, sizeof(char*) * count); // extend string array by one
        length += strlen(string);
        strings[count-1] = string;
        string = va_arg(ap, char*);
    } while ( string != NULL );
    
    va_end(ap);
    
    // join them
    char* result = rsMalloc(sizeof(char) * length);
    size_t alreadyCopied = 0;
    
    for ( size_t i=0; i<count; i++ ) {
        sprintf(&result[alreadyCopied], "%s", strings[i]);
        alreadyCopied += strlen(strings[i]);
    }
    
    return result;
}

/*
 * Case-insensitive version of strcmp
 */
int rsStringCompareCaseInsensitive(char const *a, char const *b)
{
    for (;; a++, b++) {
        int d = tolower(*a) - tolower(*b);
        if (d != 0 || !*a)
            return d;
    }
}
