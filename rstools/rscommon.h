#define RSTOOLS_XSTR(x) #x
#define RSTOOLS_STR(x) RSTOOLS_XSTR(x)
#define RSTOOLS_VERSION_HASH_STRING RSTOOLS_STR(RSTOOLS_VERSION_HASH)
#define RSTOOLS_VERSION "0.7.1 (" RSTOOLS_VERSION_HASH_STRING ")"
#define RSTOOLS_VERSION_LABEL "RSTools " RSTOOLS_VERSION