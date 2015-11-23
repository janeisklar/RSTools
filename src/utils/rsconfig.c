#include <stdio.h>
#include "rscommon.h"
#include "rsmemory.h"
#include "rsstring.h"
#include "rsconfig.h"

char *rsConfigGetConfigSetting(const char* preamble, const char* description)
{
    const char *confPath = CONFIG_PATH"/rstools.conf";
    FILE *config = fopen(confPath, "r");

    if (config == NULL) {
        fprintf(stderr, "Could not read rstools config file to retrieve %s (%s)\n", description, confPath);
        return NULL;
    }

    int length;
    char *line = (char*)rsMalloc(sizeof(char)*1000);
    char *value = NULL;

    while (rsReadline(config, line, &length)) {
        if (!rsStringStartsWith(line, preamble)) {
            continue;
        }
        value = rsString(&line[strlen(preamble)]);
    }

    fclose(config);
    rsFree(line);

    if (value == NULL) {
        fprintf(stderr, "The rstools config file did not contain %s (%s)\n", description, confPath);
        return NULL;
    }

    return value;
}

char *rsConfigGetANTsPath()
{
    return rsConfigGetConfigSetting("ANTSPATH=", "path to ANTs");
}

const char *rsConfigGetRSToolsExecutablesPath()
{
    return RSTOOLS_EXECUTABLES_PATH;
}
