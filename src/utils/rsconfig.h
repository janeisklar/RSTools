#ifndef RSTOOLS_RSCONFIG_H
#define RSTOOLS_RSCONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

char *rsConfigGetConfigSetting(const char* preamble, const char* description);
char *rsConfigGetANTsPath();
const char *rsConfigGetRSToolsExecutablesPath();

#ifdef __cplusplus
}
#endif

#endif //RSTOOLS_RSCONFIG_H
