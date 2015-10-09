#include <stddef.h>
#include <sys/stat.h>
#include <src/utils/rsstring.h>
#include "dicomfinder.h"
#include "rsniftiutils.h"
#include "headerinfo.h"

typedef struct {
    char file[4096+1];
    char patientId[64+1];
    char studyId[16+1];
    char seriesNumber[12+1];
    char instanceNumber[12+1];
    char hash[99+1];
    short group;
} BasicDICOMInfo;

typedef struct {
    char preamble[128];
    char prefix[4];
} BasicDICOMHeader;

BOOL rsDicomFinderIsDirectory(const char *path);
char **rsListDicomsInDirectory(const char *path, size_t *nFiles);
BOOL rsIsDicomFile(const char *path);
void rsReadBasicDicomInfo(const char *path, BasicDICOMInfo* info);
void rsCreateDicomHash(BasicDICOMInfo *info);
int rsShortComparator(const void *a, const void *b);

/**
 * This method searches for dicom series in a path and chooses the one with the
 * most volumes. Of this series a volume in the middle is returned which is assumed
 * to contain reliable information. Due to a bug, some Siemens header information
 * are somewhat off in the first dicom which is why we want to avoid extracting
 * information from this one.
 */
char *rsDicomFinderFindReliableDicomInPath(const char *path) {

    if (!rsDicomFinderIsDirectory(path)) {
        // input is not a directory, assume we were already given a dicom
        char *result = (char *)rsMalloc(sizeof(char) * (strlen(path) + 1));
        sprintf(result, "%s", path);
        return result;
    }

    // find dicoms
    size_t tmp;
    char **dicoms = rsListDicomsInDirectory(path, &tmp);
    const size_t nDicoms = tmp;

    // read in essential dicom header information
    BasicDICOMInfo **dicomInfos = (BasicDICOMInfo**)rsMalloc(sizeof(BasicDICOMInfo*)*nDicoms);
    for (size_t i=0; i<nDicoms; i++) {
        dicomInfos[i] = (BasicDICOMInfo*)rsMalloc(sizeof(BasicDICOMInfo));
        rsReadBasicDicomInfo(dicoms[i], dicomInfos[i]);
    }

    // form groups
    size_t hashlen = 99;
    size_t nScanGroups = 0;
    char **groupHashes = (char**)rsMalloc(sizeof(char*)*nDicoms);
    const size_t max_size_t = (size_t)-1;

    for (size_t i=0; i<nDicoms; i++) {
        BasicDICOMInfo *info = dicomInfos[i];
        info->group = max_size_t; // initialize with max size_t

        // check if it can be assigned to existing group
        for (short g=0; g<nScanGroups; g++) {
            const char* groupHash = groupHashes[g];
            const int comp = strcmp(info->hash, groupHash);
            if (comp == 0) {
                info->group = g;
                break;
            }
        }

        // create new group if none was found
        if (dicomInfos[i]->group == max_size_t) {
            groupHashes[nScanGroups] =  (char*)rsMalloc(sizeof(char)*(hashlen+1));
            sprintf(groupHashes[nScanGroups], "%s", info->hash);
            nScanGroups++;
        }
    }

    // count number of occurrences of every group
    size_t *nGroupCount = (size_t *)rsMalloc(nScanGroups*sizeof(size_t));
    for (short g=0; g<nScanGroups; g++) {
        nGroupCount[g] = 0; // initialize with zero
    }

    for (size_t i = 0; i < nDicoms; i++) {
        const BasicDICOMInfo *info = dicomInfos[i];

        for (short g=0; g<nScanGroups; g++) {
            if (strcmp(info->hash, groupHashes[g])== 0) {
                nGroupCount[g]++;
                break;
            }
        }
    }

    // choose group with highest number of dicoms
    short highestIndex = 0;
    size_t highestCount = 0;
    for (short g=0; g<nScanGroups; g++) {
        if (highestCount <= nGroupCount[g]) {
            highestCount = nGroupCount[g];
            highestIndex = g;
        }
    }

    // collect instance numbers for the biggest group
    short j = 0;
    short *instanceNumbers = (short*)rsMalloc(highestCount*sizeof(short));
    for (size_t i = 0; i < nDicoms; i++) {
        const BasicDICOMInfo *info = dicomInfos[i];
        if (strcmp(info->hash, groupHashes[highestIndex]) == 0) {
            instanceNumbers[j] = (short)atoi(&info->instanceNumber[0]);
            j++;
        }
    }

    // sort instance numbers
    qsort(instanceNumbers, highestCount, sizeof(short), rsShortComparator);

    // choose one instance number in the middle
    const short referenceInstance = (short)ceil((highestCount-1)/2);

    // find that instance
    BasicDICOMInfo *referenceInfo = NULL;
    for (size_t i = 0; i < nDicoms; i++) {
        BasicDICOMInfo *info = dicomInfos[i];
        short instanceNumber = (short)atoi(&info->instanceNumber[0]);

        if (instanceNumber != referenceInstance) {
            continue;
        }

        if (strcmp(info->hash, groupHashes[highestIndex]) == 0) {
            referenceInfo = info;
            break;
        }
    }

    // save resulting path
    char *result = NULL;
    if (referenceInfo != NULL) {
        result = (char *)rsMalloc(sizeof(char) * (strlen(referenceInfo->file) + 1));
        sprintf(result, "%s", referenceInfo->file);
    }

    // cleanup
    rsFree(instanceNumbers);
    for (size_t i = 0; i < nDicoms; i++) {
        rsFree(dicomInfos[i]);
    }
    rsFree(dicomInfos);

    for (short g=0; g<nScanGroups; g++) {
        rsFree(groupHashes[g]);
    }
    rsFree(groupHashes);
    rsFree(nGroupCount);

    return result;
}

BOOL rsDicomFinderIsDirectory(const char *path)
{
    struct stat s;
    return stat(path, &s) == 0 && (s.st_mode & S_IFDIR);
}

char **rsListDicomsInDirectory(const char *path, size_t *nFiles)
{
    *nFiles = 0;
    DIR *dir;
    struct dirent *ent;

    // remove trailing slash if present
    const size_t pathlen = strlen(path);
    char dirpath[pathlen+1];
    sprintf(&dirpath[0], "%s", path);
    if (dirpath[pathlen-1] == '/') {
        dirpath[pathlen-1] = '\0';
    }

    // count number of dicoms in the directory
    if ((dir = opendir(dirpath)) == NULL) return NULL;

    while ((ent = readdir(dir)) != NULL) {
        // skip '.', '..' and any hidden files
        if (ent->d_name[0] == '.') continue;

        // test if dicom
        char *fullpath = rsStringConcat(dirpath, "/", ent->d_name);
        if (rsIsDicomFile(fullpath)) {
            *nFiles = *nFiles + 1L;
        }
        rsFree(fullpath);
    }
    closedir(dir);

    // save dicom paths in array
    char **files = (char**)rsMalloc(sizeof(char**) * *nFiles);

    if ((dir = opendir(path)) == NULL) return NULL;

    size_t i=0;
    while ((ent = readdir(dir)) != NULL) {
        // skip '.', '..' and any hidden files
        if (ent->d_name[0] == '.') continue;

        // test if dicom
        char *fullpath = rsStringConcat(dirpath, "/", ent->d_name);
        if (rsIsDicomFile(fullpath)) {
            files[i] = (char*)rsMalloc(sizeof(char)*(strlen(fullpath)+1L));
            sprintf(files[i], "%s", fullpath);
            i++;
        }
        rsFree(fullpath);
    }
    closedir(dir);

    return &files[0];
}

BOOL rsIsDicomFile(const char *path)
{
    // open dicom
    FILE *file = fopen(path, "r");
    if (!file) return FALSE;

    // read in header
    BasicDICOMHeader header;
    BOOL read = fread(&header, sizeof(BasicDICOMHeader), 1, file) > 0;
    fclose(file);

    if (!read) return FALSE;

    // check if DICOM
    return header.prefix[0] == 'D'
        && header.prefix[1] == 'I'
        && header.prefix[2] == 'C'
        && header.prefix[3] == 'M';
}

void rsReadBasicDicomInfo(const char *path, BasicDICOMInfo* info)
{
    sprintf(&info->file[0], "%s", path);

    // open dicom
    FILE *f;
    f = fopen(path, "r");

    // determine size
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    // read into buffer
    char* buffer = (char*)malloc(size * sizeof(char));
    fread(buffer, sizeof(char), size, f);

    // close dicom
    fclose(f);

    size_t lengthPreamble = 128 * sizeof(char);
    size_t lengthPrefix   =   4 * sizeof(char);
    size_t offsetData     = lengthPreamble + lengthPrefix;

    size_t i=offsetData;
    size_t headerLength = 0;

    // walk through dicom tags
    do {
        const rsDicomElement* currentElement = (rsDicomElement*)&buffer[i];

        // determine the length of the value length field
        size_t valueLength = rsNiftiGetDicomValueLength(currentElement->valueRepresentation);

        THIS_UINT32 length = 0;

        if ( valueLength == 2 ) {
            length = *((THIS_UINT16*)(&(currentElement->valueRepresentation[0]) + 2));
        } else {
            length = *((THIS_UINT32*)(&(currentElement->valueRepresentation[0]) + 4) );
        }

        // calculate next tag offset
        const size_t lengthTag   = 4;
        const size_t lengthVR    = valueLength==2 ? 2 : 4;
        const size_t lengthVL    = valueLength;
        const size_t lengthV     = length;
        const size_t offsetValue = lengthTag + lengthVR + lengthVL;
        const size_t offset      = offsetValue + lengthV;

        // check if we reached the data field
        if ( currentElement->tagGroup == 0x7fe0 && currentElement->tagElement == 0x10 ) {
            break;
        }

        char* valueBuffer = &buffer[i+offsetValue];

        // save field if required
        switch (currentElement->tagGroup) {
            case 0x0010:
                if (currentElement->tagElement == 0x0020 && lengthV < sizeof(info->patientId)) {
                    char *buf = (char*)rsMalloc(sizeof(info->patientId)*sizeof(char));
                    buf[0] = '\0';
                    strncat(buf, valueBuffer, lengthV);
                    sprintf(info->patientId, "%s", rsTrimString(buf));
                    rsFree(buf);
                }
                break;
            case 0x0020:
                if (currentElement->tagElement == 0x0010 && lengthV < sizeof(info->studyId)) {
                    char *buf = (char*)rsMalloc(sizeof(info->studyId)*sizeof(char));
                    buf[0] = '\0';
                    strncat(buf, valueBuffer, lengthV);
                    sprintf(&info->studyId[0], "%s", rsTrimString(buf));
                    rsFree(buf);
                } else if  (currentElement->tagElement == 0x0011 && lengthV < sizeof(info->seriesNumber)) {
                    char *buf = (char*)rsMalloc(sizeof(info->seriesNumber)*sizeof(char));
                    buf[0] = '\0';
                    strncat(buf, valueBuffer, lengthV);
                    sprintf(&info->seriesNumber[0], "%s", rsTrimString(buf));
                    rsFree(buf);
                } else if  (currentElement->tagElement == 0x0013 && lengthV < sizeof(info->instanceNumber)) {
                    char *buf = (char*)rsMalloc(sizeof(info->instanceNumber)*sizeof(char));
                    buf[0] = '\0';
                    strncat(buf, valueBuffer, lengthV);
                    sprintf(&info->instanceNumber[0], "%s", rsTrimString(buf));
                    rsFree(buf);
                }
                break;
        }

        i+= offset;

    } while (i<size);

    rsCreateDicomHash(info);

    rsFree(buffer);
}

void rsCreateDicomHash(BasicDICOMInfo *info)
{
    sprintf(&info->hash[0], "%s__%s__%s", &info->patientId[0], &info->studyId[0], &info->seriesNumber[0]);
}

int rsShortComparator(const void *a, const void *b)
{
    const short *da = (const short *) a;
    const short *db = (const short *) b;

    return (*da > *db) - (*da < *db);
}