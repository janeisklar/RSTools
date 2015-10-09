#ifndef RSTOOLS_DICOMFINDER_H
#define RSTOOLS_DICOMFINDER_H

/**
 * This method searches for dicom series in a path and chooses the one with the
 * most volumes. Of this series a volume in the middle is returned which is assumed
 * to contain reliable information. Due to a bug, some Siemens header information
 * are somewhat off in the first dicom which is why we want to avoid extracting
 * information from this one.
 */
char * rsDicomFinderFindReliableDicomInPath(const char *path);

#endif //RSTOOLS_DICOMFINDER_H
