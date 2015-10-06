#ifndef rstools_niftiutils_headerinfo_h
#define rstools_niftiutils_headerinfo_h

#include "rscommon.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    short  RSToolsHeaderVersion;
    char   ImplementationVersionName[16+1];     // (0002,0013)
    char   AcquisitionDate[8+1];                // (0008,0022)
    char   InstitutionName[64+1];               // (0008,0080)
    char   InstitutionAddress[1024+1];          // (0008,0081)
    char   StudyDescription[64+1];              // (0008,1030)
    char   SeriesDescription[64+1];             // (0008,103e)
    char   OperatorsName[5*64+5];               // (0008,1070)
    char   ManufacturerModelName[64+1];         // (0008,1090)
    char   PatientName[5*64+5];                 // (0010,0010)
    char   PatientID[64+1];                     // (0010,0020)
    char   PatientBirthDate[8+1];               // (0010,0030)
    char   PatientSex[16+1];                    // (0010,0040)
    char   PatientAge[4+1];                     // (0010,1010)
    char   PatientWeight[16+1];                 // (0010,1030)
    char   SequenceName[16+1];                  // (0018,0024)
    char   SliceThickness[16+1];                // (0018,0050)
    char   RepetitionTime[16+1];                // (0018,0080)
    char   EchoTime[16+1];                      // (0018,0081)
    char   MagneticFieldStrength[16+1];         // (0018,0087)
    char   SpacingBetweenSlices[16+1];          // (0018,0088)
    char   NumberOfPhaseEncodingSteps[12+1];    // (0018,0089)
    char   PixelBandwidth[16+1];                // (0018,0095)
    char   SoftwareVersions[64+1];              // (0018,1020)
    char   ProtocolName[64+1];                  // (0018,1030)
    char   TransmitCoilName[16+1];              // (0018,1251)
    char   InPlanePhaseEncodingDirection[16+1]; // (0018,1312)
    char   PatientPosition[16+1];               // (0018,5100)
    char   SeriesNumber[12+1];                  // (0020,0011)
    char   ImageComments[1024+1];               // (0020,4000) ; potentially shortened (max size: 10240)
    char   MatrixSize[64+1];                    // (0051,100b)
    char   FieldOfView[64+1];                   // (0051,100c)
    unsigned short Rows;                        // (0028,0010)
    unsigned short Columns;                     // (0028,0011)
    char   GrappaFactor[10+1];                  // (0029,1020)/sPat.lAccelFactPE (siemens extended information)
    char   PhaseEncodingLines[10+1];            // (0029,1020)/sKSpace.lPhaseEncodingLines (siemens extended information)
    double DwellTime;                           // = 1/(BandwidthPerPixelPhaseEncode * PhaseEncodingLines)
    double BandwidthPerPixelPhaseEncode;        // (0019,1028)
    double MosaicRefAcqTimes[1024];             // (0019,1029)
} rsNiftiExtendedHeaderInformation;

typedef struct {
    const short tagGroup;
    const short tagElement;

    const char valueRepresentation[2];
} rsDicomElement;

rsNiftiExtendedHeaderInformation* rsNiftiInitializeExtendedHeaderInformation();
void rsNiftiPrintExtendedHeaderInformation(rsNiftiExtendedHeaderInformation* info);
void rsNiftiAddExtendedHeaderInformation(rsNiftiExtendedHeaderInformation* info, const rsDicomElement* dicomElement, char* buffer, size_t length);
rsNiftiExtendedHeaderInformation* rsNiftiFindExtendedHeaderInformation(nifti_image *nim);

#ifdef __cplusplus
}
#endif

#endif
