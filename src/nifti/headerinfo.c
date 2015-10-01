#include <src/utils/rsstring.h>
#include "headerinfo.h"

void rsNiftiReadTagFromSiemensExtraInformationHeader(char* value, const char* buffer, const size_t bufferLength, const char* entryTagName, const size_t maxExpectedLength);

rsNiftiExtendedHeaderInformation* rsNiftiInitializeExtendedHeaderInformation()
{
    rsNiftiExtendedHeaderInformation* info = (rsNiftiExtendedHeaderInformation*)rsMalloc(sizeof(rsNiftiExtendedHeaderInformation));
    info->RSToolsHeaderVersion = 1;

    // initialize all other fields with a default value
    info->ImplementationVersionName[0]     = '\0';
    info->AcquisitionDate[0]               = '\0';
    info->InstitutionName[0]               = '\0';
    info->InstitutionAddress[0]            = '\0';
    info->StudyDescription[0]              = '\0';
    info->SeriesDescription[0]             = '\0';
    info->OperatorsName[0]                 = '\0';
    info->ManufacturerModelName[0]         = '\0';
    info->PatientName[0]                   = '\0';
    info->PatientID[0]                     = '\0';
    info->PatientBirthDate[0]              = '\0';
    info->PatientSex[0]                    = '\0';
    info->PatientAge[0]                    = '\0';
    info->PatientWeight[0]                 = '\0';
    info->SequenceName[0]                  = '\0';
    info->SliceThickness[0]                = '\0';
    info->RepetitionTime[0]                = '\0';
    info->EchoTime[0]                      = '\0';
    info->MagneticFieldStrength[0]         = '\0';
    info->SpacingBetweenSlices[0]          = '\0';
    info->NumberOfPhaseEncodingSteps[0]    = '\0';
    info->PixelBandwidth[0]                = '\0';
    info->SoftwareVersions[0]              = '\0';
    info->ProtocolName[0]                  = '\0';
    info->TransmitCoilName[0]              = '\0';
    info->InPlanePhaseEncodingDirection[0] = '\0';
    info->PatientPosition[0]               = '\0';
    info->SeriesNumber[0]                  = '\0';
    info->ImageComments[0]                 = '\0';
    info->MatrixSize[0]                    = '\0';
    info->FieldOfView[0]                   = '\0';
    info->GrappaFactor[0]                  = '\0';
    info->BandwidthPerPixelPhaseEncode     = log(-1);
    info->DwellTime                        = log(-1);
    info->Rows                             = log(-1);
    info->Columns                          = log(-1);

    // initialize slice acquisition times with NaN
    for (short i=0; i<sizeof(info->MosaicRefAcqTimes)/sizeof(double); i++) {
        info->MosaicRefAcqTimes[i] = log(-1.0);
    }

    return info;
}

void rsNiftiPrintExtendedHeaderInformation(rsNiftiExtendedHeaderInformation* info) {
    fprintf(stdout, "ImplementationVersionName    : %s\n", info->ImplementationVersionName);
    fprintf(stdout, "AcquisitionDate              : %s\n", info->AcquisitionDate);
    fprintf(stdout, "InstitutionName              : %s\n", info->InstitutionName);
    fprintf(stdout, "InstitutionAddress           : %s\n", info->InstitutionAddress);
    fprintf(stdout, "StudyDescription             : %s\n", info->StudyDescription);
    fprintf(stdout, "SeriesDescription            : %s\n", info->SeriesDescription);
    fprintf(stdout, "OperatorsName                : %s\n", info->OperatorsName);
    fprintf(stdout, "ManufacturerModelName        : %s\n", info->ManufacturerModelName);
    fprintf(stdout, "PatientName                  : %s\n", info->PatientName);
    fprintf(stdout, "PatientID                    : %s\n", info->PatientID);
    fprintf(stdout, "PatientBirthDate             : %s\n", info->PatientBirthDate);
    fprintf(stdout, "PatientSex                   : %s\n", info->PatientSex);
    fprintf(stdout, "PatientAge                   : %s\n", info->PatientAge);
    fprintf(stdout, "PatientWeight                : %s\n", info->PatientWeight);
    fprintf(stdout, "SequenceName                 : %s\n", info->SequenceName);
    fprintf(stdout, "SliceThickness               : %s\n", info->SliceThickness);
    fprintf(stdout, "RepetitionTime               : %s\n", info->RepetitionTime);
    fprintf(stdout, "EchoTime                     : %s\n", info->EchoTime);
    fprintf(stdout, "MagneticFieldStrength        : %s\n", info->MagneticFieldStrength);
    fprintf(stdout, "SpacingBetweenSlices         : %s\n", info->SpacingBetweenSlices);
    fprintf(stdout, "NumberOfPhaseEncodingSteps   : %s\n", info->NumberOfPhaseEncodingSteps);
    fprintf(stdout, "PixelBandwidth               : %s\n", info->PixelBandwidth);
    fprintf(stdout, "SoftwareVersions             : %s\n", info->SoftwareVersions);
    fprintf(stdout, "ProtocolName                 : %s\n", info->ProtocolName);
    fprintf(stdout, "TransmitCoilName             : %s\n", info->TransmitCoilName);
    fprintf(stdout, "InPlanePhaseEncodingDirection: %s\n", info->InPlanePhaseEncodingDirection);
    fprintf(stdout, "PatientPosition              : %s\n", info->PatientPosition);
    fprintf(stdout, "BandwidthPerPixelPhaseEncode : %.15f\n", info->BandwidthPerPixelPhaseEncode);
    fprintf(stdout, "SeriesNumber                 : %s\n", info->SeriesNumber);
    fprintf(stdout, "ImageComments                : %s\n", info->ImageComments);
    fprintf(stdout, "MatrixSize                   : %s\n", info->MatrixSize);
    fprintf(stdout, "FieldOfView                  : %s\n", info->FieldOfView);
    fprintf(stdout, "GrappaFactor                 : %s\n", info->GrappaFactor);
    fprintf(stdout, "PhaseEncodingLines           : %s\n", info->PhaseEncodingLines);
    fprintf(stdout, "DwellTime                    : %.15f\n", info->DwellTime);
    fprintf(stdout, "Rows                         : %d\n", (int)info->Rows);
    fprintf(stdout, "Columns                      : %d\n", (int)info->Columns);

    short nSlicesTimes = 0;
    for (nSlicesTimes=0; (nSlicesTimes<sizeof(info->MosaicRefAcqTimes)/sizeof(double) && !isnan(info->MosaicRefAcqTimes[nSlicesTimes])); nSlicesTimes++);;
    fprintf(stdout, "MosaicRefAcqTimes            : [%d values]\n", (int)nSlicesTimes);
    for (short i=0; i<nSlicesTimes; i++) {
        fprintf(stdout, "    %.1f\n", info->MosaicRefAcqTimes[i]);
    }
}

void rsNiftiAddExtendedHeaderInformation(rsNiftiExtendedHeaderInformation* info, const rsDicomElement* dicomElement, char* buffer, size_t length) {
    switch (dicomElement->tagGroup) {
        case 0x0002:
            if (dicomElement->tagElement == 0x0013 && length < sizeof(info->ImplementationVersionName)) {
                strncat(&info->ImplementationVersionName[0], buffer, length);
            }
            break;
        case 0x0008:
            switch (dicomElement->tagElement) {
                case 0x0022:
                    if (length < sizeof(info->AcquisitionDate))
                        strncat(&info->AcquisitionDate[0], buffer, length);
                    break;
                case 0x0080:
                    if (length < sizeof(info->InstitutionName))
                        strncat(&info->InstitutionName[0], buffer, length);
                    break;
                case 0x0081:
                    if (length < sizeof(info->InstitutionAddress))
                        strncat(&info->InstitutionAddress[0], buffer, length);
                    break;
                case 0x1030:
                    if (length < sizeof(info->StudyDescription))
                        strncat(&info->StudyDescription[0], buffer, length);
                    break;
                case 0x103e:
                    if (length < sizeof(info->SeriesDescription))
                        strncat(&info->SeriesDescription[0], buffer, length);
                    break;
                case 0x1070:
                    if (length < sizeof(info->OperatorsName))
                        strncat(&info->OperatorsName[0], buffer, length);
                    break;
                case 0x1090:
                    if (length < sizeof(info->ManufacturerModelName))
                        strncat(&info->ManufacturerModelName[0], buffer, length);
                    break;
            }
            break;
        case 0x0010:
            switch (dicomElement->tagElement) {
                case 0x0010:
                    if (length < sizeof(info->PatientName))
                        strncat(&info->PatientName[0], buffer, length);
                    break;
                case 0x0020:
                    if (length < sizeof(info->PatientID))
                        strncat(&info->PatientID[0], buffer, length);
                    break;
                case 0x0030:
                    if (length < sizeof(info->PatientBirthDate))
                        strncat(&info->PatientBirthDate[0], buffer, length);
                    break;
                case 0x0040:
                    if (length < sizeof(info->PatientSex))
                        strncat(&info->PatientSex[0], buffer, length);
                    break;
                case 0x1010:
                    if (length < sizeof(info->PatientAge))
                        strncat(&info->PatientAge[0], buffer, length);
                    break;
                case 0x1030:
                    if (length < sizeof(info->PatientWeight))
                        strncat(&info->PatientWeight[0], buffer, length);
                    break;
            }
            break;
        case 0x0018:
            switch (dicomElement->tagElement) {
                case 0x0024:
                    if (length < sizeof(info->SequenceName))
                        strncat(&info->SequenceName[0], buffer, length);
                    break;
                case 0x0050:
                    if (length < sizeof(info->SliceThickness))
                        strncat(&info->SliceThickness[0], buffer, length);
                    break;
                case 0x0080:
                    if (length < sizeof(info->RepetitionTime))
                        strncat(&info->RepetitionTime[0], buffer, length);
                    break;
                case 0x0081:
                    if (length < sizeof(info->EchoTime))
                        strncat(&info->EchoTime[0], buffer, length);
                    break;
                case 0x0087:
                    if (length < sizeof(info->MagneticFieldStrength))
                        strncat(&info->MagneticFieldStrength[0], buffer, length);
                    break;
                case 0x0088:
                    if (length < sizeof(info->SpacingBetweenSlices))
                        strncat(&info->SpacingBetweenSlices[0], buffer, length);
                    break;
                case 0x0089:
                    if (length < sizeof(info->NumberOfPhaseEncodingSteps))
                        strncat(&info->NumberOfPhaseEncodingSteps[0], buffer, length);
                    break;
                case 0x0095:
                    if (length < sizeof(info->PixelBandwidth))
                        strncat(&info->PixelBandwidth[0], buffer, length);
                    break;
                case 0x1020:
                    if (length < sizeof(info->SoftwareVersions))
                        strncat(&info->SoftwareVersions[0], buffer, length);
                    break;
                case 0x1030:
                    if (length < sizeof(info->ProtocolName))
                        strncat(&info->ProtocolName[0], buffer, length);
                    break;
                case 0x1251:
                    if (length < sizeof(info->TransmitCoilName))
                        strncat(&info->TransmitCoilName[0], buffer, length);
                    break;
                case 0x1312:
                    if (length < sizeof(info->InPlanePhaseEncodingDirection))
                        strncat(&info->InPlanePhaseEncodingDirection[0], buffer, length);
                    break;
                case 0x5100:
                    if (length < sizeof(info->PatientPosition))
                        strncat(&info->PatientPosition[0], buffer, length);
                    break;
            }
            break;
        case 0x0019:
            if (dicomElement->tagElement == 0x1028) {
                const double* bandwidthPerPixelPhaseEncode = (double*)buffer;
                info->BandwidthPerPixelPhaseEncode = *bandwidthPerPixelPhaseEncode;
            } else if (dicomElement->tagElement == 0x1029) {
                const short maxNumberOfSlices = sizeof(info->MosaicRefAcqTimes) / sizeof(double);
                const short presentNumberOfSlices = length / sizeof(double);
                const short applicableSlices = fmin(maxNumberOfSlices, presentNumberOfSlices);
                const double* sliceTimings = (double*)buffer;
                for (short i=0; i<applicableSlices; i++) {
                    info->MosaicRefAcqTimes[i] = sliceTimings[i];
                }
            }
            break;
        case 0x0020:
            if (dicomElement->tagElement == 0x0011 && length < sizeof(info->SeriesNumber)) {
                strncat(&info->SeriesNumber[0], buffer, length);
            } else if (dicomElement->tagElement == 0x4000) {
                strncat(&info->ImageComments[0], buffer, fmin(length, sizeof(info->ImageComments)-1)); // truncate
            }
            break;
        case 0x0028:
            if (dicomElement->tagElement == 0x0010) {
                info->Rows = *((unsigned short*)buffer);
            } else if (dicomElement->tagElement == 0x0011) {
                info->Columns = *((unsigned short*)buffer);
            }
            break;
        case 0x0029:
            if (dicomElement->tagElement == 0x1020) {
                rsNiftiReadTagFromSiemensExtraInformationHeader(&info->GrappaFactor[0], &buffer[0], length, "sPat.lAccelFactPE", sizeof(info->GrappaFactor)-1);
                rsNiftiReadTagFromSiemensExtraInformationHeader(&info->PhaseEncodingLines[0], &buffer[0], length, "sKSpace.lPhaseEncodingLines", sizeof(info->GrappaFactor)-1);
            }
            break;
        case 0x0051:
            if (dicomElement->tagElement == 0x100b && length < sizeof(info->MatrixSize)) {
                strncat(&info->MatrixSize[0], buffer, length);
            } else if (dicomElement->tagElement == 0x100c && length < sizeof(info->FieldOfView)) {
                strncat(&info->FieldOfView[0], buffer, length);
            }
            break;
    }
}

void rsNiftiReadTagFromSiemensExtraInformationHeader(char* value, const char* buffer, const size_t bufferLength, const char* entryTagName, const size_t maxExpectedLength)
{
    // find the siemens extra info awesomeness
    const char *startToken = "### ASCCONV BEGIN";
    const char *endToken = "### ASCCONV END";
    const char *seiStart = memmem(&buffer[0], bufferLength, startToken, strlen(startToken));
    if (!seiStart) return;

    const char *seiEnd = memmem(&seiStart[0], bufferLength-(buffer-seiStart), endToken, strlen(endToken));
    if (!seiEnd) return;

    // if we got here we've found a proper siemens extra info definition, so let's search for
    // the info we're interested in, e.g.
    // entry.tag.name                        = value

    char *entryStart = memmem(&seiStart[0], bufferLength-(buffer-seiStart), entryTagName, strlen(entryTagName));
    if (!entryStart) return;

    // find '='
    entryStart = memmem(&entryStart[0], bufferLength-(buffer- entryStart), "=", strlen("="));
    if (!entryStart) return;

    // forward to the number/value
    entryStart += 2;

    // find the end of it
    const char *entryEnd = memmem(&entryStart[0], bufferLength-(buffer- entryStart), "\n", strlen("\n"));
    if (!entryEnd) return;

    // copy the result
    const size_t entryLength = entryEnd - entryStart;
    const size_t resultLength = fmin(entryLength, maxExpectedLength);
    strncat(&value[0], entryStart, resultLength);
}