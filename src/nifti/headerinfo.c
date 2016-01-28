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
    info->PhaseEncodingDirection[0]        = '\0';
    info->BandwidthPerPixelPhaseEncode     = log(-1);
    info->DwellTime                        = log(-1);
    info->Rows                             = 0;
    info->Columns                          = 0;

    // initialize slice acquisition times with NaN
    for (short i=0; i<sizeof(info->MosaicRefAcqTimes)/sizeof(double); i++) {
        info->MosaicRefAcqTimes[i] = log(-1.0);
    }

    return info;
}

void* rsNiftiExtendendHeaderInformationGet(void* self) {
    return ((rsNiftiExtendedHeaderInformationEntry *)self)->data;
}

void rsNiftiExtendendHeaderInformationSetString(void* self, const void* value) {
    sprintf((char*)((rsNiftiExtendedHeaderInformationEntry *)self)->data, "%s", (char*)value);
}

void rsNiftiExtendendHeaderInformationSetDouble(void* self, const void* value) {
    const double v = *((double*)value);
    *((double*)((rsNiftiExtendedHeaderInformationEntry *)self)->data) = v;
}

void rsNiftiExtendendHeaderInformationParseDouble(void* self, const void* value) {
    const char* vString = (char*)value;
    *((double*)((rsNiftiExtendedHeaderInformationEntry *)self)->data) = strtod(vString, NULL);
}

void rsNiftiExtendendHeaderInformationSetUnsignedShort(void* self, const void* value) {
    const unsigned short v = *((unsigned short*)value);
    *((unsigned short*)((rsNiftiExtendedHeaderInformationEntry *)self)->data) = v;
}

void rsNiftiExtendendHeaderInformationSetShort(void* self, const void* value) {
    const short v = *((short*)value);
    *((short*)((rsNiftiExtendedHeaderInformationEntry *)self)->data) = v;
}

void rsNiftiExtendendHeaderInformationParseUnsignedShort(void* self, const void* value) {
    const char* vString = (char*)value;
    *((unsigned short*)((rsNiftiExtendedHeaderInformationEntry *)self)->data) = (unsigned short)atoi(vString);
}

void rsNiftiExtendendHeaderInformationParseShort(void* self, const void* value) {
    const char* vString = (char*)value;
    *((short*)((rsNiftiExtendedHeaderInformationEntry *)self)->data) = (short)atoi(vString);
}

void rsNiftiExtendendHeaderInformationSetSliceTimings(void* self, const void* value) {
    double* values = (double*)(((rsNiftiExtendedHeaderInformationEntry *)self)->data);
    const double* newValues = *((double**)value);

    short nSlicesTimes = 0;
    for (nSlicesTimes=0; (nSlicesTimes<1024 && !isnan(newValues[nSlicesTimes])); nSlicesTimes++) {
        values[nSlicesTimes] = newValues[nSlicesTimes];
    }

    if (nSlicesTimes<1024) {
        values[nSlicesTimes] = log(-1.0); // terminate with nan
    }
}

void rsNiftiExtendendHeaderInformationParseSliceTimings(void* self, const void* value) {
    char* vString = rsString((char*)value);
    double* values = (double*)(((rsNiftiExtendedHeaderInformationEntry *)self)->data);

    char *tmp = vString;
    char *current;
    size_t i = 0;

    // reset current values
    for (i=0; i<1024; i++) {
        values[i] = log(-1.0);
    }

    // parse and set given values
    i = 0;
    do {
        current = strtok(tmp, ",");
        tmp = NULL;

        if (current == NULL) {
            break;
        }

        values[i++] = strtod(current, NULL);
    } while (TRUE);
}

char* rsNiftiExtendendHeaderInformationFormatString(void* self, BOOL shorten) {
    return rsString((char*)(((rsNiftiExtendedHeaderInformationEntry *)self)->data));
}

void* rsNiftiExtendendHeaderInformationFormatDouble(void* self, BOOL shorten) {
    const double* value = (double*)(((rsNiftiExtendedHeaderInformationEntry *)self)->data);
    char *string = (char*)rsMalloc(sizeof(char)*17L);
    sprintf(string, "%.15f", *value);
    return string;
}

void* rsNiftiExtendendHeaderInformationFormatUnsignedShort(void* self, BOOL shorten) {
    const unsigned short value = *((unsigned short*)(((rsNiftiExtendedHeaderInformationEntry *)self)->data));
    char *string = (char*)rsMalloc(sizeof(char)*20L);
    sprintf(string, "%d", value);
    return string;
}

void* rsNiftiExtendendHeaderInformationFormatShort(void* self, BOOL shorten) {
    const short value = *((short*)(((rsNiftiExtendedHeaderInformationEntry *)self)->data));
    char *string = (char*)rsMalloc(sizeof(char)*20L);
    sprintf(string, "%d", value);
    return string;
}

void* rsNiftiExtendendHeaderInformationFormatSliceTimes(void* self, BOOL shorten) {

    const rsNiftiExtendedHeaderInformationEntry* entry  = (rsNiftiExtendedHeaderInformationEntry *)self;
    const double* values = (double*)(entry->data);

    short nSlicesTimes = 0;
    for (nSlicesTimes=0; (nSlicesTimes<1024 && !isnan(values[nSlicesTimes])); nSlicesTimes++)
        ;;

    if (shorten) {
        char *string = (char*)rsMalloc(sizeof(char)*14L);
        sprintf(string, "[%d values]", (int)nSlicesTimes);
        return string;
    }

    char *string = (char*)rsMalloc(sizeof(char)*9L*((size_t)nSlicesTimes));
    char *stringEnd = string;
    for (int i=0; i<nSlicesTimes; i++) {
        stringEnd += sprintf(stringEnd, "%.1f\n", values[i]);
    }
    return string;
}

rsNiftiExtendedHeaderInformationEntry * rsNiftiExtendendHeaderInformationListCreateEntry(char* key, rsNiftiExntededHeaderInformationGetterFunc get, rsNiftiExntededHeaderInformationSetterFunc set, rsNiftiExntededHeaderInformationSetterFunc parse, rsNiftiExntededHeaderInformationFormatterFunc format, void *data) {
    rsNiftiExtendedHeaderInformationEntry * entry = (rsNiftiExtendedHeaderInformationEntry *)rsMalloc(sizeof(rsNiftiExtendedHeaderInformationEntry));
    entry->key = key;
    entry->format = format;
    entry->data = data;
    entry->get = get;
    entry->set = set;
    entry->parse = parse;
    return entry;
}

rsNiftiExtendedHeaderInformationEntry ** rsNiftiExtendendHeaderInformationListCreateEntryMap(rsNiftiExtendedHeaderInformation* info) {
    rsNiftiExntededHeaderInformationGetterFunc get = (rsNiftiExntededHeaderInformationGetterFunc)rsNiftiExtendendHeaderInformationGet;
    rsNiftiExntededHeaderInformationSetterFunc setString = (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationSetString;
    rsNiftiExntededHeaderInformationSetterFunc setDouble = (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationSetDouble;
    rsNiftiExntededHeaderInformationSetterFunc setUShort = (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationSetShort;
    rsNiftiExntededHeaderInformationSetterFunc setShort = (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationSetUnsignedShort;
    rsNiftiExntededHeaderInformationSetterFunc setSliceTimings= (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationSetSliceTimings;
    rsNiftiExntededHeaderInformationSetterFunc parseDouble = (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationParseDouble;
    rsNiftiExntededHeaderInformationSetterFunc parseUShort = (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationParseUnsignedShort;
    rsNiftiExntededHeaderInformationSetterFunc parseShort = (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationParseShort;
    rsNiftiExntededHeaderInformationSetterFunc parseSliceTimings = (rsNiftiExntededHeaderInformationSetterFunc)rsNiftiExtendendHeaderInformationParseSliceTimings;
    rsNiftiExntededHeaderInformationFormatterFunc formatString = (rsNiftiExntededHeaderInformationFormatterFunc)rsNiftiExtendendHeaderInformationFormatString;
    rsNiftiExntededHeaderInformationFormatterFunc formatDouble = (rsNiftiExntededHeaderInformationFormatterFunc)rsNiftiExtendendHeaderInformationFormatDouble;
    rsNiftiExntededHeaderInformationFormatterFunc formatUShort= (rsNiftiExntededHeaderInformationFormatterFunc)rsNiftiExtendendHeaderInformationFormatUnsignedShort;
    rsNiftiExntededHeaderInformationFormatterFunc formatShort= (rsNiftiExntededHeaderInformationFormatterFunc)rsNiftiExtendendHeaderInformationFormatShort;
    rsNiftiExntededHeaderInformationFormatterFunc formatSliceTimings= (rsNiftiExntededHeaderInformationFormatterFunc)rsNiftiExtendendHeaderInformationFormatSliceTimes;

    rsNiftiExtendedHeaderInformationEntry * entries[] = {
        rsNiftiExtendendHeaderInformationListCreateEntry("RSToolsHeaderVersion",          get, setShort,  parseShort,  formatShort,  &info->RSToolsHeaderVersion),
        rsNiftiExtendendHeaderInformationListCreateEntry("ImplementationVersionName",     get, setString, setString,   formatString, &info->ImplementationVersionName),
        rsNiftiExtendendHeaderInformationListCreateEntry("AcquisitionDate",               get, setString, setString,   formatString, &info->AcquisitionDate),
        rsNiftiExtendendHeaderInformationListCreateEntry("InstitutionName",               get, setString, setString,   formatString, &info->InstitutionName),
        rsNiftiExtendendHeaderInformationListCreateEntry("InstitutionAddress",            get, setString, setString,   formatString, &info->InstitutionAddress),
        rsNiftiExtendendHeaderInformationListCreateEntry("StudyDescription",              get, setString, setString,   formatString, &info->StudyDescription),
        rsNiftiExtendendHeaderInformationListCreateEntry("SeriesDescription",             get, setString, setString,   formatString, &info->SeriesDescription),
        rsNiftiExtendendHeaderInformationListCreateEntry("OperatorsName",                 get, setString, setString,   formatString, &info->OperatorsName),
        rsNiftiExtendendHeaderInformationListCreateEntry("ManufacturerModelName",         get, setString, setString,   formatString, &info->ManufacturerModelName),
        rsNiftiExtendendHeaderInformationListCreateEntry("PatientName",                   get, setString, setString,   formatString, &info->PatientName),
        rsNiftiExtendendHeaderInformationListCreateEntry("PatientID",                     get, setString, setString,   formatString, &info->PatientID),
        rsNiftiExtendendHeaderInformationListCreateEntry("PatientBirthDate",              get, setString, setString,   formatString, &info->PatientBirthDate),
        rsNiftiExtendendHeaderInformationListCreateEntry("PatientSex",                    get, setString, setString,   formatString, &info->PatientSex),
        rsNiftiExtendendHeaderInformationListCreateEntry("PatientAge",                    get, setString, setString,   formatString, &info->PatientAge),
        rsNiftiExtendendHeaderInformationListCreateEntry("PatientWeight",                 get, setString, setString,   formatString, &info->PatientWeight),
        rsNiftiExtendendHeaderInformationListCreateEntry("SequenceName",                  get, setString, setString,   formatString, &info->SequenceName),
        rsNiftiExtendendHeaderInformationListCreateEntry("SliceThickness",                get, setString, setString,   formatString, &info->SliceThickness),
        rsNiftiExtendendHeaderInformationListCreateEntry("RepetitionTime",                get, setString, setString,   formatString, &info->RepetitionTime),
        rsNiftiExtendendHeaderInformationListCreateEntry("EchoTime",                      get, setString, setString,   formatString, &info->EchoTime),
        rsNiftiExtendendHeaderInformationListCreateEntry("MagneticFieldStrength",         get, setString, setString,   formatString, &info->MagneticFieldStrength),
        rsNiftiExtendendHeaderInformationListCreateEntry("SpacingBetweenSlices",          get, setString, setString,   formatString, &info->SpacingBetweenSlices),
        rsNiftiExtendendHeaderInformationListCreateEntry("NumberOfPhaseEncodingSteps",    get, setString, setString,   formatString, &info->NumberOfPhaseEncodingSteps),
        rsNiftiExtendendHeaderInformationListCreateEntry("PixelBandwidth",                get, setString, setString,   formatString, &info->PixelBandwidth),
        rsNiftiExtendendHeaderInformationListCreateEntry("SoftwareVersions",              get, setString, setString,   formatString, &info->SoftwareVersions),
        rsNiftiExtendendHeaderInformationListCreateEntry("ProtocolName",                  get, setString, setString,   formatString, &info->ProtocolName),
        rsNiftiExtendendHeaderInformationListCreateEntry("TransmitCoilName",              get, setString, setString,   formatString, &info->TransmitCoilName),
        rsNiftiExtendendHeaderInformationListCreateEntry("InPlanePhaseEncodingDirection", get, setString, setString,   formatString, &info->InPlanePhaseEncodingDirection),
        rsNiftiExtendendHeaderInformationListCreateEntry("PatientPosition",               get, setString, setString,   formatString, &info->PatientPosition),
        rsNiftiExtendendHeaderInformationListCreateEntry("BandwidthPerPixelPhaseEncode",  get, setDouble, parseDouble, formatDouble, &info->BandwidthPerPixelPhaseEncode),
        rsNiftiExtendendHeaderInformationListCreateEntry("SeriesNumber",                  get, setString, setString,   formatString, &info->SeriesNumber),
        rsNiftiExtendendHeaderInformationListCreateEntry("ImageComments",                 get, setString, setString,   formatString, &info->ImageComments),
        rsNiftiExtendendHeaderInformationListCreateEntry("MatrixSize",                    get, setString, setString,   formatString, &info->MatrixSize),
        rsNiftiExtendendHeaderInformationListCreateEntry("FieldOfView",                   get, setString, setString,   formatString, &info->FieldOfView),
        rsNiftiExtendendHeaderInformationListCreateEntry("GrappaFactor",                  get, setString, setString,   formatString, &info->GrappaFactor),
        rsNiftiExtendendHeaderInformationListCreateEntry("PhaseEncodingLines",            get, setString, setString,   formatString, &info->PhaseEncodingLines),
        rsNiftiExtendendHeaderInformationListCreateEntry("PhaseEncodingDirection",        get, setString, setString,   formatString, &info->PhaseEncodingDirection),
        rsNiftiExtendendHeaderInformationListCreateEntry("DwellTime",                     get, setDouble, parseDouble, formatDouble, &info->DwellTime),
        rsNiftiExtendendHeaderInformationListCreateEntry("Rows",                          get, setUShort, parseUShort, formatUShort, &info->Rows),
        rsNiftiExtendendHeaderInformationListCreateEntry("Columns",                       get, setUShort, parseUShort, formatUShort, &info->Columns),
        rsNiftiExtendendHeaderInformationListCreateEntry("MosaicRefAcqTimes",             get, setSliceTimings, parseSliceTimings, formatSliceTimings, &info->MosaicRefAcqTimes[0]),
        NULL
    };

    size_t nEntries = sizeof(entries) / sizeof(double*);

    rsNiftiExtendedHeaderInformationEntry** finalEntries = (rsNiftiExtendedHeaderInformationEntry**)rsMalloc(nEntries * sizeof(rsNiftiExtendedHeaderInformationEntry*));

    for (int i=0; i<nEntries; i++) {
        finalEntries[i] = entries[i];
    }

    return finalEntries;
}

void rsNiftiExtendendHeaderInformationListDestroyEntryMap(rsNiftiExtendedHeaderInformationEntry** entryMap) {
    for (int i=0; entryMap[i] != NULL; i++) {
        rsFree(entryMap[i]);
    }
    rsFree(entryMap);
}

void rsNiftiPrintExtendedHeaderInformation(rsNiftiExtendedHeaderInformation* info) {
    rsNiftiExtendedHeaderInformationEntry **entries = rsNiftiExtendendHeaderInformationListCreateEntryMap(info);

    for (int i=0; entries[i] != NULL; i++) {
        char *value = entries[i]->format(entries[i], TRUE);
        fprintf(stdout, "% 30s : %s\n", entries[i]->key, value);
        rsFree(value);
    }

    rsNiftiExtendendHeaderInformationListDestroyEntryMap(entries);
}

void rsNiftiCopyTrimmedValue(char *dest, char*buffer, size_t length)
{
    // copy raw value into a buffer
    char *buf = (char*)rsMalloc((length+1)*sizeof(char));
    memcpy(buf, buffer, length);

    // null-terminate it
    buf[length] = '\0';

    // copy trimmed value to destination
    sprintf(dest, "%s", rsTrimString(buf));

    rsFree(buf);
}

void rsNiftiAddExtendedHeaderInformation(rsNiftiExtendedHeaderInformation* info, const rsDicomElement* dicomElement, char* buffer, size_t length) {
    switch (dicomElement->tagGroup) {
        case 0x0002:
            if (dicomElement->tagElement == 0x0013 && length < sizeof(info->ImplementationVersionName)) {
                rsNiftiCopyTrimmedValue(&info->ImplementationVersionName[0], buffer, length);
            }
            break;
        case 0x0008:
            switch (dicomElement->tagElement) {
                case 0x0022:
                    if (length < sizeof(info->AcquisitionDate))
                        rsNiftiCopyTrimmedValue(&info->AcquisitionDate[0], buffer, length);
                    break;
                case 0x0080:
                    if (length < sizeof(info->InstitutionName))
                        rsNiftiCopyTrimmedValue(&info->InstitutionName[0], buffer, length);
                    break;
                case 0x0081:
                    if (length < sizeof(info->InstitutionAddress))
                        rsNiftiCopyTrimmedValue(&info->InstitutionAddress[0], buffer, length);
                    break;
                case 0x1030:
                    if (length < sizeof(info->StudyDescription))
                        rsNiftiCopyTrimmedValue(&info->StudyDescription[0], buffer, length);
                    break;
                case 0x103e:
                    if (length < sizeof(info->SeriesDescription))
                        rsNiftiCopyTrimmedValue(&info->SeriesDescription[0], buffer, length);
                    break;
                case 0x1070:
                    if (length < sizeof(info->OperatorsName))
                        rsNiftiCopyTrimmedValue(&info->OperatorsName[0], buffer, length);
                    break;
                case 0x1090:
                    if (length < sizeof(info->ManufacturerModelName))
                        rsNiftiCopyTrimmedValue(&info->ManufacturerModelName[0], buffer, length);
                    break;
                default:
                    break;
            }
            break;
        case 0x0010:
            switch (dicomElement->tagElement) {
                case 0x0010:
                    if (length < sizeof(info->PatientName))
                        rsNiftiCopyTrimmedValue(&info->PatientName[0], buffer, length);
                    break;
                case 0x0020:
                    if (length < sizeof(info->PatientID))
                        rsNiftiCopyTrimmedValue(&info->PatientID[0], buffer, length);
                    break;
                case 0x0030:
                    if (length < sizeof(info->PatientBirthDate))
                        rsNiftiCopyTrimmedValue(&info->PatientBirthDate[0], buffer, length);
                    break;
                case 0x0040:
                    if (length < sizeof(info->PatientSex))
                        rsNiftiCopyTrimmedValue(&info->PatientSex[0], buffer, length);
                    break;
                case 0x1010:
                    if (length < sizeof(info->PatientAge))
                        rsNiftiCopyTrimmedValue(&info->PatientAge[0], buffer, length);
                    break;
                case 0x1030:
                    if (length < sizeof(info->PatientWeight))
                        rsNiftiCopyTrimmedValue(&info->PatientWeight[0], buffer, length);
                    break;
                default:
                    break;
            }
            break;
        case 0x0018:
            switch (dicomElement->tagElement) {
                case 0x0024:
                    if (length < sizeof(info->SequenceName))
                        rsNiftiCopyTrimmedValue(&info->SequenceName[0], buffer, length);
                    break;
                case 0x0050:
                    if (length < sizeof(info->SliceThickness))
                        rsNiftiCopyTrimmedValue(&info->SliceThickness[0], buffer, length);
                    break;
                case 0x0080:
                    if (length < sizeof(info->RepetitionTime))
                        rsNiftiCopyTrimmedValue(&info->RepetitionTime[0], buffer, length);
                    break;
                case 0x0081:
                    if (length < sizeof(info->EchoTime))
                        rsNiftiCopyTrimmedValue(&info->EchoTime[0], buffer, length);
                    break;
                case 0x0087:
                    if (length < sizeof(info->MagneticFieldStrength))
                        rsNiftiCopyTrimmedValue(&info->MagneticFieldStrength[0], buffer, length);
                    break;
                case 0x0088:
                    if (length < sizeof(info->SpacingBetweenSlices))
                        rsNiftiCopyTrimmedValue(&info->SpacingBetweenSlices[0], buffer, length);
                    break;
                case 0x0089:
                    if (length < sizeof(info->NumberOfPhaseEncodingSteps))
                        rsNiftiCopyTrimmedValue(&info->NumberOfPhaseEncodingSteps[0], buffer, length);
                    break;
                case 0x0095:
                    if (length < sizeof(info->PixelBandwidth))
                        rsNiftiCopyTrimmedValue(&info->PixelBandwidth[0], buffer, length);
                    break;
                case 0x1020:
                    if (length < sizeof(info->SoftwareVersions))
                        rsNiftiCopyTrimmedValue(&info->SoftwareVersions[0], buffer, length);
                    break;
                case 0x1030:
                    if (length < sizeof(info->ProtocolName))
                        rsNiftiCopyTrimmedValue(&info->ProtocolName[0], buffer, length);
                    break;
                case 0x1251:
                    if (length < sizeof(info->TransmitCoilName))
                        rsNiftiCopyTrimmedValue(&info->TransmitCoilName[0], buffer, length);
                    break;
                case 0x1312:
                    if (length < sizeof(info->InPlanePhaseEncodingDirection))
                        rsNiftiCopyTrimmedValue(&info->InPlanePhaseEncodingDirection[0], buffer, length);
                    break;
                case 0x5100:
                    if (length < sizeof(info->PatientPosition))
                        rsNiftiCopyTrimmedValue(&info->PatientPosition[0], buffer, length);
                    break;
                default:
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
                for (short i=0; i<maxNumberOfSlices; i++) {
                    info->MosaicRefAcqTimes[i] = log(-1.0);
                }
                for (short i=0; i<applicableSlices; i++) {
                    info->MosaicRefAcqTimes[i] = sliceTimings[i];
                }
            }
            break;
        case 0x0020:
            if (dicomElement->tagElement == 0x0011 && length < sizeof(info->SeriesNumber)) {
                rsNiftiCopyTrimmedValue(&info->SeriesNumber[0], buffer, length);
            } else if (dicomElement->tagElement == 0x4000) {
                rsNiftiCopyTrimmedValue(&info->ImageComments[0], buffer, (size_t)fmin(length, sizeof(info->ImageComments)-1)); // truncate
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
                rsNiftiCopyTrimmedValue(&info->MatrixSize[0], buffer, length);
            } else if (dicomElement->tagElement == 0x100c && length < sizeof(info->FieldOfView)) {
                rsNiftiCopyTrimmedValue(&info->FieldOfView[0], buffer, length);
            }
            break;
        default:
            break;
    }
}

void *rsFindSubstring(const void *haystack, size_t hlen, const void *needle, size_t nlen)
{
    int needle_first;
    const void *p = haystack;
    size_t plen = hlen;

    if (!nlen)
        return NULL;

    needle_first = *(unsigned char *)needle;

    while (plen >= nlen && (p = memchr(p, needle_first, plen - nlen + 1)))
    {
        if (!memcmp(p, needle, nlen))
            return (void *)p;

        p++;
        plen = hlen - (p - haystack);
    }

    return NULL;
}

void rsNiftiReadTagFromSiemensExtraInformationHeader(char* value, const char* buffer, const size_t bufferLength, const char* entryTagName, const size_t maxExpectedLength)
{
    // find the siemens extra info awesomeness
    const char *startToken = "### ASCCONV BEGIN";
    const char *endToken = "### ASCCONV END";
    const char *seiStart = rsFindSubstring(&buffer[0], bufferLength, startToken, strlen(startToken));
    if (!seiStart) return;

    size_t remainingBytes = bufferLength-(seiStart-buffer);
    const char *seiEnd = rsFindSubstring(&seiStart[0], remainingBytes, endToken, strlen(endToken));
    if (!seiEnd) return;

    // if we got here we've found a proper siemens extra info definition, so let's search for
    // the info we're interested in, e.g.
    // entry.tag.name                        = value

    char *entryStart = rsFindSubstring(&seiStart[0], remainingBytes, entryTagName, strlen(entryTagName));
    if (!entryStart) return;

    // find '='
    remainingBytes = bufferLength-(entryStart-buffer);
    entryStart = rsFindSubstring(&entryStart[0], remainingBytes, "=", strlen("="));
    if (!entryStart) return;

    // forward to the number/value
    entryStart += 2;

    // find the end of it
    const char *entryEnd = rsFindSubstring(&entryStart[0], remainingBytes, "\n", strlen("\n"));
    if (!entryEnd) return;

    // copy the result
    const size_t entryLength = entryEnd - entryStart;
    const size_t resultLength = (size_t)fmin(entryLength, maxExpectedLength);
    rsNiftiCopyTrimmedValue(&value[0], entryStart, resultLength);
}

rsNiftiExtendedHeaderInformation* rsNiftiFindExtendedHeaderInformation(nifti_image *nim)
{
    if( nim->num_ext <= 0 || nim->ext_list == NULL ){
        return NULL;
    }

    nifti1_extension *ext = nim->ext_list;

    for ( int c = 0; c < nim->num_ext; c++ ){
        if ( ext->ecode == NIFTI_ECODE_JIMDIMINFO && ext->edata != NULL ) {
            return (rsNiftiExtendedHeaderInformation *) ext->edata;
        }
        ext++;
    }

    return NULL;
}

size_t rsNiftiGetDicomValueLength(const char *const valueRepresentation)
{
    const char longValueRepresentaions[][2] = {
        "OB", "OW", "OF", "SQ", "UT", "UN"
    };

    size_t length = 6;

    for ( size_t i=0; i<length; i++ ) {
        const char* currentRepresentation = longValueRepresentaions[i];

        if ( currentRepresentation[0] == valueRepresentation[0] && currentRepresentation[1] == valueRepresentation[1] ) {
            return 4;
        }
    }

    return 2;
}
