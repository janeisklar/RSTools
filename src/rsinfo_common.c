#include <src/nifti/rsniftiutils.h>
#include <src/nifti/headerinfo.h>
#include "rsinfo_common.h"
#include "utils/rsio.h"
#include "rsinfo_ui.h"

BOOL rsInfoPrintInfoForKey(rsNiftiExtendedHeaderInformation* info, const char* key);

void rsInfoInit(rsInfoParameters* p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->dicompath,
        RSIO_LASTFILE
    });
    
    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }    

    // open input file (header-only)
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_NONE);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as input (%s) could not be read.\n", p->inputpath);
        return;
    }

    // prepare dicom file if required
    if ( p->dicompath != NULL ) {
        p->dicom = fopen(p->dicompath, "w");
    }

    // enable showInfo and showComments if none specified
    if (p->dicom == NULL && !p->showComments && !p->showInfo && p->infoKey == NULL) {
        p->showComments = TRUE;
        p->showInfo = TRUE;
    }
    
    p->parametersValid = TRUE;    
}

void rsInfoRun(rsInfoParameters *p)
{
    p->parametersValid = FALSE;

    // check extensions
    nifti_image *nim = p->input->fslio->niftiptr;

    if( nim->num_ext <= 0 || nim->ext_list == NULL ){
        fprintf(stderr, "File does not contain any RSTools header information.\n");
        return;
    }

    // find extensions
    nifti1_extension *ext = nim->ext_list;
    nifti1_extension *commentExt = NULL;
    nifti1_extension *dicomExt = NULL;
    nifti1_extension *infoExt = NULL;

    for ( int c = 0; c < nim->num_ext; c++ ){
        if ( ext->ecode == NIFTI_ECODE_COMMENT && ext->edata != NULL ) {
            commentExt = ext;
        } else if ( ext->ecode == NIFTI_ECODE_DICOM && ext->edata != NULL ) {
            dicomExt = ext;
        } else if ( ext->ecode == NIFTI_ECODE_JIMDIMINFO && ext->edata != NULL ) {
            infoExt = ext;
        }
        ext++;
    }

    if ( commentExt == NULL && infoExt == NULL) {
        fprintf(stderr, "File does not contain any RSTools header information.\n");
        return;
    }

    // read out comment extension
    if (p->showComments) {
        if (commentExt == NULL) {
            fprintf(stderr, "File does not contain any comments");
            return;
        } else {
            const int size = commentExt->esize;
            char data[size + 1];
            strncpy(data, commentExt->edata, size);
            data[size] = '\0';
            fprintf(stdout, "File comments:\n%s\n\n", data);
        }
    }

    // read out info extension
    if (p->showInfo || p->infoKey != NULL) {
        if (infoExt == NULL) {
            fprintf(stderr, "File does not contain any extended header info");
            return;
        } else {
            const size_t size = infoExt->esize;
            if (size < sizeof(rsNiftiExtendedHeaderInformation)) {
                fprintf(stderr, "Found illegal extra header information!\n");
            } else if (p->infoKey != NULL) {
                rsNiftiExtendedHeaderInformation *info = (rsNiftiExtendedHeaderInformation *) infoExt->edata;
                BOOL success = rsInfoPrintInfoForKey(info, p->infoKey);
                if (!success) {
                    return;
                }
            } else {
                fprintf(stdout, "Extra header information:\n");
                rsNiftiExtendedHeaderInformation *info = (rsNiftiExtendedHeaderInformation *) infoExt->edata;
                rsNiftiPrintExtendedHeaderInformation(info);
                fprintf(stdout, "\n");
            }
        }
    }

    // write out dicom header if requested
    if (p->dicom != NULL) {
        if (dicomExt == NULL) {
            fprintf(stderr, "File does not contain a dicom header!");
            return;
        } else {
            const int size = dicomExt->esize;
            fwrite(dicomExt->edata, sizeof(char), size-8, p->dicom);
        }
    }

    p->parametersValid = TRUE;
}

void rsInfoDestroy(rsInfoParameters* p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFile(p->input, FALSE);
        p->input = NULL;
    }

    if ( p->dicom != NULL ) {
        fclose(p->dicom);
        p->input = NULL;
    }

    rsInfoFreeParams(p);
}

BOOL rsInfoPrintInfoForKey(rsNiftiExtendedHeaderInformation* info, const char* key)
{
    // make all caps-down
    char* k = rsString(key);
    for(short i = 0; i<strlen(key); i++){
        k[i] = tolower(k[i]);
    }

    if ( ! strcmp(k, "rstoolsheaderversion") ) {
        fprintf(stdout, "%d", (int)info->RSToolsHeaderVersion);
    } else if ( ! strcmp(k, "implementationversionname") ) {
        fprintf(stdout, "%s", info->ImplementationVersionName);
    } else if ( ! strcmp(k, "acquisitiondate") ) {
        fprintf(stdout, "%s", info->AcquisitionDate);
    } else if ( ! strcmp(k, "institutionname") ) {
        fprintf(stdout, "%s", info->InstitutionName);
    } else if ( ! strcmp(k, "institutionaddress") ) {
        fprintf(stdout, "%s", info->InstitutionAddress);
    } else if ( ! strcmp(k, "studydescription") ) {
        fprintf(stdout, "%s", info->StudyDescription);
    } else if ( ! strcmp(k, "seriesdescription") ) {
        fprintf(stdout, "%s", info->SeriesDescription);
    } else if ( ! strcmp(k, "operatorsname") ) {
        fprintf(stdout, "%s", info->OperatorsName);
    } else if ( ! strcmp(k, "manufacturermodelname") ) {
        fprintf(stdout, "%s", info->ManufacturerModelName);
    } else if ( ! strcmp(k, "patientname") ) {
        fprintf(stdout, "%s", info->PatientName);
    } else if ( ! strcmp(k, "patientid") ) {
        fprintf(stdout, "%s", info->PatientID);
    } else if ( ! strcmp(k, "patientbirthdate") ) {
        fprintf(stdout, "%s", info->PatientBirthDate);
    } else if ( ! strcmp(k, "patientsex") ) {
        fprintf(stdout, "%s", info->PatientSex);
    } else if ( ! strcmp(k, "patientage") ) {
        fprintf(stdout, "%s", info->PatientAge);
    } else if ( ! strcmp(k, "patientweight") ) {
        fprintf(stdout, "%s", info->PatientWeight);
    } else if ( ! strcmp(k, "sequencename") ) {
        fprintf(stdout, "%s", info->SequenceName);
    } else if ( ! strcmp(k, "slicethickness") ) {
        fprintf(stdout, "%s", info->SliceThickness);
    } else if ( ! strcmp(k, "repetitiontime") ) {
        fprintf(stdout, "%s", info->RepetitionTime);
    } else if ( ! strcmp(k, "echotime") ) {
        fprintf(stdout, "%s", info->EchoTime);
    } else if ( ! strcmp(k, "magneticfieldstrength") ) {
        fprintf(stdout, "%s", info->MagneticFieldStrength);
    } else if ( ! strcmp(k, "spacingbetweenslices") ) {
        fprintf(stdout, "%s", info->SpacingBetweenSlices);
    } else if ( ! strcmp(k, "numberofphaseencodingsteps") ) {
        fprintf(stdout, "%s", info->NumberOfPhaseEncodingSteps);
    } else if ( ! strcmp(k, "pixelbandwidth") ) {
        fprintf(stdout, "%s", info->PixelBandwidth);
    } else if ( ! strcmp(k, "softwareversions") ) {
        fprintf(stdout, "%s", info->SoftwareVersions);
    } else if ( ! strcmp(k, "protocolname") ) {
        fprintf(stdout, "%s", info->ProtocolName);
    } else if ( ! strcmp(k, "transmitcoilname") ) {
        fprintf(stdout, "%s", info->TransmitCoilName);
    } else if ( ! strcmp(k, "inplanephaseencodingdirection") ) {
        fprintf(stdout, "%s", info->InPlanePhaseEncodingDirection);
    } else if ( ! strcmp(k, "patientposition") ) {
        fprintf(stdout, "%s", info->PatientPosition);
    } else if ( ! strcmp(k, "bandwidthperpixelphaseencode") ) {
        fprintf(stdout, "%.15f", info->BandwidthPerPixelPhaseEncode);
    } else if ( ! strcmp(k, "seriesnumber") ) {
        fprintf(stdout, "%s", info->SeriesNumber);
    } else if ( ! strcmp(k, "imagecomments") ) {
        fprintf(stdout, "%s", info->ImageComments);
    } else if ( ! strcmp(k, "matrixsize") ) {
        fprintf(stdout, "%s", info->MatrixSize);
    } else if ( ! strcmp(k, "fieldofview") ) {
        fprintf(stdout, "%s", info->FieldOfView);
    } else if ( ! strcmp(k, "grappafactor") ) {
        fprintf(stdout, "%s", info->GrappaFactor);
    } else if ( ! strcmp(k, "phaseencodinglines") ) {
        fprintf(stdout, "%s", info->PhaseEncodingLines);
    } else if ( ! strcmp(k, "dwelltime") ) {
        fprintf(stdout, "%.15f", info->DwellTime);
    } else if ( ! strcmp(k, "rows") ) {
        fprintf(stdout, "%d", (int)info->Rows);
    } else if ( ! strcmp(k, "columns") ) {
        fprintf(stdout, "%d", (int)info->Columns);
    } else if ( ! strcmp(k, "mosaicrefacqtimes") ) {
        for (short i=0; (i<sizeof(info->MosaicRefAcqTimes)/sizeof(double) && !isnan(info->MosaicRefAcqTimes[i])); i++) {
            fprintf(stdout, "%.1f\n", info->MosaicRefAcqTimes[i]);
        }
    } else {
        return FALSE;
    }

    return TRUE;
}
