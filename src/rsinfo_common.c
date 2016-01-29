#include <src/nifti/rsniftiutils.h>
#include <src/nifti/headerinfo.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/fcntl.h>
#include <externals/fslio/fslio.h>
#include <nifti1_io.h>
#include "rsinfo_common.h"
#include "utils/rsio.h"
#include "rsinfo_ui.h"

BOOL rsInfoPrintInfoForKey(rsNiftiExtendedHeaderInformation* info, const char* key);
BOOL rsInfoParseHeaderModificationArgument(char **key, char **value, const char *arg);
BOOL rsInfoSetValueForKey(rsNiftiExtendedHeaderInformation* info, const char* key, const char* value);

void rsInfoInit(rsInfoParameters* p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        (const char*)p->extensionSource,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->dicompath,
        (const char*)p->outputpath,
        RSIO_LASTFILE
    });
    
    if (!inputsReadable || !outputsWritable) {
        return;
    }    

    // open input file (header-only)
    p->input = rsOpenNiftiFile(p->inputpath, p->outputpath==NULL ? RSNIFTI_OPEN_NONE : RSNIFTI_OPEN_READ);

    if (!p->input->readable) {
        fprintf(stderr, "\nError: The nifti file that was supplied as input (%s) could not be read.\n", p->inputpath);
        return;
    }

    // if an extension source is set the input file needs to be copied over to a new file

    // open extension source if required
    if (p->extensionSource != NULL) {
        p->extensionInput = rsOpenNiftiFile(p->extensionSource, RSNIFTI_OPEN_NONE);

        if (!p->extensionInput->readable) {
            fprintf(stderr, "\nError: The nifti file that was supplied as source for the extensions to be copied (%s) could not be read.\n", p->extensionSource);
            return;
        }
    }

    if (p->extensionSource != NULL && p->outputpath == NULL) {
        fprintf(stderr, "\nSpecifying an extension source requires an output paht to be set as the file will have to be re-created (header length varies).\n");
        return;
    }

    // prepare dicom file if required
    if (p->dicompath != NULL) {
        p->dicom = fopen(p->dicompath, "w");
    }

    // enable showInfo and showComments if none specified
    if (p->dicom == NULL && !p->showComments && !p->showInfo && p->infoKey == NULL) {
        p->showComments = TRUE;
        p->showInfo = TRUE;
    }
    
    p->parametersValid = TRUE;    
}

int rsNiftiWriteExtensions(znzFile fp, nifti_image *nim)
{
    nifti1_extension * list;
    char               extdr[4] = { 0, 0, 0, 0 };
    int                c, size, ok = 1;

    if( znz_isnull(fp) || !nim || nim->num_ext < 0 ){
        fprintf(stderr,"** nifti_write_extensions, bad params\n");
        return -1;
    }

    /* if invalid extension list, clear num_ext */
    if (!valid_nifti_extensions(nim)) {
        nim->num_ext = 0;
    }

    /* write out extender block */
    if (nim->num_ext > 0) {
        extdr[0] = 1;
    }

    if (nifti_write_buffer(fp, extdr, 4) != 4){
        fprintf(stderr,"** failed to write extender\n");
        return -1;
    }

    list = nim->ext_list;
    for ( c = 0; c < nim->num_ext; c++ ){
        size = (int)nifti_write_buffer(fp, &list->esize, sizeof(int));
        ok = (size == (int)sizeof(int));

        if (ok) {
            size = (int)nifti_write_buffer(fp, &list->ecode, sizeof(int));
            ok = (size == (int)sizeof(int));
        }

        if (ok) {
            size = (int)nifti_write_buffer(fp, list->edata, list->esize - 8);
            ok = (size == list->esize - 8);
        }

        if (!ok) {
            fprintf(stderr,"** failed while writing extension #%d\n",c);
            return -1;
        }

        list++;
    }

    return nim->num_ext;
}

void rsInfoRun(rsInfoParameters *p)
{
    p->parametersValid = FALSE;
    nifti_image *nim = p->input->fslio->niftiptr;

    BOOL modifyHeaderOnly = p->outputpath == NULL;
    BOOL showOutput = p->modArgs < 1 && p->extensionSource == NULL;

    // create output file if required
    if (p->outputpath != NULL) {
        rsNiftiFile *ref = p->extensionInput == NULL ? p->input : p->extensionInput;
        p->output = rsCloneNiftiFileWithNewDimensions(p->outputpath, ref, RSNIFTI_OPEN_NONE, p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
        if (!p->output->readable) {
            fprintf(stderr, "\nError: The nifti file containing the output (%s) could not be created.\n", p->outputpath);
            return;
        }
        // we copied the extension source to copy the extensions along with the header we therefore have to
        // make sure that all other nifti properties are copied over from the input nifti and not from the extension source
        p->output->intent_code = p->input->intent_code;
        p->output->intent_p1   = p->input->intent_p1;
        p->output->intent_p2   = p->input->intent_p2;
        p->output->intent_p3   = p->input->intent_p3;
        p->output->dt          = p->input->dt;
        p->output->pixtype     = p->input->pixtype;
        p->output->inter       = p->input->inter;
        p->output->slope       = p->input->slope;
        FslSetIntent(p->output->fslio, p->output->intent_code, p->output->intent_p1, p->output->intent_p2, p->output->intent_p3);
        size_t dim = 4;
        FslGetDimensionality(p->input->fslio, &dim);
        FslSetDimensionality(p->output->fslio, dim);
        FslSetDataType(p->output->fslio, p->output->dt);
        float tmp[4];
        FslGetVoxDim(p->input->fslio, &tmp[0], &tmp[1], &tmp[2], &tmp[3]);
        FslSetVoxDim(p->output->fslio, tmp[0], tmp[1], tmp[2], tmp[3]);
        mat44 tmpMat;
        FslGetStdXform(p->input->fslio, &tmpMat);
        FslSetStdXform(p->output->fslio, p->input->fslio->niftiptr->sform_code, tmpMat);
        FslGetRigidXform(p->input->fslio, &tmpMat);
        FslSetRigidXform(p->output->fslio, p->input->fslio->niftiptr->qform_code, tmpMat);
        short tmpDt;
        FslGetDataType(p->input->fslio, &tmpDt);
        FslSetDataType(p->output->fslio, tmpDt);
        nim = p->output->fslio->niftiptr;
    }

    // check extensions
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

    unsigned int numModArgs = p->modArgs == NULL ? 0 : g_strv_length(p->modArgs);

    // read out comment extension
    if (p->showComments && showOutput) {
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
    rsNiftiExtendedHeaderInformation *info;
    if (p->showInfo || p->infoKey != NULL || numModArgs > 0) {
        if (infoExt == NULL) {
            fprintf(stderr, "File does not contain any extended header info");
            return;
        } else {
            const size_t size = infoExt->esize;
            if (size < sizeof(rsNiftiExtendedHeaderInformation)) {
                fprintf(stderr, "Found illegal extra header information!\n");
            } else if (p->infoKey != NULL) {
                info = (rsNiftiExtendedHeaderInformation *) infoExt->edata;
                p->parametersValid = rsInfoPrintInfoForKey(info, p->infoKey);
                return;
            } else {
                info = (rsNiftiExtendedHeaderInformation *) infoExt->edata;
                if (showOutput) {
                    fprintf(stdout, "Extra header information:\n");
                    rsNiftiPrintExtendedHeaderInformation(info);
                    fprintf(stdout, "\n");
                }
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

    // modify header values if requested
    if (numModArgs > 0) {
        for(int i = 0; i< numModArgs; i++){
            const char *argument = p->modArgs[i];
            char *key;
            char *value;

            if (!rsInfoParseHeaderModificationArgument(&key, &value, argument)) {
                fprintf(stderr, "Failed to parse header modification argument '%s'. Expected notation is key=value.\n", argument);
                return;
            }

            if (!rsInfoSetValueForKey(info, key, value)) {
                fprintf(stderr, "Key '%s' is invalid!\n", key);
                return;
            }
        }
    }

    // write out modified header only
    if (modifyHeaderOnly) {

        // write out nifti header
        int imageType = FslGetFileType(p->input->fslio);
        znzFile output = znzopen(p->inputpath, "r+", FslIsCompressedFileType(imageType));

        if (znz_isnull(output)) {
            fprintf(stderr, "failed to open %s for writing\n", p->inputpath);
            return;
        }

        nifti_1_header nhdr = nifti_convert_nim2nhdr(p->input->fslio->niftiptr);
        size_t ss = znzwrite(&nhdr, 1, sizeof(nhdr), output);

        if (ss < sizeof(nhdr)) {
            fprintf(stderr, "failed writing header to %s\n", p->inputpath);
            znzclose(output);
            return;
        }

        // write out nifti header extensions
        ss = rsNiftiWriteExtensions(output, p->input->fslio->niftiptr);

        if (ss < p->input->fslio->niftiptr->num_ext) {
            fprintf(stderr, "failed writing all nifti extensions (%ld out of %ld) to %s\n", ss, p->input->fslio->niftiptr->num_ext, p->inputpath);
            znzclose(output);
            return;
        }

        znzclose(output);
    } else {
        // rewrite the complete file (outputpath)
        rsWriteNiftiHeader(p->output->fslio, NULL);
        FslWriteVolumes(p->output->fslio, p->input->data, p->output->vDim);
    }

    p->parametersValid = TRUE;
}

void rsInfoDestroy(rsInfoParameters* p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFile(p->input, FALSE);
        p->input = NULL;
    }

    if ( p->output != NULL ) {
        rsCloseNiftiFile(p->output, FALSE);
        p->output = NULL;
    }

    if ( p->extensionInput != NULL ) {
        rsCloseNiftiFile(p->extensionInput, TRUE);
        p->output = NULL;
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

    rsNiftiExtendedHeaderInformationEntry **entries = rsNiftiExtendendHeaderInformationListCreateEntryMap(info);

    // find key in header information map
    BOOL found = FALSE;
    for (int i=0; entries[i] != NULL; i++) {
        char* key = rsString(entries[i]->key);

        for (short j = 0; j<strlen(key); j++){
            key[j] = tolower(key[j]);
        }

        if (!strcmp(k, key)) {
            char *value = entries[i]->format(entries[i], FALSE);
            fprintf(stdout, "%s", value);
            rsFree(value);
            found = TRUE;
        }
        rsFree(key);

        if (found) {
            break;
        }
    }

    rsNiftiExtendendHeaderInformationListDestroyEntryMap(entries);
    return found;
}

BOOL rsInfoSetValueForKey(rsNiftiExtendedHeaderInformation* info, const char* key, const char* value)
{
    // make all caps-down
    char* k = rsString(key);
    for(short i = 0; i<strlen(key); i++){
        k[i] = tolower(k[i]);
    }

    rsNiftiExtendedHeaderInformationEntry **entries = rsNiftiExtendendHeaderInformationListCreateEntryMap(info);

    // find key in header information map
    BOOL found = FALSE;
    for (int i=0; entries[i] != NULL; i++) {
        char* key = rsString(entries[i]->key);

        for (short j = 0; j<strlen(key); j++){
            key[j] = tolower(key[j]);
        }

        if (!strcmp(k, key)) {
            entries[i]->parse(entries[i], value);
            found = TRUE;
        }
        rsFree(key);

        if (found) {
            break;
        }
    }

    rsNiftiExtendendHeaderInformationListDestroyEntryMap(entries);
    return found;
}

BOOL rsInfoParseHeaderModificationArgument(char **key, char **value, const char *arg)
{
    char *argCopy = rsString(arg);

    const char* keyPtr = strtok(argCopy, "=");
    const char* valPtr = strtok(NULL, "=");

    BOOL success;

    if ( strtok(NULL,"=") != NULL || keyPtr == NULL || valPtr == NULL || strlen(keyPtr) < 1) {
        success = FALSE;
    } else {
        *key = rsString(keyPtr);
        *value = rsString(valPtr);
        success = TRUE;
    }

    rsFree(argCopy);
    return success;
}