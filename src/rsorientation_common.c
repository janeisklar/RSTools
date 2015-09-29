#include "nifti/rsniftiutils.h"
#include "nifti/headerinfo.h"
#include "rsorientation_common.h"
#include "utils/rsio.h"
#include "rsorientation_ui.h"

rsOrientationTransformations* rsOrientationDetermineRequiredTransformations(mat44 in, const char *orientationCode);
const char rsOrientationPrintTransformationsDim(const short dim);
const char rsOrientationPrintTransformationsSign(const short dim, const rsOrientationTransformations* t);
void rsOrientationPrintTransformations(FILE *stream, const rsOrientationTransformations* t);
gsl_matrix* rsOrientationCreateGSLCameraMatrix(const mat44 in);
gsl_matrix* rsOrientationApplyTransformations(const rsOrientationTransformations* t, const gsl_matrix* mat, const short xhIn, const short yhIn, const short zhIn);
mat44 rsOrientationCameraMatrixToMat44(const gsl_matrix* m);
Point3D* rsOrientationTransformDataPoint(const rsOrientationTransformations* t, const rsNiftiFile* file, const Point3D* pointIn);
rsNiftiExtendedHeaderInformation* rsOrientationAttachDICOMInfo(const char* dicompath, rsNiftiFile* nifti);

void rsOrientationInit(rsOrientationParameters *p)
{
    p->parametersValid = FALSE;

    /* verify accessibility of inputs/outputs */
    BOOL inputsReadable = rsCheckInputs((const char*[]){
        (const char*)p->inputpath,
        (const char*)p->dicompath,
        RSIO_LASTFILE
    });
    
    BOOL outputsWritable = rsCheckOutputs((const char*[]){
        (const char*)p->outputpath,
        RSIO_LASTFILE
    });
    
    if ( ! inputsReadable || ! outputsWritable ) {
        return;
    }
        
    /* open input file */
    p->input = rsOpenNiftiFile(p->inputpath, RSNIFTI_OPEN_READ);

    if ( ! p->input->readable ) {
        fprintf(stderr, "\nError: The nifti file that was supplied as an input (%s) could not be read.\n", p->inputpath);
        return;
    }
	
    /* output the most important parameters to the user */
    if ( p->verbose ) {
        fprintf(stdout, "Input file: %s\n", p->inputpath);
        fprintf(stdout, "DICOM file: %s\n", p->dicompath);
        fprintf(stdout, "Output file: %s\n", p->outputpath);
        fprintf(stdout, "Output orientation: %s\n", p->orientation);
        fprintf(stdout, "Input Dim: %d %d %d (%d Volumes)\n", p->input->xDim, p->input->yDim, p->input->zDim, p->input->vDim);
    }
	
    p->parametersValid = TRUE;
}

void rsOrientationRun(rsOrientationParameters *p)
{
    p->parametersValid = FALSE;

	// get voxel spacing
	float xSpacing, ySpacing, zSpacing, tr;
	FslGetVoxDim(p->input->fslio, &xSpacing, &ySpacing, &zSpacing, &tr);

    // create output file
    p->output = rsCloneNiftiFile(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, p->input->vDim);

    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the output (%s) could not be created.\n", p->outputpath);
        return;
    }

	// load sform matrix
    nifti_image* inputNifti = p->input->fslio->niftiptr;
	gsl_matrix* qform = rsOrientationCreateGSLCameraMatrix(inputNifti->qto_xyz);
	gsl_matrix* sform = rsOrientationCreateGSLCameraMatrix(inputNifti->sto_xyz);
	
	int qform_i, qform_j, qform_k, sform_i, sform_j, sform_k;
	nifti_mat44_to_orientation(inputNifti->qto_xyz, &qform_i, &qform_j, &qform_k);
	nifti_mat44_to_orientation(inputNifti->sto_xyz, &sform_i, &sform_j, &sform_k);

    const int sform_code = inputNifti->sform_code,
              qform_code = inputNifti->qform_code;
	
	rsOrientationTransformations* qformTransformations = rsOrientationDetermineRequiredTransformations(inputNifti->qto_xyz, p->orientation);
	rsOrientationTransformations* sformTransformations = rsOrientationDetermineRequiredTransformations(inputNifti->sto_xyz, p->orientation);
	rsOrientationTransformations* dataTransformations = sform_code!=NIFTI_XFORM_UNKNOWN ? sformTransformations : qformTransformations;
	
	if ( p->verbose ) {
		fprintf(stdout, "\nqform camera matrix:\n");
		rs_gsl_matrix_fprintf(stdout, qform, " %+.6f");
		fprintf(stdout, "x=%s, y=%s, z=%s\n", nifti_orientation_string(qform_i), nifti_orientation_string(qform_j), nifti_orientation_string(qform_k));
		fprintf(stdout, "New orientation: ");
		rsOrientationPrintTransformations(stdout, qformTransformations);
		fprintf(stdout, "\n");

		fprintf(stdout, "\nsform camera matrix:\n");
		rs_gsl_matrix_fprintf(stdout, sform, " %+.6f");
		fprintf(stdout, "x=%s, y=%s, z=%s\n", nifti_orientation_string(sform_i), nifti_orientation_string(sform_j), nifti_orientation_string(sform_k));
		fprintf(stdout, "New orientation: ");
		rsOrientationPrintTransformations(stdout, sformTransformations);
		fprintf(stdout, "\n");
	}

    // prepare output file
    short dimsIn[3]  = {p->input->xDim, p->input->yDim, p->input->zDim};
    short dimsOut[3] = {dimsIn[dataTransformations->xDim], dimsIn[dataTransformations->yDim], dimsIn[dataTransformations->zDim]};
    p->output = rsCloneNiftiFile(p->outputpath, p->input, RSNIFTI_OPEN_ALLOC, 0);

    if ( ! p->output->readable ) {
        fprintf(stderr, "\nError: The nifti file containing the output (%s) could not be created.\n", p->outputpath);
        return;
    }

    p->output->xDim = dimsOut[0];
    p->output->yDim = dimsOut[1];
    p->output->zDim = dimsOut[2];
    FslSetDim(p->output->fslio, p->output->xDim, p->output->yDim, p->output->zDim, p->output->vDim);

    // transform camera matrices
    gsl_matrix *qform_transformed = rsOrientationApplyTransformations(qformTransformations, qform, p->input->xDim, p->input->yDim, p->input->zDim);
    gsl_matrix *sform_transformed = rsOrientationApplyTransformations(sformTransformations, sform, p->input->xDim, p->input->yDim, p->input->zDim);

    if ( p->verbose ) {
        fprintf(stdout, "\nnew qform camera matrix:\n");
        rs_gsl_matrix_fprintf(stdout, qform_transformed, " %+.6f");

        fprintf(stdout, "\nnew sform camera matrix:\n");
        rs_gsl_matrix_fprintf(stdout, sform_transformed, " %+.6f");
    }

    mat44 qform_mat = rsOrientationCameraMatrixToMat44(qform_transformed);
    mat44 sform_mat = rsOrientationCameraMatrixToMat44(qform_transformed);
    FslSetRigidXform(p->output->fslio, qform_code, qform_mat);
    FslSetStdXform(p->output->fslio, sform_code, sform_mat);
    float voxDimsIn[4];
    FslGetVoxDim(p->input->fslio,  &voxDimsIn[0], &voxDimsIn[1], &voxDimsIn[2], &voxDimsIn[3]);
    float voxDimsOut[4] = {voxDimsIn[dataTransformations->xDim], voxDimsIn[dataTransformations->yDim], voxDimsIn[dataTransformations->zDim], voxDimsIn[3]};
    FslSetVoxDim(p->output->fslio, voxDimsOut[0], voxDimsOut[1], voxDimsOut[2], voxDimsOut[3]);

    // attach DICOM header information to the nifti
    if (p->dicompath) {
        rsNiftiExtendedHeaderInformation *info = rsOrientationAttachDICOMInfo(p->dicompath, p->output);
        if (p->verbose) {
            fprintf(stdout, "Attached the following header information to the nifti file:\n");
            rsNiftiPrintExtendedHeaderInformation(info);
        }
    }

    // write nifti header
	rsWriteNiftiHeader(p->output->fslio, p->callString);

    // prepare the output file's content
    unsigned short x,y,z;
    Point3D *pointIn;
    Point3D *pointOut;

    for (z=0; z<p->input->zDim; z++) {
        for (y = 0; y < p->input->yDim; y++) {
            for (x = 0; x < p->input->xDim; x++) {

                pointIn  = rsMakePoint3D(x, y, z);
                pointOut = rsOrientationTransformDataPoint(dataTransformations, p->input, pointIn);

                rsCopyTimecourseFromInBufferToOutBuffer(
                    p->input->dt,     // datatype
                    p->output->data,  // output specification
                    pointOut,
                    p->output->xDim,
                    p->output->yDim,
                    p->output->zDim,
                    p->output->vDim,
                    p->input->data,   // input specification
                    pointIn,
                    p->input->xDim,
                    p->input->yDim,
                    p->input->zDim
                );

                rsFree(pointIn);
                rsFree(pointOut);
            }
        }
    }

    // write transformed file
    FslWriteVolumes(p->output->fslio, p->output->data, p->output->vDim);

	// cleanup
	gsl_matrix_free(sform);
	gsl_matrix_free(qform);

    p->parametersValid = TRUE;
}

void rsOrientationDestroy(rsOrientationParameters *p)
{
    if ( p->input != NULL ) {
        rsCloseNiftiFileAndFree(p->input);
        p->input = NULL;
    }
    
    if ( p->output != NULL ) {
        rsCloseNiftiFileAndFree(p->output);
        p->output = NULL;
    }

    rsOrientationFreeParams(p);
}

#pragma mark helper functions: orientation

gsl_matrix* rsOrientationCreateGSLCameraMatrix(const mat44 in)
{
	gsl_matrix* camera = gsl_matrix_calloc(4, 4);
	
	for (unsigned short i=0; i<4; i++) {
		for (unsigned short j=0; j<4; j++) {
			gsl_matrix_set(camera, i, j, in.m[i][j]);
		}
	}
	
	return camera;
}

mat44 rsOrientationCameraMatrixToMat44(const gsl_matrix* m)
{
    mat44 mat;
    for (unsigned short i=0; i<4; i++) {
        for (unsigned short j=0; j<4; j++) {
            mat.m[i][j] = gsl_matrix_get(m, i, j);
        }
    }
    return mat;
}

gsl_matrix* rsOrientationApplyTransformations(const rsOrientationTransformations* t, const gsl_matrix* mat, const short xhIn, const short yhIn, const short zhIn)
{
    const BOOL debug = FALSE;

    gsl_matrix *tmp = gsl_matrix_calloc(4, 4);
    gsl_matrix* result = gsl_matrix_alloc(4, 4);
    gsl_matrix_memcpy(result, mat);

    // create swap matrix
    gsl_matrix *swapMat = gsl_matrix_calloc(4, 4);
    const unsigned short dimMap[4] = {t->xDim, t->yDim, t->zDim, 3};
    for (unsigned short i=0; i<3; i++) {
        gsl_matrix_set(swapMat, i, dimMap[i], 1);
    }
    gsl_matrix_set(swapMat, 3, 3, 1);

    // create flip matrix
    const BOOL flipRequired[3] = {t->xFlip, t->yFlip, t->zFlip};
    const short sizesIn[3]  = {xhIn, yhIn, zhIn};
    const short sizesOut[3] = {sizesIn[dimMap[0]], sizesIn[dimMap[1]], sizesIn[dimMap[2]]};
    gsl_matrix *flipMat = gsl_matrix_calloc(4, 4);
    gsl_matrix_set_identity(flipMat);

    for (unsigned short dim=0; dim<3; dim++) {
        if (flipRequired[dim]) {
            gsl_matrix_set(flipMat, dim, dim, -1.0); // mirror
            gsl_matrix_set(flipMat, dim, 3, sizesOut[dim]-1); // correct translation
        }
    }

    // create combined transformation matrix
    int signum;
    gsl_matrix *transMat = gsl_matrix_calloc(4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, flipMat, swapMat, 0.0, transMat); // transMat <- flipMat * swapMat
    gsl_matrix *invTransMat = gsl_matrix_calloc(4, 4); // invTransMat <- transMat^-1
    gsl_permutation * perm = gsl_permutation_alloc(4);
    gsl_linalg_LU_decomp(transMat, perm, &signum);
    gsl_linalg_LU_invert(transMat, perm, invTransMat);
    gsl_permutation_free(perm);
    gsl_matrix_free(transMat);
    gsl_matrix_free(flipMat);
    gsl_matrix_free(swapMat);

    // apply transformation
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, result, invTransMat, 0.0, tmp); // result <- mat * invTransMat
    gsl_matrix_memcpy(result, tmp);

    gsl_matrix_free(invTransMat);
    gsl_matrix_free(tmp);

    return result;
}

rsOrientationTransformations* rsOrientationDetermineRequiredTransformations(mat44 in, const char *orientationCode)
{
	int orientations[3];
	rsOrientationTransformations* transformations = (rsOrientationTransformations*) rsMalloc(sizeof(rsOrientationTransformations));
	nifti_mat44_to_orientation(in, &orientations[0], &orientations[1], &orientations[2]);

    short lrDim, paDim, siDim;
    BOOL lrFlipped, paFlipped, siFlipped;
    for (short i=0; i<3; i++) {
        const char code = orientationCode[i];
        if (code=='L' || code=='R') {
            lrDim = i;
            lrFlipped = code=='R';
        } else if (code=='P' || code=='A') {
            paDim = i;
            paFlipped = code=='A';
        } else if (code=='S' || code=='I') {
            siDim = i;
            siFlipped = code=='I';
        }
    }

	for (short dim=0; dim<3; dim++) {
        short resultingDim;
        BOOL doFlip;
		switch (orientations[dim]) {
			case NIFTI_L2R:
                resultingDim = lrDim;
                doFlip = lrFlipped;
                break;
			case NIFTI_R2L:
                resultingDim = lrDim;
                doFlip = !lrFlipped;
				break;
			case NIFTI_P2A:
                resultingDim = paDim;
                doFlip = paFlipped;
                break;
			case NIFTI_A2P:
                resultingDim = paDim;
                doFlip = !paFlipped;
                break;
			case NIFTI_S2I:
                resultingDim = siDim;
                doFlip = siFlipped;
                break;
            case NIFTI_I2S:
                resultingDim = siDim;
                doFlip = !siFlipped;
                break;
            default:
                // we should never get here, let's do something nasty and divide by zero. that'll teach them!
                resultingDim =  1 / 0;
		}

        switch (resultingDim) {
            case 0:
                transformations->xDim = dim;
                transformations->xFlip = doFlip;
                break;
            case 1:
                transformations->yDim = dim;
                transformations->yFlip = doFlip;
                break;
            case 2:
                transformations->zDim = dim;
                transformations->zFlip = doFlip;
                break;
        }
	}
	
	return transformations;
}

Point3D* rsOrientationTransformDataPoint(const rsOrientationTransformations* t, const rsNiftiFile* file, const Point3D* pointIn)
{
    const unsigned int coordsIn[3]  = {pointIn->x, pointIn->y, pointIn->z};

    // swap dimensions
    unsigned int coordsOut[3] = {coordsIn[t->xDim], coordsIn[t->yDim], coordsIn[t->zDim]};

    // flip where necessary
    const BOOL flipRequired[3] = {t->xFlip, t->yFlip, t->zFlip};
    const short dimsIn[3]  = {file->xDim, file->yDim, file->zDim};
    const short dimsOut[3] = {dimsIn[t->xDim], dimsIn[t->yDim], dimsIn[t->zDim]};

    for (short dim=0; dim<3; dim++) {
        if (!flipRequired[dim]) {
            continue;
        }
        coordsOut[dim] = dimsOut[dim] - coordsOut[dim] - 1;
    }

    return rsMakePoint3D(coordsOut[0], coordsOut[1], coordsOut[2]);
}

void rsOrientationPrintTransformations(FILE *stream, const rsOrientationTransformations* t)
{
	const char xDim = rsOrientationPrintTransformationsDim(t->xDim);
	const char yDim = rsOrientationPrintTransformationsDim(t->yDim);
	const char zDim = rsOrientationPrintTransformationsDim(t->zDim);
	const char xSign = rsOrientationPrintTransformationsSign(0, t);
	const char ySign = rsOrientationPrintTransformationsSign(1, t);
	const char zSign = rsOrientationPrintTransformationsSign(2, t);
	
	fprintf(stream, "[%c%c,%c%c,%c%c]", xSign, xDim, ySign, yDim, zSign, zDim);
}

const char rsOrientationPrintTransformationsDim(const short dim) {
	switch (dim) {
		case 0:
			return 'x';
		case 1:
			return 'y';
		case 2:
			return 'z';
		default:
			return 'E';
	}
}

const char rsOrientationPrintTransformationsSign(const short dim, const rsOrientationTransformations* t) {
    switch (dim) {
    	case 0:
			return t->xFlip ? '-' : '\0';
    	case 1:
			return t->yFlip ? '-' : '\0';
		case 2:
			return t->zFlip ? '-' : '\0';
		default:
			return 'E';
    }
}

#pragma mark helper functions: dicom header

size_t rsOrientationGetDicomValueLength(const char valueRepresentation[2])
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

rsNiftiExtendedHeaderInformation* rsOrientationAttachDICOMInfo(const char* dicompath, rsNiftiFile* nifti)
{
    // open dicom
    FILE *f;
    f = fopen(dicompath, "r");

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

    rsNiftiExtendedHeaderInformation *info = rsNiftiInitializeExtendedHeaderInformation();

    // walk through dicom tags
    do {
        const rsDicomElement* currentElement = (rsDicomElement*)&buffer[i];

        // determine the length of the value length field
        size_t valueLength = rsOrientationGetDicomValueLength(currentElement->valueRepresentation);

        uint32_t length = 0;

        if ( valueLength == 2 ) {
            length = *((uint16_t*)(&(currentElement->valueRepresentation[0]) + 2));
        } else {
            length = *((uint32_t*)(&(currentElement->valueRepresentation[0]) + 4) );
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
            headerLength = i-offsetData;
            break;
        }

        // save field to nifti header if applicable
        rsNiftiAddExtendedHeaderInformation(info, currentElement, &buffer[i+offsetValue], lengthV);

        i+= offset;

    } while (i<size);

    // attach header info to nifti
    nifti_add_extension(nifti->fslio->niftiptr, (char*)info, sizeof(rsNiftiExtendedHeaderInformation), NIFTI_ECODE_JIMDIMINFO);

    // attach DICOM header to nifti
    if (headerLength > 0) {
        nifti_add_extension(nifti->fslio->niftiptr, &buffer[0], headerLength, NIFTI_ECODE_DICOM);
    }

    return info;
}