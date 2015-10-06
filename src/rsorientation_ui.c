#include "rsorientation_ui.h"

rsOrientationParameters *rsOrientationInitParameters() {
    rsOrientationParameters *p = (rsOrientationParameters*)rsMalloc(sizeof(rsOrientationParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->dicompath            = NULL;
    p->orientation          = NULL;
    p->phaseencdir          = NULL;
    p->callString           = NULL;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->dicom                = NULL;
    p->output               = NULL;
    p->parametersValid      = FALSE;
    
    return p;
}

void rsOrientationFreeParams(rsOrientationParameters *p) {
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->dicompath);
    rsFree(p->orientation);
    rsFree(p->phaseencdir);
    rsFree(p->input);
    rsFree(p->dicom);
    rsFree(p->output);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

rsOrientationParameters *rsOrientationParseParams(int argc, char * argv[]) {

    rsOrientationParameters *p = rsOrientationInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsOrientationBuildInterface(p);
    
    // parse
    BOOL parsingSuccessful = rsUIParse(p->interface, argc, argv, (void*)p);
    
    if ( ! parsingSuccessful ) {
        return p;
    }
    
    // check if the required arguments have been provided
    if ( p->inputpath == NULL ) {
        fprintf(stderr, "No input volume specified!\n");
        return p;
    }
    
    if ( p->outputpath == NULL ) {
        fprintf(stderr, "An output volume must be specified!\n");
        return p;
    }

    if (p->orientation == NULL) {
        p->orientation = "LPI";
    }

    BOOL validOrientation = strlen(p->orientation) == 3;

    // validate and normalize orientation code
    if (validOrientation) {

        BOOL LRpresent = FALSE;
        BOOL PApresent = FALSE;
        BOOL ISpresent = FALSE;

        // make all uppercase
        p->orientation = rsString(p->orientation);
        for(short i = 0; p->orientation[i]; i++){
            p->orientation[i] = toupper(p->orientation[i]);
        }

        // check if we got a letter for each direction
        for (short i=0; i<3; i++) {
            const char code = p->orientation[i];
            if (code=='L' || code=='R') {
                LRpresent = TRUE;
            } else if (code=='P' || code=='A') {
                PApresent = TRUE;
            } else if (code=='I' || code=='S') {
                ISpresent = TRUE;
            }
        }
        validOrientation = LRpresent && PApresent && ISpresent;
    }

    if (!validOrientation) {
        fprintf(stderr, "The parameter --orientation (-a) needs to be a three letter code containing the letters (L, R, P, A, I, S)!\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsOrientationBuildInterface(rsOrientationParameters *p)
{  
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Ensures an LPI orientation (unless specified differently) of the output nifti, optionally attaches the header of the DICOM to the nifti header and furthermore saves some properties from the Siemens DICOM header in the nifti which can be used for further processing (such as the TR, multiband-factor, multiband slice-timing information, grappa factor, bandwith, etc.). The attached header information are accesible through rsinfo.";

    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the 4D volume to processed";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "output";
    o->shorthand           = 'o';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->outputpath;
    o->cli_description     = "the volume in which the result will be saved";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
	
    o = rsUINewOption();
    o->name                = "dicom";
    o->shorthand           = 'd';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->dicompath;
    o->cli_description     = "the DICOM file from which the header information will be taken";
    o->cli_arg_description = "<dicom>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "orientation";
    o->shorthand           = 'a';
    o->type                = G_OPTION_ARG_STRING;
    o->storage             = &p->orientation;
    o->cli_description     = "the desired orientation of the output as a 3 letter-code, such as LPI. The following letters are possible: R(Right-to-left), L(Left-to-right), A(Anterior-to-posterior), P(Posterior-to-anterior, I(Inferior-to-superior), S(Superior-to-inferior)";
    o->cli_arg_description = "<ori-code>";
    rsUIAddOption(p->interface, o);

    o = rsUINewOption();
    o->name                = "phaseencdir";
    o->shorthand           = 'p';
    o->type                = G_OPTION_ARG_STRING;
    o->storage             = &p->phaseencdir;
    o->cli_description     = "the phase encoding direction of the input volume. If specified, the phase encdoing direction will be written to the resulting nifti output with respect to the new orientation and will be available using rsinfo. The expected format is a two letter string such as: y-, x+, etc.";
    o->cli_arg_description = "<direction>";
    rsUIAddOption(p->interface, o);
	
    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
}
