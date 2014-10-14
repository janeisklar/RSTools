#include "rstimecourse_ui.h"

rsTimecourseParameters* rsTimecourseInitParameters()
{
    rsTimecourseParameters *p = (rsTimecourseParameters*)rsMalloc(sizeof(rsTimecourseParameters));
    
    p->inputpath            = NULL;
    p->outputpath           = NULL;
    p->maskpath             = NULL;
    p->mask2path            = NULL;
    p->savemaskpath         = NULL;
    p->eigenvaluespath      = NULL;
    p->spatialmappath       = NULL;
    p->algorithm            = RSTOOLS_TIMECOURSE_ALGORITHM_MEAN;
    p->useStandardScores    = FALSE;
    p->minVariance          = 1.0;
    p->nComponents          = -1;
    p->input                = NULL;
    p->output               = NULL;
    p->mask                 = NULL;
    p->mask2                = NULL;
    p->point                = NULL;

    p->parametersValid      = FALSE;
    p->verbose              = FALSE;
    p->interface            = NULL;
    p->callString           = NULL;
    p->threads              = 1;

    return p;
}

rsTimecourseParameters* rsTimecourseParseParams(int argc, char * argv[])
{
    rsTimecourseParameters *p = rsTimecourseInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    rsTimecourseBuildInterface(p);

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

    if ( p->maskpath == NULL && p->point == NULL ) {
        fprintf(stderr, "Either a binary mask or a voxel coordinate must be specified!\n");
        return p;
    }

    if ( p->mask2path == NULL && p->algorithm == RSTOOLS_TIMECOURSE_ALGORITHM_CSP ) {
        fprintf(stderr, "A second binary mask must be supplied to run CSP! (use --mask2)\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsTimecourseBuildInterface(rsTimecourseParameters *p)
{
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Given a 4D-Nifti, this tool extracts the time course for a single voxel or aggregates the timecourse in a mask using different algorithms";
    //p->interface->helpIndent    = 33;
    
    GOptionArgFunc cbAlgorithm = (GOptionArgFunc)rsTimecourseParseAlgorithm;
    GOptionArgFunc cbPoint     = (GOptionArgFunc)rsTimecourseParsePoint;
    
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the volume from which the timecourse is to be extracted";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "output";
    o->shorthand           = 'o';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->outputpath;
    o->cli_description     = "file in which the resulting timecourse will be saved. If omitted the result will be printed directly to stdout";
    o->gui_description     = "file in which the resulting timecourse will be saved";
    o->cli_arg_description = "<*.txt>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "algorithm";
    o->shorthand           = 'a';
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbAlgorithm;
    o->cli_description     = "the algorithm used to aggregate the data within the specified mask(s), e.g. mean, stddev, spca, tpca or csp";
    o->cli_arg_description = "<algo>";
    o->defaultValue        = "mean";
    rsUIOptionValue allowedValues[] = {
      {"mean",   "for every volume in the 4D input file the intensity values in the mask region are meaned, resulting in a meaned timecourse for that region"},
      {"stddev", "like 'mean', but instead of meaning the mask region, the standard deviation is computed"},
      {"tpca",   "performs a PCA on the temporal dimension of the 4D input in the mask region"},
      {"spca",   "performs a PCA on the spatial dimension of the 4D input limited to voxels in the mask region"},
      {"csp",    "performs the common spatial pattern procedure on the 4D input in the mask region for both masks. The first half of the returned components will maximize the variance for mask1 while the variance for mask2 is minimal. The second half will be exactly the opposite (var(mask1) minimal, var(mask2) maximal)"},
      NULL
    };
    rsUISetOptionValues(o, allowedValues);
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "mask";
    o->shorthand           = 'm';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->maskpath;
    o->cli_description     = "a mask specifying the region that the algorithm is perforned on";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "mask2";
    o->shorthand           = 'n';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->mask2path;
    o->cli_description     = "(use only with the csp algorithm) a second mask specifying the region for the second condition";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "point";
    o->shorthand           = 'p';
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbPoint;
    o->cli_description     = "speficies a voxel using nifti coordinates(0-based) from which the timecourse is to be extracted";
    o->cli_arg_description = "<x>,<y>,<z>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "retainVariance";
    o->shorthand           = 'V';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->minVariance;
    o->cli_description     = "(use only with pca) percentage of explained variance that will be retained, e.g. '0.4'. keep in mind that a higher percentage will result in more components that are to be returned.";
    o->cli_arg_description = "<float>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "retainComponents";
    o->shorthand           = 'c';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->nComponents;
    o->cli_description     = "(use only with pca or csp) number of PCA/CSP components that will be outputted. This should be a multiple of two when running CSP as the components will be distributed equally over both masks.";
    o->cli_arg_description = "<n>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "useStandardScores";
    o->shorthand           = 's';
    o->storage             = &p->useStandardScores;
    o->cli_description     = "(use only with pca) remove mean and set std. dev to 1 prior to running the pca";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "spatialMap";
    o->shorthand           = 'M';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->spatialmappath;
    o->cli_description     = "(use only with pca) store the spatial map that is created using the PCA components(4D nifti with one volume for every component)";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
        
    o = rsUINewOption();
    o->name                = "eigenvalues";
    o->shorthand           = 'e';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->eigenvaluespath;
    o->cli_description     = "(use only with pca) write out all eigenvalues to the file that is specified with this option";
    o->cli_arg_description = "<*.txt>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "threads";
    o->shorthand           = 't';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->threads;
    o->cli_description     = "number of threads used for processing";
    o->cli_arg_description = "<n>";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    o->showInGUI           = FALSE;
    rsUIAddOption(p->interface, o);
       
    // initialize the more advanced and rather unusual options   
    o = rsUINewOption();
    o->name                = "savemask";
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->savemaskpath;
    o->cli_description     = "optional path where the rescaled mask specified with -mask will be saved. The saved file with have the same dimensions as the input volume.";
    o->cli_arg_description = "<volume>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
}

void rsTimecourseFreeParams(rsTimecourseParameters* p)
{
    rsFree(p->inputpath);
    rsFree(p->outputpath);
    rsFree(p->maskpath);
    rsFree(p->mask2path);
    rsFree(p->savemaskpath);
    rsFree(p->eigenvaluespath);
    rsFree(p->spatialmappath);
    rsFree(p->input);
    rsFree(p->callString);
    rsUIDestroyInterface(p->interface);
    rsFree(p);
}

gboolean rsTimecourseParsePoint(const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
    rsTimecourseParameters *p = (rsTimecourseParameters*) data;

    // copy value(const)
    size_t length = strlen(value);
    char v[length+1];
    sprintf(&v[0], "%s", value);

    // parse value
    char *strX;
    char *strY;
    char *strZ;
    strX = strtok(   v, ",");
    strY = strtok(NULL, ",");
    strZ = strtok(NULL, ",");

    int x=-1, y=-1, z=-1;

    // if we were given exactly 3 numbers separated by comma parse them 
    if ( strtok(NULL,",") == NULL && strX != NULL && strY != NULL && strZ != NULL ) {
        x = atoi(strX);
        y = atoi(strY);
        z = atoi(strZ);
    }
    
    // return success if we sucessfully received 3 numbers
    if ( x >= 0 && y >= 0 && z >= 0 ) {
        p->point = rsMakePoint3D(x,y,z);
        return TRUE;
    }
    
    // anything else should lead to an error
    g_set_error(
        error,
        G_OPTION_ERROR,
        G_OPTION_ERROR_BAD_VALUE,
        "%s: %s",
        option_name,
        "format should be 'x,y,z'"
    );
    
    return FALSE;
}

gboolean rsTimecourseParseAlgorithm(const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
    rsTimecourseParameters *p = (rsTimecourseParameters*) data;

    // parse value
    if ( ! strcmp(value, "mean") ) {
        p->algorithm = RSTOOLS_TIMECOURSE_ALGORITHM_MEAN;
    } else if ( ! strcmp(value, "stddev") ) {
        p->algorithm = RSTOOLS_TIMECOURSE_ALGORITHM_STDDEV;
    } else if ( ! strcmp(value, "spca") ) {
        p->algorithm = RSTOOLS_TIMECOURSE_ALGORITHM_SPCA;
    } else if ( ! strcmp(value, "tpca") ) {
        p->algorithm = RSTOOLS_TIMECOURSE_ALGORITHM_TPCA;
    } else if ( ! strcmp(value, "csp") ) {
        p->algorithm = RSTOOLS_TIMECOURSE_ALGORITHM_CSP;
    } else {    
        // any other value should lead to an error
        g_set_error(
            error,
            G_OPTION_ERROR,
            G_OPTION_ERROR_BAD_VALUE,
            "%s: %s",
            option_name,
            "accepted values are 'mean', 'stddev', 'spca', 'tpca' and 'csp'"
        );
        return FALSE;
    }
    
    return TRUE;
}
