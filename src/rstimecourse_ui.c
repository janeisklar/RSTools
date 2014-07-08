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
    p->context              = NULL;
    p->callString           = NULL;
    p->threads              = 1;

    return p;
}

rsTimecourseParameters* rsTimecourseParseParams(int argc, char * argv[])
{
    rsTimecourseParameters *p = rsTimecourseInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    // initialize the most common options
    GError *error = NULL;
    p->context = g_option_context_new("\n\n Given a 4D-Nifti, this tool extracts the time course for a single voxel or aggregates the timecourse in a mask using different algorithms");
    GOptionGroup *g = g_option_group_new("common", "Common options", "The most commonly used options", (void*)p, NULL);
    
    g_option_context_set_summary(p->context, RSTOOLS_VERSION_LABEL);
    
    GOptionArgFunc cbAlgorithm = (GOptionArgFunc)rsTimecourseParseAlgorithm;
    GOptionArgFunc cbPoint     = (GOptionArgFunc)rsTimecourseParsePoint;
    
     /* long, short, flags, arg, arg_data, desc, arg_desc */
    GOptionEntry entries[] = {
      { "input",             'i', 0, G_OPTION_ARG_FILENAME, &p->inputpath,         "the volume to be regressed", "<volume>" },
      { "output",            'o', 0, G_OPTION_ARG_FILENAME, &p->outputpath,        "file to which the resulting timecourse will be writting to. if omitted the result will be printed directly to stdout", "<*.txt>" },
      { "algorithm",         'a', 0, G_OPTION_ARG_CALLBACK, cbAlgorithm,           "<algo> the algorithm used to aggregate the data within the specified mask(s), e.g. mean, stddev, spca, tpca or csp\n\n\t'mean'\t\t- for every volume in the 4D input file the\n\t\t\t  intensity values in the mask region are \n\t\t\t  meaned, resulting in a meaned timecourse\n\t\t\t  for that region\n\t'stddev'\t- like 'mean', but instead of meaning the\n\t\t\t  mask region, the standard deviation is\n\t\t\t  computed\n\t'tpca'\t\t- performs a PCA on the temporal dimension\n\t\t\t  of the 4D input in the mask region\n\t'spca'\t\t- performs a PCA on the spatial dimension of\n\t\t\t  the 4D input limited to voxels in the mask\n\t\t\t  region\n\t'csp'\t\t- performs the common spatial pattern\n\t\t\t  procedure on the 4D input in the mask\n\t\t\t  region for both masks. The first half of\n\t\t\t  the returned components will maximize the\n\t\t\t  variance for mask1 while the variance for\n\t\t\t  mask2 is minimal. The second half will be\n\t\t\t  exactly the opposite (var(mask1) minimal,\n\t\t\t  var(mask2) maximal)\n", "<algo>" },
      { "mask",              'm', 0, G_OPTION_ARG_FILENAME, &p->maskpath,          "a mask specifying the region that the algorithm is perforned on", "<volume>" },
      { "mask2",             'n', 0, G_OPTION_ARG_FILENAME, &p->mask2path,         "(use only with csp) a second mask specifying the region for the second condition", "<volume>" },
      { "point",             'p', 0, G_OPTION_ARG_CALLBACK, cbPoint,               "speficies a voxel using nifti coordinates(0-based) from which the timecourse is to be extracted", "<x> <y> <z>" },
      { "retainVariance",      0, 0, G_OPTION_ARG_DOUBLE,   &p->minVariance,       "(use only with pca) percentage of explained variance that will be retained, e.g. '0.4'. keep in mind that a higher percentage will result in more components that are to be returned.", "<float>" },
      { "retainComponents",  'c', 0, G_OPTION_ARG_INT,      &p->nComponents,       "(use only with pca or csp) number of PCA/CSP components that will be outputted. This should be a multiple of two when running CSP as the components will be distributed equally over both masks.", "<int>" },
      { "useStandardScores", 's', 0, G_OPTION_ARG_NONE,     &p->useStandardScores, "(use only with pca) remove mean and set std. dev to 1 prior to running pca", NULL },
      { "spatialMap",          0, 0, G_OPTION_ARG_FILENAME, &p->spatialmappath,    "(use only with pca) store the spatial map that is created using the PCA components(4D nifti with one volume for every component)", "<volume>" },
      { "eigenvalues",       'e', 0, G_OPTION_ARG_FILENAME, &p->eigenvaluespath,   "(use only with pca) write out all eigenvalues to the file that is specified with this option", "<txt>" },
      { "threads",           't', 0, G_OPTION_ARG_INT,      &p->threads,           "number of threads used for processing", "<n>" },
      { "verbose",           'v', 0, G_OPTION_ARG_NONE,     &p->verbose,           "show debug information", NULL},
      { NULL }
    };
    
    g_option_context_set_main_group(p->context, g);
    g_option_group_add_entries(g, entries);

    // initialize the more advanced and rather unusual options
    g = g_option_group_new("extended", "Extended options", "Additional options that are rarely going to be used", (void*)p, NULL);
    g_option_context_add_group(p->context, g);

    GOptionEntry extended_entries[] = {
      { "savemask",          0, 0, G_OPTION_ARG_FILENAME, &p->savemaskpath,        "optional path where the rescaled mask specified with -mask will be saved. The saved file with have the same dimensions as the input volume.", "volume" },
      { NULL }
    };

    g_option_group_add_entries(g, extended_entries);

    // check if parameters are valid
    if ( ! g_option_context_parse(p->context, &argc, &argv, &error) ) {
        fprintf(stderr, "option parsing failed: %s\n", error->message);
        return p;
    }

    p->parametersValid = TRUE;
    return p;
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
    g_option_context_free(p->context);
    rsFree(p);
}

void rsTimecoursePrintHelp(rsTimecourseParameters* p)
{
    fprintf(stdout, "%s\n", g_option_context_get_help(p->context, TRUE, NULL));
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
