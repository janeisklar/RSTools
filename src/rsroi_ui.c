#include "rsroi_ui.h"

rsRoiParameters *rsRoiInitParameters()
{
    rsRoiParameters *p = (rsRoiParameters*)rsMalloc(sizeof(rsRoiParameters));

    p->inputpath                = NULL;
    p->maskpath                 = NULL;
    p->callString               = NULL;
    p->sphereradius             = -1.0;
    p->center                   = rsMakeFloatPoint3D(-9999.9, -9999.9, -9999.9);
    p->cubeDim                  = rsMakeFloatPoint3D(-9999.9, -9999.9, -9999.9);
    p->keepVolume               = FALSE;
    p->nSamples                 = -1;
    p->roiValue                 = 1;
    p->useImageSpaceCoordinates = FALSE;
    p->verbose                  = FALSE;
    p->parametersValid          = FALSE;
    p->context                  = NULL;
    p->input                    = NULL;
    p->mask                     = NULL;
    p->inputEqualsOutput        = FALSE;

    return p;
}

rsRoiParameters *rsRoiParseParams(int argc, char * argv[])
{
    rsRoiParameters *p = rsRoiInitParameters();
    p->callString = rsMergeStringArray(argc, argv);

    // initialize the most common options
    GError *error = NULL;
    p->context = g_option_context_new("\n\nGiven a 4D-Nifti that will be cloned this tool will create a binary mask for the specified region(sphere or cube). Already existing masks may be extended using -keepVolume.");
    GOptionGroup *g = g_option_group_new("common", "Common options", "The most commonly used options", (void*)p, NULL);

    g_option_context_set_summary(p->context, RSTOOLS_VERSION_LABEL);

    GOptionArgFunc cbPoint     = (GOptionArgFunc)rsRoiParsePoint;

    /* long, short, flags, arg, arg_data, desc, arg_desc */
    GOptionEntry entries[] = {
      { "input",           'i', 0, G_OPTION_ARG_FILENAME, &p->inputpath,                "the volume from which the header, dimension(except for temporal)) and alignment infos will be taken. The resulting mask will thus be coregistered to this volume.", "<volume>" },
      { "mask",            'm', 0, G_OPTION_ARG_FILENAME, &p->maskpath,                 "the mask that results from the ROI operations below", "<volume>" },
      { "cube",            'u', 0, G_OPTION_ARG_CALLBACK, cbPoint,                      "use this option if the ROI that is to be added is a cube. The height/weight/depth shoulö be given in mm unless when used with --useImageSpace.", "<w>,<h>,<d>" },
      { "sphere",          's', 0, G_OPTION_ARG_DOUBLE,   &p->sphereradius,             "use this option if the ROI that is to be added is a sphere. Its radius should be given in mm unless when used with --useImageSpace.", "<float>" },
      { "center",          'C', 0, G_OPTION_ARG_CALLBACK, cbPoint,                      "this option specifies the center of the sphere/cube in MNI space(mm) unless when used with --useImageSpace.", "<x>,<y>,<z>" },
      { "roiValue",        'I', 0, G_OPTION_ARG_DOUBLE,   &p->roiValue,                 "intensity value that will be used for the ROI", "<float>" },
      { "useImageSpace",   'V', 0, G_OPTION_ARG_NONE,     &p->useImageSpaceCoordinates, "values supplied with --sphere, --cube and --center will be interpreted as image space coordinates rather than mm", NULL },
      { "keepVolume",      'k', 0, G_OPTION_ARG_NONE,     &p->keepVolume,               "keep values from the input volume. Use this option if you want to add ROIs to a mask that you specified as an input.", NULL },
      { "randomsample",    'n', 0, G_OPTION_ARG_INT,      &p->nSamples,                 "randomly sample <n> voxels from the file specified using --input(through --keepVolume) and/or the other ROI commands(--sphere, --cube)", "<int>" },
      { "verbose",         'v', 0, G_OPTION_ARG_NONE,     &p->verbose,                  "show debug information", NULL},
      { NULL }
    };

    g_option_context_set_main_group(p->context, g);
    g_option_group_add_entries(g, entries);

    // check if parameters are valid
    if ( ! g_option_context_parse(p->context, &argc, &argv, &error) ) {
        fprintf(stderr, "option parsing failed: %s\n", error->message);
        return p;
    }

    p->parametersValid = TRUE;
    return p;
}


void rsRoiFreeParams(rsRoiParameters *p)
{
    rsFree(p->inputpath);
    rsFree(p->maskpath);
    rsFree(p->callString);
    rsFree(p->center);
    rsFree(p->cubeDim);
    rsFree(p->input);
    rsFree(p->mask);
    g_option_context_free(p->context);
    rsFree(p);
}

void rsRoiPrintHelp(rsRoiParameters *p)
{
    fprintf(stdout, "%s\n", g_option_context_get_help(p->context, TRUE, NULL));
}

gboolean rsRoiParsePoint(const gchar *option_name, const gchar *value, gpointer data, GError **error)
{
    rsRoiParameters *p = (rsRoiParameters*) data;

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

    float x=-9999, y=-9999, z=-9999;

    // if we were given exactly 3 numbers separated by comma parse them
    if ( strtok(NULL,",") == NULL && strX != NULL && strY != NULL && strZ != NULL ) {
        x = atof(strX);
        y = atof(strY);
        z = atof(strZ);
    }

    // return success if we sucessfully received 3 numbers
    if ( x > -9999 && y > -9999 && z > -9999 ) {
        FloatPoint3D *point = rsMakeFloatPoint3D(x,y,z);
        if ( ! strcmp(option_name, "--center") ) {
            p->center = point;
        } else {
            p->cubeDim = point;
        }
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
