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
    p->input                    = NULL;
    p->mask                     = NULL;
    p->inputEqualsOutput        = FALSE;
    p->interface                = NULL;

    return p;
}

rsRoiParameters *rsRoiParseParams(int argc, char * argv[])
{
    rsRoiParameters *p = rsRoiInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    rsRoiBuildInterface(p);
    
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

    if ( p->maskpath == NULL ) {
        fprintf(stderr, "A  path for the resulting mask must be specified!\n");
        return p;
    }

    if ( p->center->x < -9999.0 && p->nSamples < 0 ) {
        fprintf(stderr, "ROI center needs to be specified!(-center)!\n");
        return p;
    }

    if ( p->sphereradius <= 0 && p->cubeDim->x < 0 && p->nSamples < 0 ) {
        fprintf(stderr, "ROI sphere radius or cube dimensions needs to be specified!(--sphere, --cube)!\n");
        return p;
    }
    
    p->parametersValid = parsingSuccessful;
    return p;
}

void rsRoiBuildInterface(rsRoiParameters *p)
{
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description       = "Given a 4D-Nifti that will be cloned this tool will create a binary mask for the specified region(sphere or cube). Already existing masks may be extended using -keepVolume.";
    p->interface->gui_description   = "Given a 4D-Nifti that will be cloned this tool will create a binary mask for the specified region(sphere or cube). Already existing masks may be extended using the 'keepVolume' option.";
    p->interface->helpIndent        = 31;
        
    GOptionArgFunc cbPoint = (GOptionArgFunc)rsRoiParsePoint;

    // initialize the most common options
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the volume from which the header, dimension(except for temporal)) and alignment infos will be taken. The resulting mask will thus be coregistered to this volume.";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "mask";
    o->shorthand           = 'm';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->maskpath;
    o->cli_description     = "the mask that results from the ROI operations below";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "cube";
    o->shorthand           = 'c';
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbPoint;
    o->cli_description     = "use this option if the ROI that is to be added is a cube. The height/weight/depth should be given in mm unless when used with --useImageSpace.";
    o->gui_description     = "Use this option if the ROI that is to be added is a cube. The height/weight/depth should be given in mm unless when used with the 'useImageSpace' option.";
    o->cli_arg_description = "<w>,<h>,<d>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "sphere";
    o->shorthand           = 's';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->sphereradius;
    o->cli_description     = "use this option if the ROI that is to be added is a sphere. Its radius should be given in mm unless when used with --useImageSpace.";
    o->gui_description     = "Use this option if the ROI that is to be added is a sphere. Its radius should be given in mm unless when used with the 'useImageSpace' option.";
    o->cli_arg_description = "<float>";
    rsUIAddOption(p->interface, o);
        
    o = rsUINewOption();
    o->name                = "center";
    o->shorthand           = 'C';
    o->type                = G_OPTION_ARG_CALLBACK;
    o->storage             = cbPoint;
    o->cli_description     = "this option specifies the center of the sphere/cube in MNI space(mm) unless when used with --useImageSpace.";
    o->gui_description     = "this option specifies the center of the sphere/cube in MNI space(mm) unless when used with the 'useImageSpace' option.";
    o->cli_arg_description = "<w>,<h>,<d>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "roiValue";
    o->shorthand           = 'V';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->roiValue;
    o->cli_description     = "intensity value that will be used for the ROI";
    o->cli_arg_description = "<float>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "useImageSpace";
    o->shorthand           = 'I';
    o->storage             = &p->useImageSpaceCoordinates;
    o->cli_description     = "values supplied with --sphere, --cube and --center will be interpreted as image space coordinates rather than mm";
    o->gui_description     = "values supplied with the 'sphere', 'cube' and 'center' option will be interpreted as image space coordinates rather than mm";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "keepVolume";
    o->shorthand           = 'k';
    o->storage             = &p->keepVolume;
    o->cli_description     = "keep values from the input volume. Use this option if you want to add ROIs to a mask that you specified as an input.";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "randomsample";
    o->shorthand           = 'S';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->nSamples;
    o->cli_description     = "randomly sample <n> voxels from the file specified using --input(through --keepVolume) and/or the other ROI commands(--sphere, --cube)";
    o->gui_description     = "randomly sample <n> voxels from the file specified using the 'input' combined with the 'keepVolume' option and/or the other ROI options: 'sphere', cube'";
    o->cli_arg_description = "<int>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    rsUIAddOption(p->interface, o);
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
    rsUIDestroyInterface(p->interface);
    rsFree(p);
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
