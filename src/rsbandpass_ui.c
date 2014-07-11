#include "rsbandpass_ui.h"

rsBandpassParameters *rsBandpassInitParameters() {
    rsBandpassParameters *p = (rsBandpassParameters*)rsMalloc(sizeof(rsBandpassParameters));
    
    p->inputpath            = NULL;
    p->maskpath             = NULL;
    p->savemaskpath         = NULL;
    p->saveFilteredPath     = NULL;
    p->saveAttenuationPath  = NULL;
    p->callString           = NULL;
    p->fftParams            = NULL;
    p->paddedT              = 0;
    p->freqLow              = -1.0;
    p->freqHigh             = -1.0;
    p->TR                   = -1.0;
    p->verbose              = FALSE;
    p->input                = NULL;
    p->filteredOutput       = NULL;
    p->parametersValid      = FALSE;
    p->mask                 = NULL;
    p->threads              = 1;
    p->rolloff_method       = RSFFTFILTER_CUTOFF;
    p->rolloff              = 10.0;
    p->keepMean             = FALSE;
    p->interface            = NULL;
    p->progressCallback     = NULL;
    
    return p;
}

void rsBandpassFreeParams(rsBandpassParameters *p) {
    free(p->inputpath);
    free(p->maskpath);
    free(p->savemaskpath);
    free(p->saveFilteredPath);
    free(p->saveAttenuationPath);
    free(p->input);
    free(p->filteredOutput);
    free(p->callString);
    free(p->fftParams);
    rsUIDestroyInterface(p->interface);
    free(p);
}

rsBandpassParameters *rsBandpassParseParams(int argc, char * argv[]) {

    rsBandpassParameters *p = rsBandpassInitParameters();
    p->callString = rsMergeStringArray(argc, argv);
    
    // initialize the most common options
    rsUIOption *o;
    p->interface = rsUINewInterface();
    p->interface->description   = "Given a 4D-Nifti and a frequency band this tool will apply FFT-based temporal filtering in the specified frequency range onto the dataset.";
    
    o = rsUINewOption();
    o->name                = "input";
    o->shorthand           = 'i';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->inputpath;
    o->cli_description     = "the input volume that is going to be filtered";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "filtered";
    o->shorthand           = 'f';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->saveFilteredPath;
    o->cli_description     = "the output volume in which the filtered data will be saved";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "f1";
    o->shorthand           = 'l';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->freqLow;
    o->cli_description     = "the lower frequency of the bandpass filter";
    o->cli_arg_description = "<frequency in Hz>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "f2";
    o->shorthand           = 'u';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->freqHigh;
    o->cli_description     = "the upper frequency of the bandpass filter";
    o->cli_arg_description = "<frequency in Hz>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "TR";
    o->shorthand           = 'r';
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->TR;
    o->cli_description     = "the time to repeat (1/sampling frequency) that is used for the FFT";
    o->cli_arg_description = "<rate in s>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "mask";
    o->shorthand           = 'm';
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->maskpath;
    o->cli_description     = "a mask specifying the region that the filter is applied on (may be specified for improved performance)";
    o->cli_arg_description = "<volume>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "keepMean";
    o->shorthand           = 'k';
    o->storage             = &p->keepMean;
    o->cli_description     = "retains the first bin of the FFT (the mean) independent of it being included in the selected frequency range";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "threads";
    o->shorthand           = 't';
    o->type                = G_OPTION_ARG_INT;
    o->storage             = &p->threads;
    o->cli_description     = "number of threads used for processing";
    o->cli_arg_description = "<n>";
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "verbose";
    o->shorthand           = 'v';
    o->storage             = &p->verbose;
    o->cli_description     = "show debug information";
    rsUIAddOption(p->interface, o);
        
    // initialize the more advanced and rather unusual options
    o = rsUINewOption();
    o->name                = "savemask";
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->savemaskpath;
    o->cli_description     = "optional path where the rescaled mask that was specified with the mask option will be saved. The saved file with have the same dimensions as the input volume.";
    o->cli_arg_description = "<volume>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);

#if RS_FFTW_ENABLED == 1
    o = rsUINewOption();
    o->name                = "fftw";
    o->storage             = &p->fftw;
    o->cli_description     = "show debug information";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
#endif  
    
    o = rsUINewOption();
    o->name                = "saveattenuation";
    o->type                = G_OPTION_ARG_FILENAME;
    o->storage             = &p->saveAttenuationPath;
    o->cli_description     = "save txt file that contains the bin's frequencies and the corresponding attenuation weight that was used.";
    o->cli_arg_description = "<txt-file>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
    
    o = rsUINewOption();
    o->name                = "sigmoidrolloff";
    o->type                = G_OPTION_ARG_DOUBLE;
    o->storage             = &p->rolloff;
    o->cli_description     = "uses a sigmoid function for rolling off the passband. The specified number controls how fast it is rolled off with higher numbers corresponding to a quicker rolloff. A good starting point would be 10, then double-check by saving the attenuation file. (does not work in combination with FFTW3)";
    o->cli_arg_description = "<double>";
    o->group               = RS_UI_GROUP_EXTENDED;
    rsUIAddOption(p->interface, o);
    
    // parse
    p->parametersValid = rsUIParse(p->interface, argc, argv);
    return p;
}
