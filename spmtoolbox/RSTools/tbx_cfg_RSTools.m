function rstools = tbx_cfg_RSTools
% Configuration file for toolbox 'RSTools'
%_______________________________________________________________________

% André Hoffmann

% ---------------------------------------------------------------------
% data images to extract the ROI time-course from
% ---------------------------------------------------------------------
data           = cfg_files;
data.tag       = 'data';
data.name      = 'Images to extract the ROI time-course from';
data.filter    = 'image';
data.ufilter   = '.*';
data.num       = [1 Inf];
data.help      = {['Select the images from which the ROI will be extracted ']};

% ---------------------------------------------------------------------
% alg_sum sum algorithm
% ---------------------------------------------------------------------
alg_sum        = cfg_menu;
alg_sum.tag    = 'SumAlgorithm';
alg_sum.name   = 'Summation Algorithm';
alg_sum.help   = {'The algorithm that will be used to summarize the data in the ROI. For further information see {spmPath}/toolbox/marsbar/@maroi/get_marsy.m.'};
alg_sum.labels = {
                    'mean'
                    'median'
                    'eigen1'
                    'wtmean'
}';
alg_sum.values = {
                    'mean'
                    'median'
                    'eigen1'
                    'wtmean'
}';
%alg_sum.def    = @(val)spm_get_defaults('rstools.roiextraction.sumalgorithm', val{:});

% ---------------------------------------------------------------------
% roi ROI binary mask
% ---------------------------------------------------------------------
roi            = cfg_files;
roi.tag        = 'roi';
roi.name       = 'Binary ROI Volume Mask';
roi.filter     = 'image';
roi.ufilter    = '.*';
roi.num        = [1 1];
roi.help       = {['Select a binary roi volume mask.']};

% ---------------------------------------------------------------------
% prefix prefix of the output file
% ---------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames. Default prefix is ''tc_roi''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'tc_roi'};

% ---------------------------------------------------------------------
% ROIExtraction ROI Extraction
% ---------------------------------------------------------------------
ROIExtraction         = cfg_exbranch;
ROIExtraction.tag     = 'ROIExtraction';
ROIExtraction.name    = 'ROI Time-course Extraction';
ROIExtraction.val     = {data alg_sum roi prefix };
ROIExtraction.help    = {'Batch wrapper for MarsBaR'};
ROIExtraction.prog    = @spm_local_roiExtraction;
ROIExtraction.vout    = @vout_roiExtraction;

% ---------------------------------------------------------------------
% regressor Regressor specification
% ---------------------------------------------------------------------
regressor            = cfg_files;
regressor.tag        = 'regressor';
regressor.name       = 'Regressor(s)';
regressor.filter     = {'mat','timecourse'};
regressor.ufilter    = '.*';
regressor.num        = [1 1];
regressor.help       = {['Select a regressor that will be merged with the others.']};

% ---------------------------------------------------------------------
% regressors List of several regressors
% ---------------------------------------------------------------------
regressors           = cfg_repeat;
regressors.tag       = 'regressors';
regressors.name      = 'Multiple Regressors';
regressors.help      = {'Add as many regressors as should be merged together.'};
regressors.values    = {regressor };
regressors.num       = [1 Inf];

% ---------------------------------------------------------------------
% filename filename of the output file
% ---------------------------------------------------------------------
filename         = cfg_entry;
filename.tag     = 'filename';
filename.name    = 'Filename';
filename.help    = {'Specify the name of file to which the merged regressors will be written to.'};
filename.strtype = 's';
filename.num     = [1 Inf];
filename.val     = {'merged_regressors.txt'};

% ---------------------------------------------------------------------
% RegressorMerging Merging of Regressors
% ---------------------------------------------------------------------
RegressorMerging         = cfg_exbranch;
RegressorMerging.tag     = 'RegressorMerging';
RegressorMerging.name    = 'Merge Multiple Regressors';
RegressorMerging.val     = {regressors filename };
RegressorMerging.help    = {'Merges multiple regressors like realignment parameters or timecourses together into one, so that it can be used as a dependency in the fMRI model specification routine.'};
RegressorMerging.prog    = @spm_local_regressorMerging;
RegressorMerging.vout    = @vout_regressorMerging;

% ---------------------------------------------------------------------
% dataLR data used for the linear regression
% ---------------------------------------------------------------------
dataLR           = cfg_files;
dataLR.tag       = 'dataLR';
dataLR.name      = 'Data used for the linear regression';
dataLR.filter    = 'image';
dataLR.ufilter   = '.*';
dataLR.num       = [1 Inf];
dataLR.help      = {['Select the images that will be used for the linear regression']};

% ---------------------------------------------------------------------
% dataLR data used for the linear regression
% ---------------------------------------------------------------------
maskLR           = cfg_files;
maskLR.tag       = 'maskLR';
maskLR.name      = 'Explicit mask that is applied to the data before the regression';
maskLR.filter    = 'image';
maskLR.ufilter   = '.*';
maskLR.num       = [0 1];
maskLR.help      = {['Select a binary mask that will be applied to the data before regression such as tpm/brainmask.nii.']};

% ---------------------------------------------------------------------
% prefix prefix of the output file
% ---------------------------------------------------------------------
prefixLR         = cfg_entry;
prefixLR.tag     = 'prefixLR';
prefixLR.name    = 'Filename Prefix';
prefixLR.help    = {'Specify the string to be prepended to the filenames. Default prefix is ''ResI_''.'};
prefixLR.strtype = 's';
prefixLR.num     = [1 Inf];
prefixLR.val     = {'ResI_'};

% ---------------------------------------------------------------------
% LinearRegression Simple linear regression
% ---------------------------------------------------------------------
LinearRegression         = cfg_exbranch;
LinearRegression.tag     = 'LinearResgression';
LinearRegression.name    = 'Simple linear regression';
LinearRegression.val     = {regressor dataLR prefixLR maskLR};
LinearRegression.help    = {'This is mainly a wrapper for Matlab''s regress() function. As opposed to SPM''s model specification, no high pass-filter or any other additional corrections are being applied.'};
LinearRegression.prog    = @spm_local_linearRegression;
LinearRegression.vout    = @vout_linearRegression;

% ---------------------------------------------------------------------
% dataDF data used for the linear regression
% ---------------------------------------------------------------------
dataDF           = cfg_files;
dataDF.tag       = 'dataDF';
dataDF.name      = 'Data that is to be filtered';
dataDF.filter    = 'image';
dataDF.ufilter   = '.*';
dataDF.num       = [1 Inf];
dataDF.help      = {['Select the volumes on which the filter will be applied on']};

% ---------------------------------------------------------------------
% prefixDF prefix of the output file
% ---------------------------------------------------------------------
prefixDF         = cfg_entry;
prefixDF.tag     = 'prefixDF';
prefixDF.name    = 'Filename Prefix';
prefixDF.help    = {'Specify the string to be prepended to the filenames. Default prefix is ''Filtered_''.'};
prefixDF.strtype = 's';
prefixDF.num     = [1 Inf];
prefixDF.val     = {'Filtered_'};

% ---------------------------------------------------------------------
% filterSpecification matlab code that creates the filter
% ---------------------------------------------------------------------
filterSpecification         = cfg_entry;
filterSpecification.tag     = 'filterSpecification';
filterSpecification.name    = 'Filter specification';
filterSpecification.help    = {'Specify the filter that will be applied on the data by creating a Matlab filter object. You can use fdatool to generate the corresponding code for your filter. Then transform it in the following format before copying it into this field: "design(fdesign.lowpass(0.08, 0.16, 1, 80, 1/1.8), ''butter'', ''MatchExactly'', ''stopband'')"'};
filterSpecification.strtype = 's';
filterSpecification.num     = [1 Inf];
filterSpecification.val     = {'design(fdesign.lowpass(0.08, 0.16, 1, 80, 1/1.8), ''butter'', ''MatchExactly'', ''stopband'')'};

% ---------------------------------------------------------------------
% RT TR time
% ---------------------------------------------------------------------
RT         = cfg_entry;
RT.tag     = 'RT';
RT.name    = 'TR in s';
RT.help    = {'Please specify the TR time, which will be needed to correct for the delay of the filter.'};
RT.strtype = 'r';
RT.num     = [1 1];
RT.val     = {1.8};

% ---------------------------------------------------------------------
%  DigitalFilter Apply digital filter
% ---------------------------------------------------------------------
DigitalFilter         = cfg_exbranch;
DigitalFilter.tag     = 'DigitalFilter';
DigitalFilter.name    = 'Digital filter';
DigitalFilter.val     = {dataDF prefixDF filterSpecification RT};
DigitalFilter.help    = {'This batch task will apply a filter that can be specified on every voxel in the given dataset. Please keep in mind that the filter will have a delay which will result in less output volumes. The data will be cut off at the beginning of the recording. This delay depends on the specified filter and will be saved in a ''-delay.txt''-file with the same name as the output volume which contains the delay in seconds and a second ''-cutoff.txt'' that holds the number of dropped volumes).'};
DigitalFilter.prog    = @spm_local_digitalFilter;
DigitalFilter.vout    = @vout_digitalFilter;
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% RT TR time
% ---------------------------------------------------------------------
RT         = cfg_entry;
RT.tag     = 'RT';
RT.name    = 'TR in s';
RT.help    = {'The TR time.'};
RT.strtype = 'r';
RT.num     = [1 1];
RT.val     = {1.8};

% ---------------------------------------------------------------------
% f Frequency range
% ---------------------------------------------------------------------
frequencyRange         = cfg_entry;
frequencyRange.tag     = 'frequencyRange';
frequencyRange.name    = 'Frequency range in Hz';
frequencyRange.help    = {'Provide the lower and upper bounds of the frequency range that should be filtered. For a low-pass filter use [0 f], for a high-pass filter [f Inf] and for a band-pass filter [f1 f2].'};
frequencyRange.strtype = 'r';
frequencyRange.num     = [1 2];
frequencyRange.val     = {[0.009 0.08]};

% ---------------------------------------------------------------------
%  FFTFilter Apply filter using FFT
% ---------------------------------------------------------------------
FFTFilter         = cfg_exbranch;
FFTFilter.tag     = 'FFTFilter';
FFTFilter.name    = 'FFT filter';
FFTFilter.val     = {dataDF prefixDF frequencyRange RT};
FFTFilter.help    = {'This batch task will apply a frequency filter by transforming the data into frequency space using FFT, cutting out the provided frequency range and transforming it back using IFFT.'};
FFTFilter.prog    = @spm_local_FFTFilter;
FFTFilter.vout    = @vout_FFTFilter;
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% description a description of the command
% ---------------------------------------------------------------------
description         = cfg_entry;
description.tag     = 'description';
description.name    = 'Description of this task';
description.help    = {'This field doesn''t serve any functional purpose and exists soley for the purpose of documentation(text will be printed during execution).'};
description.strtype = 's';
description.num     = [1 Inf];
description.val     = {'>>Description of the task<<'};

% ---------------------------------------------------------------------
% command unix command to be executed
% ---------------------------------------------------------------------
command         = cfg_entry;
command.tag     = 'command';
command.name    = 'Command';
command.help    = {'The unix command to be executed'};
command.strtype = 's';
command.num     = [1 Inf];
command.val     = {'echo foo;'};

% ---------------------------------------------------------------------
% dataUC data tbat will be passed onto the next module
% ---------------------------------------------------------------------
dataUC           = cfg_files;
dataUC.tag       = 'dataUC';
dataUC.name      = 'Output Data';
dataUC.filter    = 'image';
dataUC.ufilter   = '.*';
dataUC.num       = [1 Inf];
dataUC.help      = {['Select the volumes that will be passed onto the next module']};

% ---------------------------------------------------------------------
%  UnixCommand Execute unix command
% ---------------------------------------------------------------------
UnixCommand         = cfg_exbranch;
UnixCommand.tag     = 'UnixCommand';
UnixCommand.name    = 'Execute Unix Command';
UnixCommand.val     = {description command dataUC};
UnixCommand.help    = {'This task will execute a supplied unix command.'};
UnixCommand.prog    = @spm_local_UnixCommand;
UnixCommand.vout    = @vout_UnixCommand;
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% rstools RSTools
% ---------------------------------------------------------------------
rstools         = cfg_choice;
rstools.tag     = 'RSTools';
rstools.name    = 'RS Tools';
rstools.help    = {'This is a toolbox that provides a couple of helpers for the preprocessing of RS data.'};
rstools.values  = {ROIExtraction RegressorMerging LinearRegression DigitalFilter FFTFilter UnixCommand};

%======================================================================
function out = spm_local_roiExtraction(varargin)
if ~isdeployed 
    addpath(fullfile(spm('dir'),'toolbox','RSTools'));
    addpath(fullfile(spm('dir'),'toolbox','marsbar'));
    marsbar('on');
end
    
job = varargin{1};

% convert cell to string array to get it to work with the marsbar api
max_str_len = 0;

for i=1:length(job.data)
    max_str_len = max(length(job.data{i}), max_str_len);
end

data = [];

for i=1:length(job.data)
   str  = [job.data{i} repmat(' ',1,max_str_len-length(job.data{i}))];
   data = [data; str];
end

rois = maroi_image(job.roi{1});
mY   = get_marsy(rois, data, job.SumAlgorithm);  % extract data into marsy data object
y    = summary_data(mY);

% put together new file name
basename = fliplr(regexprep(fliplr(job.data{1}), '^[^,]*,?[^\.]*\.([^/]*)/(.*)$', ['$1' fliplr(job.prefix) '/$2']));
filename = [basename '.txt'];
out.timecourse = {filename};

% write to file
dlmwrite(filename, y, 'delimiter', '');

% write roi
filename = [basename '.roi.img'];
save_as_image(rois, filename);

%======================================================================
function dep = vout_roiExtraction(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Extracted ROI Time-course';
dep(1).src_output = substruct('.','timecourse');
dep(1).tgt_spec   = cfg_findspec({{'filter','timecourse','strtype','e'}});

%======================================================================
function out = spm_local_regressorMerging(varargin)
job = varargin{1};
%regressors = job.regressor;

regressors = [];

for i = 1:length(job.regressor)
    regressor = job.regressor{i}{1};
    regressors = [regressors importdata(regressor)];
end

save(job.filename ,'regressors','-ascii');

out.multi_reg = {job.filename};

%======================================================================
function dep = vout_regressorMerging(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Multiple regressors';
dep(1).src_output = substruct('.','multi_reg');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat'}});

%======================================================================
function out = spm_local_linearRegression(varargin)
job = varargin{1};

% load job params
prefix = job.prefixLR;
regressorFile = job.regressor{1};
maskFile = job.maskLR{1};

% load volumes
volumes = spm_vol(char(job.dataLR));
volumesData = spm_read_vols(volumes);
nX=size(volumesData,1);
nY=size(volumesData,2);
nZ=size(volumesData,3);
t = size(volumesData, 4);
NaNrep = spm_type(volumes(1).dt(1),'nanrep');

% load mask
if ~isempty(maskFile)
   maskVol = spm_vol(char(maskFile));
   
   % resample mask to fit the data's dimensions
   C = spm_bsplinc(maskVol, [0 0 0 0 0 0]');
   mask = true(nX,nY,nZ);
   [x1,x2] = ndgrid(1:nX,1:nY);
   for x3 = 1:nZ
       M  = inv(volumes(1).mat\maskVol.mat);
       y1 = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
       y2 = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
       y3 = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
       mask(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
   end
end

% load regressors
regressors = importdata(regressorFile);
regressors = [regressors, ones(t,1)]; % add constant factor(required by regress())

% prepare residuals 
residuals = ones(size(volumesData));

% do linear regression on a voxel-basis
spm_progress_bar('Init',nX*nY*nZ,'Linear regression','Voxels complete');
c = 0;
for x = 1:nX
    for y = 1:nY
        for z = 1:nZ            
            if mask(x,y,z)==0
                r = NaNrep;
            else 
                voxel = volumesData(x,y,z,:);
                [b,bint,r] = regress(voxel(:),regressors);
            end
            residuals(x,y,z,:) = r;
        end
    end
    c = c + nY*nZ;
    spm_progress_bar('Set',c);
end
spm_progress_bar('Clear');

% create residual filename
basename = fliplr(regexprep(fliplr(volumes(1).fname), '^[^,]*,?[^\.]*\.([^/]*)/(.*)$', ['$1' fliplr(prefix) '/$2']));
filename = [basename '.nii'];

spm_unlink(filename);

% create residuals file
residualVolumes = volumes;
%residualDimensions = [min(residualVolumes), max(residualVolumes)]
%for k=1:size(residualVolumes, 1)
%    residualVolumes(k).fname  = spm_file(volumes(k).fname, 'prefix', prefix);
%    if isfield(residualVolumes(k),'descrip')
%        desc = [residualVolumes(k).descrip ' '];
%    else
%        desc = '';
%    end
%    residualVolumes(k).descrip = [desc 'residuals from linear regression'];
%%    residualVolumes(k).pinfo = [Inf Inf 0]';
%end
%residualVolumes = rmfield(residualVolumes, 'pinfo');
%residualVolumes = spm_create_vol(residualVolumes);

for k=1:size(residualVolumes, 1)
    vol = struct('fname', spm_file(volumes(k).fname, 'prefix', prefix), ...
         'dim',     [nX nY nZ], ...
         'dt',      [spm_type('float32') spm_platform('bigend')], ...
         'mat',     volumes(k).mat, ...
         'pinfo',   [1 0 0]', ...
         'descrip', 'residuals from linear regression', ...
         'n',       volumes(k).n ...
    );
    
    vol = spm_create_vol(vol);
    spm_write_vol(vol, double(residuals(:,:,:,k)));
         
    if k<2
        residualVolumes = vol;
    else
        residualVolumes(k) = vol;
    end
end

out.residuals = cell(t,1);
for k=1:t
    out.residuals{k} = [residualVolumes(k).fname ',' residualVolumes(k).n(1)];
end

% create mask file
maskVolume = volumes(1);
maskVolume.fname  = spm_file(maskVolume.fname, 'prefix', ['mask_' prefix]);
spm_write_vol(maskVolume, mask);

%======================================================================
function dep = vout_linearRegression(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Residuals of the linear regression';
dep(1).src_output = substruct('.','residuals');
dep(1).tgt_spec   = cfg_findspec({{'filter','image'}});

%======================================================================
function out = spm_local_digitalFilter(varargin)
job = varargin{1};

% load job params
prefix = job.prefixDF;
data = job.dataDF;
filterHandle = eval(job.filterSpecification);
RT = job.RT(1);

% load volumes
volumes = spm_vol(char(data));
volumesData = spm_read_vols(volumes);
nX=size(volumesData,1);
nY=size(volumesData,2);
nZ=size(volumesData,3);
t = size(volumesData, 4);
NaNrep = spm_type(volumes(1).dt(1),'nanrep');

% determine the coefficients of the filter
if isprop(filterHandle, 'sosMatrix')
    % IIR filter
    SOS = filterHandle.sosMatrix;
    G = filterHandle.ScaleValues;
else
    % FIR filter
    SOS = filterHandle.Numerator;
    G = 1;
end

% apply the filter on a voxel-basis
spm_progress_bar('Init',nX*nY*nZ,'Digital filtering','Voxels complete');
c = 0;
for x = 1:nX
    for y = 1:nY
        if any(volumesData(x,y,:,:) ~= NaNrep)
            line = reshape(volumesData(x,y,:,:), [nZ, t]);
            f = filtfilt(SOS, G, line');
            volumesData(x,y,:,:) = f';
        end
    end
    c = c + nY*nZ;
    spm_progress_bar('Set',c);
end
spm_progress_bar('Clear');

% create the filename for the filtered volumes
basename = fliplr(regexprep(fliplr(volumes(1).fname), '^[^,]*,?[^\.]*\.([^/]*)/(.*)$', ['$1' fliplr(prefix) '/$2']));
filename = [basename '.nii'];

spm_unlink(filename);

% create filtered volumes
filteredVolumes = volumes(1:t);

for k=1:t
    vol = struct('fname', spm_file(volumes(k).fname, 'prefix', prefix), ...
         'dim',     [nX nY nZ], ...
         'dt',      [spm_type('float32') spm_platform('bigend')], ...
         'mat',     volumes(k).mat, ...
         'pinfo',   [1 0 0]', ...
         'descrip', ['volumes filtered with the following filter: ' job.filterSpecification], ...
         'n',       volumes(k).n ...
    );
    
    vol = spm_create_vol(vol);
    spm_write_vol(vol, double(volumesData(:,:,:,k)));
         
    if k<2
        filteredVolumes = vol;
    else
        filteredVolumes(k) = vol;
    end
end

out.data = cell(t,1);
for k=1:t
    out.data{k} = [filteredVolumes(k).fname ',' filteredVolumes(k).n(1)];
end

%======================================================================
function dep = vout_digitalFilter(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Filtered data';
dep(1).src_output = substruct('.','data');
dep(1).tgt_spec   = cfg_findspec({{'filter','image'}});

%======================================================================
function out = spm_local_FFTFilter(varargin)
job = varargin{1};

% load job params
prefix = job.prefixDF
data = job.dataDF
frequencyRange = job.frequencyRange
RT = job.RT(1)

% load volumes
volumes = spm_vol(char(data));
volumesData = spm_read_vols(volumes);
nX=size(volumesData,1);
nY=size(volumesData,2);
nZ=size(volumesData,3);
t = size(volumesData, 4);
NaNrep = spm_type(volumes(1).dt(1),'nanrep');


% compute how many volumes will have to be cut off
cutoff = mean(grpdelay(filterHandle))

% compute the delay of the filter in seconds
delay = cutoff * RT

% apply the filter on a voxel-basis
spm_progress_bar('Init',nX*nY*nZ,'Digital filtering','Voxels complete');
c = 0;
for x = 1:nX
    for y = 1:nY
        for z = 1:nZ            
            if all(volumesData(x,y,z,:) == 0)
                f = NaNrep;
            else 
                voxel = volumesData(x,y,z,:);
                f = filter(filterHandle, voxel(:));
            end
            volumesData(x,y,z,:) = f;
        end
    end
    c = c + nY*nZ;
    spm_progress_bar('Set',c);
end
spm_progress_bar('Clear');

% create the filename for the filtered volumes
basename = fliplr(regexprep(fliplr(volumes(1).fname), '^[^,]*,?[^\.]*\.([^/]*)/(.*)$', ['$1' fliplr(prefix) '/$2']));
filename = [basename '.nii'];

spm_unlink(filename);

% create delay and cutoff files
save([basename '-delay.txt'] ,'delay','-ascii');
save([basename '-cutoff.txt'] ,'cutoff','-ascii');

% create filtered volumes
newT = t-cutoff
filteredVolumes = volumes{1:newT};

for i=1:newT
    k = i+cutoff;
    vol = struct('fname', spm_file(volumes(i).fname, 'prefix', prefix), ...
         'dim',     [nX nY nZ], ...
         'dt',      [spm_type('float32') spm_platform('bigend')], ...
         'mat',     volumes(i).mat, ...
         'pinfo',   [1 0 0]', ...
         'descrip', ['volumes filtered with the following filter: ' job.filterSpecification], ...
         'n',       volumes(i).n ...
    );
    
    vol = spm_create_vol(vol);
    spm_write_vol(vol, double(volumesData(:,:,:,k)));
         
    if i<2
        filteredVolumes = vol;
    else
        filteredVolumes(i) = vol;
    end
end

out.data = cell(newT,1);
for k=1:newT
    out.data{k} = [filteredVolumes(k).fname ',' filteredVolumes(k).n(1)];
end
out.delay = delay;
out.cutoff = cutoff;

%======================================================================
function dep = vout_FFTFilter(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Filtered data';
dep(1).src_output = substruct('.','data');
dep(1).tgt_spec   = cfg_findspec({{'filter','image'}});

%======================================================================
function out = spm_local_UnixCommand(varargin)

job = varargin{1};
command = job.command;
description = job.description;
out.data = job.dataUC;

disp(['Job description: ' description]);

[status, error]  = unix(command);

if status ~= 0
    throw(MException('RSTools:UnixCommand','The execution of the command ''%s'' returned an error.', command));
else
    disp('Execution successful!');
end

%======================================================================
function dep = vout_UnixCommand(job)
dep(1)            = cfg_dep;
dep(1).sname      = sprintf('Result from the execution of ''%s''', job.description);
dep(1).src_output = substruct('.','data');
dep(1).tgt_spec   = cfg_findspec({{'filter','image'}});
