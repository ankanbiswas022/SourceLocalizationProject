% EEG Source Analysis Pipeline
% Dependency: FieldTrip toolbox
% This script performs complete EEG source analysis including:
% - MRI preprocessing
% - Head model creation
% - Source model generation
% - Lead field computation
% - Time-frequency analysis
% - Source reconstruction (using DICS and eLoreta)

folderSave = 'prepDataSourceModel'; % All data is saved here
addpath(fullfile(pwd,folderSave));

%%%%%%%%%%%%%%%%%% Get boundary meshes if not available %%%%%%%%%%%%%%%%%%%
fileName = fullfile(folderSave,'bnd.mat');

if exist(fileName,'file')
    tmp = load(fileName);
    bnd = tmp.bnd;
else

    % Load and preprocess MRI data
    load mri_segmented.mat
    cfg = [];
    mri_resliced_segmented = ft_volumereslice(cfg, mri_segmented);
    seg_i = ft_datatype_segmentation(mri_resliced_segmented, 'segmentationstyle', 'indexed');

    % Configure initial visualization
    cfg = [];
    cfg.funparameter = 'seg';
    cfg.funcolormap = gray(4);     % Distinct color per tissue
    cfg.location = 'center';
    cfg.atlas = seg_i;
    cfg.resolution = 1;

    % Create boundary meshes
    cfg = [];
    cfg.tissue = {'brain', 'skull', 'scalp'};
    cfg.numvertices = [1000 1000 1000];
    bnd = ft_prepare_mesh(cfg, mri_segmented);
    save(fileName,'bnd');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Electrode Locations %%%%%%%%%%%%%%%%%%%%%%
fileName = fullfile(folderSave,'electrode.mat');

if exist(fileName,'file')
    tmp = load(fileName);
    elec = tmp.elec;
else
    % Write the code to generate the variable elec from Montage data
end

figure;
ft_plot_mesh(bnd(1), 'facecolor','none');
ft_plot_mesh(bnd(2), 'facecolor','none');
ft_plot_mesh(bnd(3), 'facecolor','none');
hold on;
ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'on');

% Electrode realignment
mesh = ft_convert_units(bnd, 'cm');
cfg = [];
cfg.method = 'interactive';
cfg.elec = elec;
cfg.headshape = bnd;
elec_realigned = ft_electroderealign(cfg);

% Create BEM head model
cfg = [];
if ispc
    cd('C:\Program Files\OpenMEEG\bin');
    cfg.method = 'openmeeg';
else
    cfg.method = 'dipoli';
end
cfg.tissue = {'brain', 'skull', 'scalp'};
headmodel_bem = ft_prepare_headmodel(cfg, bnd);
save headmodel_bem.mat headmodel_bem
save headmodel_dipoli headmodel_bem

% Visualize head model
figure;
ft_plot_headmodel(headmodel_bem, 'facecolor', 'skin', 'facealpha', 0.7);
hold on;
ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'on');
camlight;

% Create resolution-based source model
cfg = [];
cfg.grid.resolution = 5;
cfg.grid.unit = 'mm';
cfg.threshold = 0.1;
cfg.smooth = 5;
cfg.headmodel = headmodel_bem;
cfg.inwardshift = 1;
sourcemodel = ft_prepare_sourcemodel(cfg);

% Visualize source model
figure;
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:), 'vertexmarker', 'o', 'vertexsize', 5);
save sourcemodel sourcemodel;

% Create fMRI-based source model
cfg = [];
cfg.method = 'basedonmri';
cfg.grid.resolution = 10;
cfg.grid.unit = 'mm';
cfg.threshold = 0.1;
cfg.smooth = 5;
cfg.headmodel = headmodel_bem;
cfg.inwardshift = 1;
cfg.mri = mri_segmented;
sourcemodel = ft_prepare_sourcemodel(cfg);

% Load and preprocess EEG data
load('204NC.mat');
cfg = [];
cleandata = ft_preprocessing(cfg, data);

% Compute lead field
cfg = [];
cfg.sourcemodel = sourcemodel;
cfg.headmodel = headmodel_bem;
cfg.elec = elec;
cfg.channel = cleandata.label;
cfg.reducerank = 3;
cfg.normalize = 'no';
leadfield_bem = ft_prepare_leadfield(cfg);
save leadfield_bem.mat leadfield_bem

%%%%%%%%%%%%%%%%%%%% Estimate source using DICS %%%%%%%%%%%%%%%%%%%%%%%%
% Load data files if not already in workspace
required_files = {'sourcemodel', 'headmodel_bem', 'leadfield_bem', 'mri_resliced_orig', 'mri_segmented'};
for i = 1:length(required_files)
    if ~exist(required_files{i}, 'var')
        load([required_files{i}, '.mat']);
    end
end
% Prepare time windows for analysis
cfg = [];
cfg.toilim = [-0.5 0];
data_pre = ft_redefinetrial(cfg, cleandata);
cfg.toilim = [0.25 0.75];
data_post = ft_redefinetrial(cfg, cleandata);

% Perform spectral analysis
cfg = [];
cfg.output = 'powandcsd';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = 4;
cfg.foi = 28;
cfg.keeptrials = 'yes';
freq_pre = ft_freqanalysis(cfg, data_pre);
freq_post = ft_freqanalysis(cfg, data_post);

% Perform DICS beamforming
cfg = [];
cfg.method = 'dics';
cfg.grid = leadfield_bem;
cfg.headmodel = headmodel_bem;
cfg.frequency = freq_post.freq;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda = '2%';
cfg.dics.keepfilter = 'yes';
cfg.dics.realfilter = 'yes';

% Compute source analysis
freq_source_post = ft_sourceanalysis(cfg, freq_post);
cfg.sourcemodel.filter = freq_source_post.avg.filter;
freq_source_pre = ft_sourceanalysis(cfg, freq_pre);

% Calculate normalized source power
data_sourceNorm = freq_source_post;
data_sourceNorm.avg.pow = freq_source_post.avg.pow ./ freq_source_pre.avg.pow;

% Source Data Visualization
% Interpolate source data onto MRI for visualization
cfg = [];
cfg.downsample = 2;
cfg.parameter = 'pow';
data_source_interp = ft_sourceinterpolate(cfg, data_sourceNorm, mri_resliced_orig);

% Visualize interpolated source data
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'maxabs';
cfg.funcolormap = 'jet';        % Use jet colormap for better contrast
ft_sourceplot(cfg, data_source_interp);

%%%%%%%%%%%%%%%%%%%% Estimate source using eLoreta %%%%%%%%%%%%%%%%%%%%%%%%
% Load data files if not already in workspace
required_files = {'sourcemodel', 'headmodel_bem', 'leadfield_bem', 'mri_resliced_orig', 'mri_segmented'};
for i = 1:length(required_files)
    if ~exist(required_files{i}, 'var')
        load([required_files{i}, '.mat']);
    end
end
%%%%%%%%%%%%%%%%%%%% calculate covariance matrix   %%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-0.5 0]; % baseline period
tlckFC = ft_timelockanalysis(cfg, cleandata);

%%%%%%%%%%%%%%%%%%%%% define the CFG and analysze source %%%%%%%%%%%%%%%%%%
cfg             = [];
cfg.method      = 'eloreta';
cfg.sourcemodel = leadfield_bem;
cfg.headmodel   = headmodel_bem;
sourceTimePre   = ft_sourceanalysis(cfg, tlckFC);  % compute the source model

% step 1: Interpolate
cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'pow';
data_source_interp  = ft_sourceinterpolate(cfg, sourceTimePre, mri_resliced_orig);

% spatially normalize the anatomy and functional data to MNI coordinates
cfg = [];
cfg.nonlinear = 'no';
sourceDiffIntNorm = ft_volumenormalise(cfg, data_source_interp);

% step 2: Plot
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolormap    = 'jet';
ft_sourceplot(cfg, data_source_interp);
% idem, in MNI (or SPM) coordinates
ft_sourceplot(cfg, sourceDiffIntNorm);
