% Source localization methods using fieldtrip

function performSourceLocalization(data,displayResults)

model = prepareSourceModel(capType);

% Load and preprocess EEG data
load('204NC.mat');
cfg = [];
cleandata = ft_preprocessing(cfg, data);

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


%% Step 4 - Work on the EEG data       
oNum=0; refType = 'unipolar'; % 'unipolar', 'avg','bipolar','csd'
[data,layout] = getData(oNum, refType);
data.elec = elec;

%%%%%%%%% Data segmentation into baseline or stimulus epoch %%%%%%%%%%%%%%%

cfg        = [];
cfg.toilim = [-0.5 0];  % Baseline time period -0.5 to 0
data_pre   = ft_redefinetrial(cfg, data);

cfg.toilim = [0.25 0.75]; % Stimulus time period 0.25 to 0.75
data_post  = ft_redefinetrial(cfg, data);

% compute psd (power spectral density) and csd (cross spectral density) for DICS    

cfg                 = [];
cfg.output          = 'powandcsd';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.tapsmofrq       = 4;  % Smoothing over W = +- 4 Hz. For T=0.5, W=4, TW=2, 3 tapers are needed 
cfg.foi             = 28; % Perform analysis at this frequency
cfg.keeptrials      = 'yes';

freq_pre           = ft_freqanalysis(cfg, data_pre);
freq_post          = ft_freqanalysis(cfg, data_post);

% Source reconstruction using DICS

cfg              = [];
cfg.method       = 'dics';
cfg.grid         = leadfield_bem;
cfg.headmodel    = headmodel_bem;
cfg.frequency    = freq_post.freq;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '2%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.realfilter   = 'yes';

freq_source_post = ft_sourceanalysis(cfg, freq_post);
cfg.sourcemodel.filter = freq_source_post.avg.filter;
freq_source_pre = ft_sourceanalysis(cfg, freq_pre);

data_sourceNorm = freq_source_post;
data_sourceNorm.avg.pow = freq_source_post.avg.pow ./ freq_source_pre.avg.pow;

% Visualization

% step 1: Interpolate
cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'pow';
data_source_interp  = ft_sourceinterpolate(cfg, data_sourceNorm, mri_resliced_orig);

% step 2: Plot
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolormap    = 'jet';
ft_sourceplot(cfg, data_source_interp);