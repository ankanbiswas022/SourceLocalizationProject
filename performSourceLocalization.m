% Source localization methods using fieldtrip

function sourceData = performSourceLocalization(data,method,displayResults)
if ~exist('method','var');              method = 'DICS';                    end
if ~exist('displayResults','var');      displayResults = 0;                 end
if ~exist('data', 'var'),               data = load('204NC.mat').data;      end
if ~exist('capType','var');             capType = 'actiCap64';              end

% Get bad channel names from indices
badChannels = data.label(data.badElecs);

% prepare/load the model structure
model = prepareSourceModel(capType);
leadfield = model.leadfield ;
headmodel= model.headmodel;

% Preprocessing configuration
cfg = [];
cfg.demean = 'yes';         % Remove mean value
cfg.reref = 'yes';          % Re-reference the data
cfg.refmethod = 'avg';      % Average reference
cfg.refchannel = 'all';     % Use all channels for referencing

% preprocess the data for average referencing
data_reref = ft_preprocessing(cfg, data);

% Interpolate the bad electrodes
% Prepare neighbours using template method
% note: requires the neighours file witin fieldtrip/tempate directory
cfg = [];
cfg.method = 'template';    % Use template method instead of spline
cfg.template = 'elec1010';  % Using standard 10-10 template
cfg.layout = 'elec1010';    % Specify the layout file
cfg.feedback = 'no';       % Show feedback about the neighbours
neighbours = ft_prepare_neighbours(cfg);
% interpolate
cfg = [];
cfg.method = 'spline';      % Keep spline for interpolation
cfg.badchannel = badChannels;  % Using the extracted bad channel names
cfg.neighbours = neighbours;    % Use the prepared neighbours
% Perform interpolation
cleandata = ft_channelrepair(cfg, data_reref);

%%%%%%%%%%%%%%%%%%%% Estimate source using DICS %%%%%%%%%%%%%%%%%%%%%%%%
% Load data files if not already in workspace
required_files = {'mri_resliced_orig', 'mri_segmented'};
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

if strcmp(method,'DICS')
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
    cfg.grid = leadfield;
    cfg.headmodel = headmodel;
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

    if displayResults
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
    end

elseif strcmp(method,'eLoreta')
    %%%%%%%%%%%%%%%%%%%% Estimate source using eLoreta %%%%%%%%%%%%%%%%%%%%%%%%
    % Load data files if not already in workspace
    required_files = {'mri_resliced_orig', 'mri_segmented'};
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
    cfg.sourcemodel = leadfield;
    cfg.headmodel   = headmodel;
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
    if displayResults
        % step 2: Plot
        cfg              = [];
        cfg.method       = 'ortho';
        cfg.funparameter = 'pow';
        cfg.funcolormap    = 'jet';
        ft_sourceplot(cfg, data_source_interp);
        % idem, in MNI (or SPM) coordinates
        ft_sourceplot(cfg, sourceDiffIntNorm);
    end
end