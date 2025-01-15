% EEG Source Analysis Pipeline
% Dependency: FieldTrip toolbox
% This script generates the data (headmodel, sourcemodel and leadfield)
% required to perform source analysis. This includes the following:
% - MRI preprocessing
% - Head model creation
% - Source model generation
% - Lead field computation

% Each of these steps are associated with many variables (e.g. number of
% vertices, spacing between sources etc. However, only a single set of
% models is saved for each cap. If you change the variables in the code,
% please remove the folder and let the code regenerate the model variables.

function prepareSourceModel(capType,displayModelFlag)

if ~exist('capType','var');          capType = [];                      end
if ~exist('displayModelFlag','var'); displayModelFlag = 0;              end

if isempty(capType);                capType = 'actiCap64';              end

folderSave = fullfile(pwd,'prepDataSourceModel',capType); % All data is saved here
if ~exist(folderSave,'file')
    mkdir(folderSave);
end

%%%%%%%%%%%%%%%%%% Get boundary meshes if not available %%%%%%%%%%%%%%%%%%%
fileName = fullfile(folderSave,'bnd.mat');

if exist(fileName,'file')
    tmp = load(fileName);
    bnd = tmp.bnd;
else

    % Load and preprocess MRI data
    tmp = load('mri_segmented.mat');
    mri_segmented = tmp.mri_segmented;
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
    % depending on the Cap type
end

%%%%%%%%%%%%%%%%%%%%%%%%% Create a head model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thisDir = pwd;
if ispc
    cd('C:\Program Files\OpenMEEG\bin'); % OpenMEEG must be installed from https://openmeeg.github.io/ for this to work.
    method = 'openmeeg';
else
    method = 'dipoli';
end

fileName = fullfile(folderSave,['headmodel_' method '.mat']);
if exist(fileName,'file')
    tmp = load(fileName);
    headmodel = tmp.headmodel;
else
    
    cfg = [];
    cfg.method = method;
    cfg.tissue = {'brain', 'skull', 'scalp'};
    headmodel = ft_prepare_headmodel(cfg, bnd);
    save(fileName,'headmodel');
end
cd(thisDir);

% Create resolution-based source model
cfg = [];
cfg.grid.resolution = 5;
cfg.grid.unit = 'mm';
cfg.threshold = 0.1;
cfg.smooth = 5;
cfg.headmodel = headmodel;
cfg.inwardshift = 1;
sourcemodel = ft_prepare_sourcemodel(cfg);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayModelFlag
    subplot(221)
    ft_plot_mesh(bnd(1), 'facecolor','none');
    ft_plot_mesh(bnd(2), 'facecolor','none');
    ft_plot_mesh(bnd(3), 'facecolor','none');
    hold on;
    ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'on');

    subplot(222)
    ft_plot_headmodel(headmodel_bem, 'facecolor', 'skin', 'facealpha', 0.7);
    hold on;
    ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'on');
    camlight;

    subplot(223)
    ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:), 'vertexmarker', 'o', 'vertexsize', 5);
end