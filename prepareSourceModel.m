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

% Todo: the method can become string to have different methods for
% headmodel, sourcemodel generation
% and each method can have configuration structure

function model = prepareSourceModel(capType,method,displayModelFlag)

if ~exist('capType','var');          capType = 'actiCap64';                 end
if ~exist('method','var');           method = {'openmeeg','resolution'};    end
if ~exist('displayModelFlag','var'); displayModelFlag = 0;                  end

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
    elecLabel = elec.label;
else
    % Write the code to generate the variable elec from Montage data
    % depending on the Cap type
end

%%%%%%%%%%%%%%%%%%%%%%%%% Create a head model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = fullfile(folderSave,['headmodel_' method{1} '.mat']);
if exist(fileName,'file')
    tmp = load(fileName);
    headmodel = tmp.headmodel;
else
    thisDir = pwd;
    if strcmp(method{1},'openmeeg')
        cd('C:\Program Files\OpenMEEG\bin'); % OpenMEEG must be installed from https://openmeeg.github.io/ for this to work.
        method = 'openmeeg';
    elseif strcmp(method{1},'dipoli')
        method = 'dipoli';
    end
    cfg = [];
    cfg.method = method;
    cfg.tissue = {'brain', 'skull', 'scalp'};
    headmodel = ft_prepare_headmodel(cfg, bnd);
    save(fileName,'headmodel');
    cd(thisDir);
end

% Create resolution-based source model
fileName = fullfile(folderSave,['sourcemodel_' method{2} '.mat']);
if exist(fileName,'file')
    tmp = load(fileName);
    sourcemodel = tmp.sourcemodel;
else
    if strcmp(method{2},'resolution')
        cfg = [];
        cfg.grid.resolution = 10;
        cfg.grid.unit = 'mm';
        cfg.threshold = 0.1;
        cfg.smooth = 5;
        cfg.headmodel = headmodel;
        cfg.inwardshift = 1;
        sourcemodel = ft_prepare_sourcemodel(cfg);
    elseif strcmp(method{2},'fmri')
        % Create fMRI-based source model
        tmp = load('mri_segmented.mat');
        cfg = [];
        cfg.method = 'basedonmri';
        cfg.grid.resolution = 10;
        cfg.grid.unit = 'mm';
        cfg.threshold = 0.1;
        cfg.smooth = 5;
        cfg.headmodel = headmodel;
        cfg.inwardshift = 1;
        cfg.mri = tmp.mri_segmented;
        sourcemodel = ft_prepare_sourcemodel(cfg);
    end
    save(fileName,'sourcemodel');
end

% Compute lead field
fileName = fullfile(folderSave,['leadfield_' method{1} '.mat']);
if exist(fileName,'file')
    tmp = load(fileName);
    leadfield = tmp.leadfield;
else
    cfg = [];
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel = headmodel;
    cfg.elec = elec;
    cfg.channel = elecLabel;
    cfg.reducerank = 3;
    cfg.normalize = 'no';
    leadfield = ft_prepare_leadfield(cfg);
    save(fileName,'leadfield');
end

% preprae the model structure
model.sourcemodel = sourcemodel;
model.leadfield   = leadfield;
model.headmodel   = headmodel;
model.mesh        = bnd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayModelFlag
    subplot(221)
    ft_plot_mesh(bnd(1), 'facecolor','none');
    ft_plot_mesh(bnd(2), 'facecolor','none');
    ft_plot_mesh(bnd(3), 'facecolor','none');
    hold on;
    ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'on');
    title('Mesh structure');

    subplot(222)
    ft_plot_headmodel(headmodel, 'facecolor', 'skin', 'facealpha', 0.7);
    hold on;
    ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'on');
    camlight;
    title('Headmodel');

    subplot(223)
    ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:), 'vertexmarker', 'o', 'vertexsize', 5);
    title('SourceModel');
end