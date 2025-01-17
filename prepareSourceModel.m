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

% methodHeadModel - Method for headmodel. Options: openmeeg, dipoli, bemcp
% methodSourceModel - Method for sourcemodel. Options: basedonresolution, basedonmri

function model = prepareSourceModel(displayModelFlag,capType,methodHeadModel,methodSourceModel)

if ~exist('displayModelFlag','var');  displayModelFlag = 0;             end
if ~exist('capType','var');           capType = 'actiCap64';            end
if ~exist('methodHeadModel','var');   methodHeadModel = 'bemcp';                end
if ~exist('methodSourceModel','var'); methodSourceModel = 'basedonresolution';  end

folderSave = fullfile(pwd,'prepDataSourceModel');
if ~exist(folderSave,'file')
    mkdir(folderSave);
end

%%%%%%%%%%%%%%%%% Get MRI data and boundary meshes %%%%%%%%%%%%%%%%%%%%%%%%

% This portion creates mri_orig.mat, mri_segmented.mat and bnd.mat and saves it
% in the prepDataSourceModel folder. This needs to be done just once, and
% hence is commented out once these files are generated.

% 1. Generate mri_orig.mat
% We can use data already available in fieldtrip toolbox in the folder
% template/anatomy. However, this is not in standard ctf format and needs
% to be realigned. For this realignment procedure, see the tutorial:
% http://www.fieldtriptoolbox.org/workshop/baci2017/forwardproblem/

% Since this realignment is subjective, we instead download a standard MRI dataset from
% ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/Subject01.zip,
% which is in the standard ctf format. Since the raw MRI data is large, we instead save the mri data in mri.mat file

% mri = ft_read_mri('Subject01\Subject01.mri');
% save(fullfile(folderSave,'mri_orig.mat'),'mri');

% 2. Generate mri_segmented.mat
% cfg = [];
% cfg.output    = {'brain','skull','scalp'};
% mri_segmented  = ft_volumesegment(cfg, mri_orig);
% save(fullfile(folderSave,'mri_segmented.mat'),'mri_segmented');

% 3. Create boundary meshes and generate bnd.mat
% cfg = [];
% cfg.tissue = {'brain', 'skull', 'scalp'};
% cfg.numvertices = [1000 1000 1000];
% bnd = ft_prepare_mesh(cfg, mri_segmented);
% save(fullfile(folderSave,'bnd.mat'),'bnd');

%%%%%%%%%%%%%%% Load fixed entities saved above %%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = load(fullfile(folderSave,'mri_segmented.mat'));
mri_segmented = tmp.mri_segmented;

tmp = load(fullfile(folderSave,'bnd.mat'));
bnd = tmp.bnd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Electrode Locations %%%%%%%%%%%%%%%%%%%%%%
fileName = fullfile(folderSave,['electrode_' capType '.mat']);

if exist(fileName,'file')
    tmp = load(fileName);
    elec = tmp.elec;
else
    % Write the code to generate the variable elec from Montage data
    % depending on the Cap type
end

%%%%%%%%%%%%%%%%%%%%%%%%% Create a head model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileName = fullfile(folderSave,['headmodel_' methodHeadModel '.mat']);
if exist(fileName,'file')
    tmp = load(fileName);
    headmodel = tmp.headmodel;
else
    cfg = [];
    thisDir = pwd;
    if strcmp(methodHeadModel,'openmeeg')
        if ispc
            cd('C:\Program Files\OpenMEEG\bin'); % OpenMEEG must be installed from https://openmeeg.github.io/ for this to work.
        end
        cfg.tempdir = 'C:\tmp'; % Create a temporary directory here since directories with blanks do not work
    end
    
    cfg.method = methodHeadModel;
    cfg.tissue = {'brain', 'skull', 'scalp'};
    headmodel = ft_prepare_headmodel(cfg, bnd);
    cd(thisDir);
    save(fileName,'headmodel');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Create a source model %%%%%%%%%%%%%%%%%%%%%%%%%
fileName = fullfile(folderSave,['sourcemodel_' methodHeadModel '_' methodSourceModel '.mat']); % Name includes headmodel also since the headmodel is also given as an input
if exist(fileName,'file')
    tmp = load(fileName);
    sourcemodel = tmp.sourcemodel;
else
    cfg = [];
    cfg.grid.unit = 'mm';
    cfg.threshold = 0.1;
    cfg.smooth = 5;
    cfg.inwardshift = 1;
    cfg.method = methodSourceModel;
    cfg.headmodel = headmodel;

    if strcmp(methodSourceModel,'basedonresolution')
        cfg.grid.resolution = 5;
    elseif strcmp(methodSourceModel,'basedonmri')        
        cfg.mri = mri_segmented;
    end
    sourcemodel = ft_prepare_sourcemodel(cfg);
    save(fileName,'sourcemodel');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute lead field %%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = fullfile(folderSave,['leadfield_' capType '_' methodHeadModel '_' methodSourceModel '.mat']);
if exist(fileName,'file')
    tmp = load(fileName);
    leadfield = tmp.leadfield;
else
    cfg = [];
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel = headmodel;
    cfg.elec = elec;
    cfg.channel = elec.label;
    cfg.reducerank = 3;
    cfg.normalize = 'no';

    % Change directory for openmeeg
    thisDir = pwd;
    if strcmp(methodHeadModel,'openmeeg')
        if ispc
            cd('C:\Program Files\OpenMEEG\bin'); % OpenMEEG must be installed from https://openmeeg.github.io/ for this to work.
        end
        cfg.tempdir = 'C:\tmp'; % Create a temporary directory here since directories with blanks do not work
    end
    leadfield = ft_prepare_leadfield(cfg);
    cd(thisDir);
    save(fileName,'leadfield');
end

% preprae the model structure
model.mesh        = bnd;
model.elec        = elec;      
model.headmodel   = headmodel;
model.sourcemodel = sourcemodel;
model.leadfield   = leadfield;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayModelFlag
    
    subplot(121)
    ft_plot_headmodel(headmodel, 'facecolor', 'skin', 'facealpha', 0.7);
    hold on;
    ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'on');
    camlight;
    title('Headmodel and Electrodes');

    subplot(122)
    ft_plot_mesh(bnd(1), 'facecolor','none');
    ft_plot_mesh(bnd(2), 'facecolor','none');
    ft_plot_mesh(bnd(3), 'facecolor','none');
    ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:), 'vertexmarker', 'o', 'vertexsize', 5, 'vertexcolor', [1 0 0]);
    title('Meshes and SourceModel');
end