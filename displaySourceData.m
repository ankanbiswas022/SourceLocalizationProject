function displaySourceData(sourceData,method,seedVoxel)

if ~exist('method','var');          method = 'ortho';                   end
if ~exist('seedVoxel','var');       seedVoxel = [];                     end

% Source Data Visualization - ortho
% Interpolate source data onto MRI for visualization
tmp = load('mri_orig.mat');
mri_orig = tmp.mri;
cfg = [];
mri_resliced_orig = ft_volumereslice(cfg, mri_orig);

cfg = [];
cfg.downsample = 2;
cfg.parameter = 'pow';
sourceDataInterp = ft_sourceinterpolate(cfg, sourceData, mri_resliced_orig);

if strcmp(method,'ortho')
    % Visualize interpolated source data
    cfg = [];
    cfg.method = 'ortho';
    cfg.funparameter = 'pow';
    cfg.funcolorlim = 'maxabs';
    cfg.funcolormap = 'jet';        % Use jet colormap for better contrast
    
    if ~isempty(seedVoxel)
        cfg.location = seedVoxel;
        cfg.locationcoordinates = 'head';
    end

    ft_sourceplot(cfg, sourceDataInterp);

elseif strcmp(method,'surface')

    % Source Data Visualization - surface
    cfg = [];
    cfg.nonlinear = 'no';
    sourceDataSurface = ft_volumenormalise(cfg, sourceDataInterp); % converts to MNI coordinates

    cfg = [];
    cfg.method         = 'surface';
    cfg.funparameter   = 'pow';
    cfg.maskparameter  = cfg.funparameter;
    cfg.funcolormap    = 'jet';
    cfg.opacitymap     = 'rampup';
    cfg.projmethod     = 'nearest';
    cfg.surffile       = 'surface_white_both.mat';
    cfg.surfdownsample = 10;
    ft_sourceplot(cfg, sourceDataSurface);
    view ([90 0]);
end
end