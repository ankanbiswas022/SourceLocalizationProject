function displaySourceData(sourceData,method,seedVoxel,hAxis,titleStr)

if ~exist('method','var');          method = 'ortho';                   end
if ~exist('seedVoxel','var');       seedVoxel = [];                     end
if ~exist('hAxis','var');           hAxis = [];                         end
if ~exist('titleStr','var');        titleStr = '';                      end

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

cfg = [];
cfg.funparameter = 'pow';
cfg.funcolorlim = 'maxabs';
cfg.funcolormap = 'jet';        % Use jet colormap for better contrast
cfg.funcolorbar = 'yes';

cfg.method = method;

if ~isempty(hAxis)
    cfg.figure = gcf;
    axes(hAxis);
end
if ~isempty(titleStr)
    cfg.title = titleStr;
end

if strcmp(method,'ortho')
    if ~isempty(seedVoxel)
        cfg.location = seedVoxel;
        cfg.locationcoordinates = 'head';
    end

    ft_sourceplot(cfg, sourceDataInterp);

elseif strcmp(method,'slice')
    if ~isempty(seedVoxel)
        cfg.slicerange = [seedVoxel(3) seedVoxel(3)];
        cfg.nslices = 1;
    else
        cfg.nslices = 9;
    end

    ft_sourceplot(cfg, sourceDataInterp);

elseif strcmp(method,'surface')
    cfg2 = [];
    cfg2.nonlinear = 'no';
    sourceDataSurface = ft_volumenormalise(cfg2, sourceDataInterp); % converts to MNI coordinates

    cfg.maskparameter  = cfg.funparameter;
    cfg.opacitymap     = 'rampup';
    cfg.projmethod     = 'nearest';
    cfg.surffile       = 'surface_white_both.mat';
    cfg.surfdownsample = 10;

    ft_sourceplot(cfg, sourceDataSurface);
    view ([-90 0]);
end
end