% This program returns average power binned as a function of distance from
% given (seed) voxels. Also fits an exponential to the binned data.

% sourceData - contains positions of sources and also change in power at
% each source location. Returned by performSourceLocalization.

% seedVoxels - one voxel in each hemifield. If left empty, the program
% finds the one which has the highest absolute power.

function [powerVsDistance,distanceBins,coeffs] = getPowerVsDistance(sourceData,seedVoxels,displayResultsFlag)

if ~exist('seedVoxels','var');          seedVoxels = [];                end
if ~exist('displayResultsFlag','var');  displayResultsFlag = 0;         end

% Find voxel positions inside the brain and power values
goodPos = sourceData.pos(sourceData.inside,:); % All voxels inside the brain
goodPowerVals = sourceData.avg.pow(sourceData.inside,:); % Power values at these voxels

% Find voxels for the two hemifields. Midline voxels can be assigned to one
% of the two groups if needed
indexList{1} = find(goodPos(:,2)>0);
indexList{2} = find(goodPos(:,2)<0);
goodPosList{1} = goodPos(indexList{1},:);
goodPosList{2} = goodPos(indexList{2},:);

% Find seed voxels if not specified
if isempty(seedVoxels)
    tmpSeedVoxels = cell(1,2);
    for i=1:2
        vals = goodPowerVals(indexList{i});
        pos = (abs(vals)==max(abs(vals)));
        tmpSeedVoxels{i} = goodPosList{i}(pos,:);
    end
    seedVoxels = tmpSeedVoxels;
end

% Find euclidean distanceBins from seed voxels on both sides, compute binned
% power and fit exponentials
distanceList = cell(1,2);
powerVsDistance = cell(1,2);
distanceBins = cell(1,2);
coeffs = zeros(2,2);
for i=1:2
    posDifference = goodPosList{i} - seedVoxels{i};
    distanceList{i} = sqrt(posDifference(:,1).^2 + posDifference(:,2).^2 + posDifference(:,3).^2);
    [powerVsDistance{i},distanceBins{i},coeffs(i,:)] = getBinnedData(goodPowerVals(indexList{i}),distanceList{i});
end

if displayResultsFlag
    titleStr = [{'Left'} {'Right'}];
    % Figure 1 - powerVsDistance
    figure;
    for i=1:2
        subplot(1,2,i)
        plot(distanceList{i},goodPowerVals(indexList{i}),'.','color',[0.7 0.7 0.7]);
        hold on;
        plot(distanceBins{i},powerVsDistance{i},'ko','markerfacecolor','k');
    end

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

    for i=1:2
        % Visualize interpolated source data
        cfg = [];
        cfg.method = 'ortho';
        cfg.funparameter = 'pow';
        cfg.funcolorlim = 'maxabs';
        cfg.funcolormap = 'jet';        % Use jet colormap for better contrast

        cfg.title = titleStr{i};
        cfg.location = seedVoxels{i};
        cfg.locationcoordinates = 'head';

        ft_sourceplot(cfg, sourceDataInterp);
    end
end
end

function [powerVsDistance,distanceBins,coeffs] = getBinnedData(powerVals,distanceList)

% Average power across bins
binWidth = 20; binLimit = binWidth*ceil(max(distanceList)/binWidth);
binEdges = 0:binWidth:binLimit;
numBins = length(binEdges)-1;

distanceBins = zeros(1,numBins);
powerVsDistance = zeros(1,numBins);

for i=1:numBins
    pos = intersect(find(distanceList>binEdges(i)),find(distanceList<=binEdges(i+1)));
    powerVsDistance(i) = mean(powerVals(pos));
    distanceBins(i) = (binEdges(i) + binEdges(i+1))/2;
end

% Fit exponential to the binned data
coeffs = [0 0];
end