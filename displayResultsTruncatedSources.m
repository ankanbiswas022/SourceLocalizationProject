subjectName = '204NC';
seedVoxel = [-20 10 70]; % obtained from displayDataExampleSubject

tmp = load(fullfile('prepDataSourceModel',[subjectName '.mat']));
data = tmp.data;

methodSourceLocFreq = 'eloreta';
methodSourceLocTime = 'eloreta';

capType = 'actiCap64';
methodHeadModel = 'bemcp';
methodSourceModel = 'basedonresolution';
commonFilterFlag = 1;
stRange = [0.25 0.75];
freqRange = [20 34];

%%%%%%%%%%%%%%%%%%%%%%%%%% Get Source Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = prepareSourceModel(0,capType,methodHeadModel,methodSourceModel);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduce number of trials %%%%%%%%%%%%%%%%%%%%%%%
cfg                 = [];
N                   = 10;
cfg.trials          = [true(1,N) false(1,data.numTrials - N)];
dataShort           = ft_selectdata(cfg,data);
dataShort.numTrials = N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dataSegmentPre,dataSegmentPost] = prepareDataSourceLocalization(dataShort,stRange);

%%%%%%%%%%%%%%% Source localization in frequency domain %%%%%%%%%%%%%%%%%%%
[sourceDataDiffFreq,sourceDataPreFreq,sourceDataPostFreq] = performSourceLocalization(dataSegmentPre,dataSegmentPost,0,methodSourceLocFreq,capType,methodHeadModel,methodSourceModel,commonFilterFlag,freqRange);

%%%%%%%%%%%% Find power values as a function of seed voxel %%%%%%%%%%%%%%%%
seedVoxels{1} = seedVoxel;
[powerVsDistanceOriginal,distanceBinsOriginal,coeffsOriginal,distanceList,seedVoxelTruncated] = getPowerVsDistance(sourceDataDiffFreq,seedVoxels);

%%%%%%%%%%%%%%%%%% Display topo and source plots %%%%%%%%%%%%%%%%%%%%%%%%%%
displayTruncatedDataFlag = 0;

if displayTruncatedDataFlag

    figure;

    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'pow';
    cfg.taper     = 'dpss';
    cfg.tapsmofrq = diff(freqRange)/2;
    cfg.foi       = mean(freqRange);

    freqDataPre = ft_freqanalysis(cfg, dataSegmentPre);
    freqDataPost = ft_freqanalysis(cfg, dataSegmentPost);

    subplot(2,3,1);
    cfg = [];
    cfg.figure = gca;
    cfg.colormap ='jet';
    cfg.colorbar = 'yes';
    cfg.rotate = 90;
    freqOriginal = freqDataPre;
    freqOriginal.powspctrm = 10*(log10(freqDataPost.powspctrm) - log10(freqDataPre.powspctrm));
    ft_topoplotER(cfg, freqOriginal);

    subplot(2,3,2); displaySourceData(sourceDataDiffFreq,'slice',seedVoxel,gca,'Change (dB)');
    subplot(2,3,4); displaySourceData(sourceDataPreFreq,'slice',seedVoxel,gca,'baseline');
    subplot(2,3,5); displaySourceData(sourceDataPostFreq,'slice',seedVoxel,gca,'stimulus');

    subplot(1,3,3);
    plot(distanceBinsOriginal{1},powerVsDistanceOriginal{1},'ko');
    xlabel('Distance from seed voxel'); ylabel('Change in power (dB)');
end

%%%%%%%%%%%%%%% Source localization in time domain %%%%%%%%%%%%%%%%%%%%%%%%
% Get covariance matrix for the entire segment
cfg               = [];
cfg.covariance    = 'yes';
cfg.keeptrials    = 'yes';
timelockPre       = ft_timelockanalysis(cfg, dataSegmentPre);
timelockPost      = ft_timelockanalysis(cfg, dataSegmentPost);

% Source time series for individual trials
cfg                 = [];
cfg.method          = methodSourceLocTime;
cfg.grid            = model.leadfield;
cfg.headmodel       = model.headmodel;
cfg.rawtrial        = 'yes';

sourceSegmentTimePre  = ft_sourceanalysis(cfg, timelockPre);
sourceSegmentTimePost = ft_sourceanalysis(cfg, timelockPost);

%%%%%%%%%%%% Signal reconstruction from a subject of voxels %%%%%%%%%%%%%%%
cutoffDistanceList = [100 80 60 40];
nCutoffs = length(cutoffDistanceList);

signalToShow = 'post'; % 'post' or 'diff'

for i=1:nCutoffs

    %%%% Find voxels within some distance from the voxel with peak power %%%%%%
    cutoffDistance = cutoffDistanceList(i); % mm
    goodPos = find(sourceDataDiffFreq.inside);
    badVoxelPos = goodPos(distanceList{1}>cutoffDistance);
    goodVoxelPos = goodPos(distanceList{1}<=cutoffDistance);

    %%%%%%%%%%%%%%%%%%%% Restrict sources within the cutoff %%%%%%%%%%%%%%%
    if strcmp(signalToShow,'diff')
        sourceDataFreqTruncated = sourceDataDiffFreq;
    elseif strcmp(signalToShow,'pre')
        sourceDataFreqTruncated = sourceDataPreFreq;
    elseif strcmp(signalToShow,'post')
        sourceDataFreqTruncated = sourceDataPostFreq;
    end
    sourceDataFreqTruncated.avg.pow(badVoxelPos) = 0;

    subplot(nCutoffs,4,4*(i-1)+1); displaySourceData(sourceDataFreqTruncated,'slice',seedVoxel,gca);

    %%%%%%%%%%%%%%%%%%%%%% Recreate signal from sources %%%%%%%%%%%%%%%%%%%%%%%

    dataSegmentTruncatedPre = dataSegmentPre;
    sensorDataPre = source2sensor(sourceSegmentTimePre,model.leadfield,goodVoxelPos);
    dataSegmentTruncatedPre.trial = sensorDataPre;
    dataSegmentTruncatedPre.label = model.leadfield.label;

    dataSegmentTruncatedPost = dataSegmentPost;
    sensorDataPost = source2sensor(sourceSegmentTimePost,model.leadfield,goodVoxelPos);
    dataSegmentTruncatedPost.trial = sensorDataPost;
    dataSegmentTruncatedPost.label = model.leadfield.label;

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'pow';
    cfg.taper     = 'dpss';
    cfg.tapsmofrq = diff(freqRange)/2;
    cfg.foi       = mean(freqRange);

    freqDataTruncatedPre = ft_freqanalysis(cfg, dataSegmentTruncatedPre);
    freqDataTruncatedPost = ft_freqanalysis(cfg, dataSegmentTruncatedPost);

    subplot(nCutoffs,4,4*(i-1)+2);
    cfg = [];
    cfg.figure = gca;
    cfg.colormap ='jet';
    cfg.colorbar = 'yes';
    cfg.rotate = 90;
    freqTruncated = freqDataTruncatedPre;

    if strcmp(signalToShow,'diff')
        freqTruncated.powspctrm = 10*(log10(freqDataTruncatedPost.powspctrm) - log10(freqDataTruncatedPre.powspctrm));
    elseif strcmp(signalToShow,'pre')
        freqTruncated.powspctrm = log10(freqDataTruncatedPre.powspctrm);
    elseif strcmp(signalToShow,'post')
        freqTruncated.powspctrm = log10(freqDataTruncatedPost.powspctrm);
    end

    ft_topoplotER(cfg, freqTruncated);

    %%%%%%%%%%%%%%%%% Source localization %%%%%%%%%%%%%%%%%%%%%
    [sourceDataDiffTruncated2,sourceDataPreTruncated2,sourceDataPostTruncated2] = performSourceLocalization(dataSegmentTruncatedPre,dataSegmentTruncatedPost,0,methodSourceLocFreq,capType,methodHeadModel,methodSourceModel,commonFilterFlag,freqRange);

    if strcmp(signalToShow,'diff')
        sourceDataFreqTruncated2 = sourceDataDiffTruncated2;
    elseif strcmp(signalToShow,'pre')
        sourceDataFreqTruncated2 = sourceDataPreTruncated2;
    elseif strcmp(signalToShow,'post')
        sourceDataFreqTruncated2 = sourceDataPostTruncated2;
    end

    subplot(nCutoffs,4,4*(i-1)+3); displaySourceData(sourceDataFreqTruncated2,'slice',seedVoxel,gca);

    subplot(nCutoffs,4,4*i);
    [powerVsDistanceTruncated,distanceBinsTruncated,coeffsTruncated] = getPowerVsDistance(sourceDataFreqTruncated,seedVoxels);
    [powerVsDistanceTruncated2,distanceBinsTruncated2,coeffsTruncated2] = getPowerVsDistance(sourceDataFreqTruncated2,seedVoxels);

    plot(distanceBinsTruncated{1},powerVsDistanceTruncated{1},'ko'); hold on;
    plot(distanceBinsTruncated2{1},powerVsDistanceTruncated2{1},'r*'); hold on;
end