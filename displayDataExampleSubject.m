% displayDataExampleSubject

stRange = [0.25 0.75];
blRange = [-diff(stRange) 0];
freqRange = [20 34]; % slow gamma freq range

% Choose subject
subjectName = '204NC'; % '463NR' - has huge fast gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = fullfile('prepDataSourceModel',[subjectName '.mat']);

if exist(fileName,'file')
    tmp = load(fileName);
    data = tmp.data;
else
    % If the file does not exist, it is created and saved. Need access to
    % the original data for this

    folderName  = 'N:\Projects\TLSAEEGProject\ftData\';

    %%%%%%%%%%%%%%%%%%%% Get signal in time domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    projectName = 'ADGammaProject';
    protocolType = 'SF_ORI';

    [expDates,protocolNames,capLayout,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);

    %%%%%%% Instead of concatenating, just find the first good protocol %%%%%%%
    iProt = 1;
    while(1)
        x = load(fullfile(folderName,projectName,protocolType,[subjectName '-' expDates{iProt} '-' protocolNames{iProt}]));
        if (x.goodProtFlag)
            data = x.data;
            break;
        else
            iProt = iProt+1;
        end
    end
    save(fileName,'data');
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot TF spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hfig = getPlotHandles(2,3,[0.1 0.1 0.8 0.8],0.05,0.05);

%%%%%%%%%%%%%%%%%%%%%%% Time-frequency representation %%%%%%%%%%%%%%%%%%%%%
windowLenS       = 0.25; % Seconds
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = (1/windowLenS):(1/windowLenS):100; 
cfg.t_ftimwin    = ones(length(cfg.foi),1).* windowLenS;   % length of time window
cfg.toi          = -0.5:0.05:0.75;

TFR = ft_freqanalysis(cfg,data);

%%%%%%%%%%%%%%%%%%% Plot TF plot(high priority electrodes) %%%%%%%%%%%%%%%%
gridType = 'EEG';
capType  = 'actiCap64';
load([capType 'Labels.mat']);
[~,~,~,~,~,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType);
selectedChannels = montageLabels(highPriorityElectrodeNums,2);

axes(hfig(1,1));

cfg = [];
cfg.figure       = gca;
cfg.baseline     = blRange;
cfg.baselinetype = 'db';
cfg.showlabels   = 'yes';
cfg.layout       = data.elec;
cfg.rotate       = 90;
cfg.box          = 'yes';
cfg.colormap     = 'jet';
cfg.channel      = selectedChannels;
cfg.xlim         = [-0.5 0.75];
cfg.ylim         = [0 100];
ft_singleplotTFR(cfg, TFR);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ch= colorbar;
ch.Label.String = '\Delta Power (db)';

%%%%%%%%%%%% Plot Topoplot (high priority electrodes) %%%%%%%%%%%%%%%%%%%%%

axes(hfig(2,1));
cfg.figure = gca;
cfg.channel = 'all';
cfg.fontsize = 8;
cfg.xlim = stRange; %st period
cfg.ylim = freqRange; 
cfg.interactive = 'no';
cfg.comment =  'auto';
cfg.shading = 'interp';
cfg.commentpos = [-0.5 0.7];
ft_topoplotTFR(cfg, TFR); 
ch = colorbar;
ch.Label.String = '\Delta Power (db)';

%%%%%%%%%%%%%%% Source localization using LORETA software %%%%%%%%%%%%%%%%%
% Note that the number of trials is more since all good protocols are used.
% 
% fileName = fullfile('prepDataSourceModel',[subjectName '_LORETAKey.mat']);
% 
% if exist(fileName,'file')
%     tmp = load(fileName);
%     sourceData = tmp.sourceData;
% else
%     dataStr = 'mid';
%     folderLORETA = 'N:\Projects\Kanishka_SourceLocalizationProject\data\sLORETA_Thres10';
%     sourceData = load(fullfile(folderLORETA,dataStr,[subjectName '.mat']));
%     save(fileName,'sourceData');
% end
% 
% freqRangePos = 2; % 1 - alpha, 2, - SG, 3 - FG
% cLims = [-1 2];
% mData = 10*(log10(sourceData.mDataST) - log10(sourceData.mDataBL));
% dataToPlot = mData(freqRangePos,:);
% 
% [posList,xyz,areaList] = getVoxelInfo;
% axes(hfig(1,3));
% scatter3(xyz(:,1),xyz(:,2),xyz(:,3),10,dataToPlot); 
% clim(cLims);
% colorbar

%%%%%%%%%%%%%%% Source localization using fieldtrip %%%%%%%%%%%%%%%%%%%%%%%
methodSourceLoc = 'eloreta';

capType = 'actiCap64';
methodHeadModel = 'bemcp';
methodSourceModel = 'basedonresolution';
commonFilterFlag = 1;

model = prepareSourceModel(0,capType,methodHeadModel,methodSourceModel);
[dataPre,dataPost] = prepareDataSourceLocalization(data,stRange);
[sourceDataDiff,sourceDataPre,sourceDataPost] = performSourceLocalization(dataPre,dataPost,0,methodSourceLoc,capType,methodHeadModel,methodSourceModel,commonFilterFlag,freqRange);

[powerVsDistanceOriginal,distanceBinsOriginal,coeffsOriginal,distanceList,seedVoxel] = getPowerVsDistance(sourceDataDiff);

displaySourceData(sourceDataPre,'slice',seedVoxel{1},hfig(1,2),'Pre');
displaySourceData(sourceDataPost,'slice',seedVoxel{1},hfig(2,2),'Post');
displaySourceData(sourceDataDiff,'slice',seedVoxel{1},hfig(1,3),'Change (dB)');

axes(hfig(2,3))
plot(distanceBinsOriginal{1},powerVsDistanceOriginal{1},'ko');

axes(hfig(1,2)); title('baseline');
axes(hfig(2,2)); title('stimulus');
axes(hfig(1,3)); title('Change (dB)');
axes(hfig(2,3)); 
xlabel('Distance from peak voxel (mm)'); 
ylabel('Change in power (dB)');