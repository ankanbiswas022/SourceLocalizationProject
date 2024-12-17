% showSimulationResults

% Choose subject
subjectName = '204NC'; % '463NR' - has huge fast gamma
folderName = 'N:\Projects\TLSAEEGProject\ftData\'; 

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

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot TF spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(441);
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

cfg = [];
cfg.baseline     = [-0.5 0];
cfg.baselinetype = 'db';
cfg.showlabels   = 'yes';
cfg.layout       = data.elec;
cfg.box          = 'yes';
cfg.colormap     = 'jet';
ft_multiplotTFR(cfg,TFR);

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Topoplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Topoplots
% noseDir = '+X';
% chanlocs = getMontageDetails(refType);
% numGroups = length(subjectNameListFinal);
% for i=1:numGroups
%     axes(hPlots(1+i)); %#ok<LAXES>
%     x = topoData{i};
%     x(isnan(x)) = 999;
%     topoplot_murty(x,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir',noseDir,'emarkercolors',x);
%     caxis(cLims);
%     title(strList{i});
% end
% displayData(hPlots(i,:),subjectNameList,strList,deltaPSD,dataForDisplay.freqVals,topoData,sourceData,dataForDisplay.rangeNames{i},refType,useMedianFlag);