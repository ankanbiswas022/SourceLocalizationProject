% showSimulationResults

% Choose subject
subjectName = '204NC'; % '463NR' - has huge fast gamma
folderName  = 'N:\Projects\TLSAEEGProject\ftData\'; 

%%%%%%%%%%%%%%%%%%%% Get signal in time domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectName = 'ADGammaProject';
protocolType = 'SF_ORI';
freqRange = [20 35]; % slow gamma freq range

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
hfig=getPlotHandles(4,4,[0.1 0.1 0.8 0.8],0.05,0.05);
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
cfg.rotate       = 90;
cfg.box          = 'yes';
cfg.colormap     = 'jet';
% ft_multiplotTFR(cfg,TFR);


%%%%%%%%%%%%%%%%%%%%%%%%%% Plot TF plot(high priority electrodes) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(hfig(1,1));
cfg.figure = gca;
gridType = 'EEG';
capType  = 'actiCap64';
load([capType 'Labels.mat']);
[~,~,~,~,~,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType);
selectedChannels = montageLabels(highPriorityElectrodeNums,2);
cfg.channel = selectedChannels;
cfg.xlim         = [-0.5 0.75];
cfg.ylim         = [0 100];
ft_singleplotTFR(cfg, TFR);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ch= colorbar;
ch.Label.String = '\Delta Power (db)';

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Topoplot (high priority electrodes) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(hfig(1,2));
cfg.figure = gca;
cfg.channel = 'all';
cfg.fontsize = 8;
cfg.xlim = [0.25 0.75]; %st period
cfg.ylim = freqRange; 
cfg.interactive = 'no';
cfg.comment =  'auto';
cfg.shading = 'interp';
cfg.commentpos = [-1 0.7];
ft_topoplotTFR(cfg, TFR); 
ch = colorbar;
ch.Label.String = '\Delta Power (db)';