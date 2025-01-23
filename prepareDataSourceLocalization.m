% This program prepares data for source localization. 

function [dataPre,dataPost] = prepareDataSourceLocalization(data,stRange)

if ~exist('stRange','var');             stRange = [0.25 0.75];          end

%%%%%%%%%%%%%%%%%%%%%%% Data Preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Interpolate the bad electrodes

% Get bad channel names from indices
badChannels = data.label(data.badElecs);

% % Prepare neighbours using distance method
% % note: requires the neighours file witin fieldtrip/tempate directory
% cfg = [];
% cfg.method = 'distance';    % Use distance method
% cfg.neighbourdist = 50;     % Maximum distance between neighbours
% cfg.elec = model.elec;
% cfg.feedback = 'no';       % Show feedback about the neighbours
% neighbours = ft_prepare_neighbours(cfg);

% Prepare neighbours using template method - gives slightly better results
% note: requires the neighours file witin fieldtrip/tempate directory
cfg = [];
cfg.method = 'template';    % Use template method instead of spline
cfg.template = 'elec1010';  % Using standard 10-10 template
cfg.layout = 'elec1010';    % Specify the layout file
cfg.feedback = 'no';       % Show feedback about the neighbours
neighbours = ft_prepare_neighbours(cfg);

% interpolate
cfg = [];
cfg.method = 'spline';      % Keep spline for interpolation
cfg.badchannel = badChannels;  % Using the extracted bad channel names
cfg.neighbours = neighbours;    % Use the prepared neighbours
data = ft_channelrepair(cfg, data);

% 2. Perform average reference
cfg = [];
cfg.demean = 'yes';         % Remove mean value
cfg.reref = 'yes';          % Re-reference the data
cfg.refmethod = 'avg';      % Average reference
cfg.refchannel = 'all';     % Use all channels for referencing
data = ft_preprocessing(cfg, data); % average reference the data

% 3. Prepare baseline and stimulus segments for analysis
cfg = [];
cfg.toilim = [-diff(stRange) 0];
dataPre = ft_redefinetrial(cfg, data);
cfg.toilim = stRange;
dataPost = ft_redefinetrial(cfg, data);

end