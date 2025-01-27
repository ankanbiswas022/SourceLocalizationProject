% This program performs source localization using fieldtrip

function [sourceData,sourcePre,sourcePost] = performSourceLocalization(dataPre,dataPost,displayResultsFlag,methodSourceLoc,capType,methodHeadModel,methodSourceModel,commonFilterFlag,freqRange)

if ~exist('displayResultsFlag','var');  displayResultsFlag = 0;         end
if ~exist('methodSourceLoc','var');     methodSourceLoc = 'dics';       end
if ~exist('capType','var');             capType = 'actiCap64';          end
if ~exist('methodHeadModel','var');     methodHeadModel = 'bemcp';                end
if ~exist('methodSourceModel','var');   methodSourceModel = 'basedonresolution';  end
if ~exist('commonFilterFlag','var');    commonFilterFlag = 1;           end
if ~exist('freqRange','var');           freqRange = [20 34];            end

%%%%%%%%%%%%%%%%%%%%%%%%%% Get Source Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = prepareSourceModel(0,capType,methodHeadModel,methodSourceModel);

% Perform cross spectral analysis
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.taper     = 'dpss';
cfg.tapsmofrq = diff(freqRange)/2;
cfg.foi       = mean(freqRange);

freqPre       = ft_freqanalysis(cfg, dataPre);
freqPost      = ft_freqanalysis(cfg, dataPost);

if commonFilterFlag
    dataAll = ft_appenddata([], dataPre, dataPost);
    freqAll = ft_freqanalysis(cfg, dataAll);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Source Localization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.method = methodSourceLoc;
cfg.grid = model.leadfield;
cfg.headmodel = model.headmodel;
cfg.frequency = freqPre.freq;

if strcmp(methodSourceLoc,'dics') % Special parameters for dics
    cfg.dics.projectnoise = 'yes';
    cfg.dics.lambda = '2%';
    cfg.dics.keepfilter = 'yes';
    cfg.dics.realfilter = 'yes';

    if commonFilterFlag  % Compute a common filter for both pre and post
        sourceAll = ft_sourceanalysis(cfg, freqAll);
        cfg.sourcemodel.filter = sourceAll.avg.filter; % Common filter
    end
end

sourcePre = ft_sourceanalysis(cfg, freqPre);
sourcePost = ft_sourceanalysis(cfg, freqPost);

% Calculate normalized source power
sourceData = sourcePost;
%   sourceData.avg.pow = (sourcePost.avg.pow - sourcePre.avg.pow) ./ sourcePre.avg.pow;
sourceData.avg.pow = 10*(log10(sourcePost.avg.pow) - log10(sourcePre.avg.pow)); % change in power in decibels

if displayResultsFlag
    displaySourceData(sourceData);
end
end