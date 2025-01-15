% ToDo
% Display
% 0. Display - show better brain maps. Currently using scatter3.
% 1. Display - 8th column - plot errorbars properly as swarm plots.
% 2. Display - the 9th column show show fraction as a function from the
% voxel that has maximum activation. 

% Analysis
% 1. Subject power matching
% 2. Find electrodes that are bad in many subjects and remove those common
% bad electrodes. Do LORETA analysis on the subset of subjects who have remaining
% good electrodes
% 3. Number of trials matching - since stats may depend on the number of
% trials, control for that by having the same number of trials for all subjects.
% 4. If there are enough subjects for Case/Control - do that comparison.

clear; clc;

thres=-inf; % only subjects who have delta power above this (in dB) are selected. Set to a low value such as -inf to take all subjects
useCommonSubjectsFlag=1; % if set to 1, only subjects for which delta power is more than threshold for all frequencies are chosen
useMedianFlag=1;

folderLORETA = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\Kanishka_SourceLocalizationProject\data\sLORETA_Thres10'; % Folder where the output of LORETA is saved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mandatory fixed options
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
stRange = [0.25 0.75];

% Choose one of these options
refType = 'unipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjects = getGoodSubjectsProjectwise(projectName,1);
uniqueSubjectNames0 = getGoodFileNamesForSubjects(goodSubjects{1});

%%%%%%%%%%%%%% Find indices for which the correct capType was used %%%%%%%%
capTypeToUse = 'actiCap64';
goodIndices = [];
for i=1:length(uniqueSubjectNames0)
    [expDates,~,capType,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,uniqueSubjectNames0{i},protocolType);
    if usableDataFlag && ~isempty(expDates) && strcmp(capType{1},capTypeToUse)
        goodIndices = cat(2,goodIndices,i);
    end
end
disp([num2str(length(goodIndices)) ' subjects with correct capType chosen for further analysis']);
uniqueSubjectNames = uniqueSubjectNames0(goodIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%% Load Power Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma1Range = [20 34]; gamma2Range = [36 66]; alphaRange = [8 12];
spatialFrequenciesToRemove=[];
dataForDisplay = combineAnalyzedData(pwd,uniqueSubjectNames,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,spatialFrequenciesToRemove,0);

% Here data is saved in the order [SG, FG, alpha]. Change to [alpha SG FG];
newOrderList = [3 1 2];
dataForDisplay.rangeNames = dataForDisplay.rangeNames(newOrderList);
dataForDisplay.powerDBAllSubjects = dataForDisplay.powerDBAllSubjects(:,newOrderList);
dataForDisplay.powerDBTopoAllSubjects = dataForDisplay.powerDBTopoAllSubjects(:,newOrderList,:);

%%%%%%%%%%%%%%%%%%%%%%%% Find Useful Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);
healthyPos = strcmp(cdrList,'HV');
ageGroup1Pos = (ageList<65) & healthyPos; strList{1} = 'mid';
ageGroup2Pos = (ageList>=65) & healthyPos; strList{2} = 'old';
    
numFreqRanges = length(dataForDisplay.rangeNames);
dataDeltaPSD = 10*(dataForDisplay.logSTPowerVsFreqAllSubjects - dataForDisplay.logBLPowerVsFreqAllSubjects);

% Select common good subjects if needed
goodSubjectPosAll = zeros(size(dataForDisplay.powerDBAllSubjects));
for i=1:numFreqRanges
    if strcmp(dataForDisplay.rangeNames{i},'Alpha')
        goodSubjectPosAll(:,i) = dataForDisplay.powerDBAllSubjects(:,i)<-thres;
    else
        goodSubjectPosAll(:,i) = dataForDisplay.powerDBAllSubjects(:,i)>thres;
    end
end
if useCommonSubjectsFlag 
    goodSubjectPosCommon = all(goodSubjectPosAll,2);
end

% Generate plots with 9 columns
% 1) deltaPSD, 2) topoBL, 3) topoST, 
% 4-5) sourceBL & sourceST - change in power
% 6-7) sourceBL & sourceST - tStat
% 8) percentSignificant - barplot
% 9) percentSignificant vs distance from peak

hPlots = getPlotHandles(numFreqRanges,9,[0.05 0.05 0.9 0.9],0.02,0.05);

for i=1:numFreqRanges
    if useCommonSubjectsFlag
        goodSubjectPos = goodSubjectPosCommon;
    else
        goodSubjectPos = goodSubjectPosAll(:,i); %#ok<*UNRCH>
    end
    
    clear goodPos subjectNameListFinal
    gp1 = ageGroup1Pos & goodSubjectPos';
    gp2 = ageGroup2Pos & goodSubjectPos';

    subjectNameListFinal{1} = uniqueSubjectNames(gp1);
    subjectNameListFinal{2} = uniqueSubjectNames(gp2);
    deltaPSD{1} = dataDeltaPSD(gp1,:);
    deltaPSD{2} = dataDeltaPSD(gp2,:);
    
    if useMedianFlag
        topoData{1} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp1),3));
        topoData{2} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp2),3));
    else
        topoData{1} = squeeze(nanmean(dataForDisplay.powerDBTopoAllSubjects(:,i,gp1),3));
        topoData{2} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp2),3));
    end
    
    [allDataBL,allDataST,alltStats,allpVals] = getLORETAData(subjectNameListFinal,strList,folderLORETA);
    for j=1:2
        sourceData(j).BL = squeeze(allDataBL{j}(:,i,:)); %#ok<*SAGROW>
        sourceData(j).ST = squeeze(allDataST{j}(:,i,:));
        sourceData(j).tStats = squeeze(alltStats{j}(:,i,:));
        sourceData(j).pVals = squeeze(allpVals{j}(:,i,:));
    end
    displayData(hPlots(i,:),subjectNameListFinal,strList,deltaPSD,dataForDisplay.freqVals,topoData,sourceData,dataForDisplay.rangeNames{i},refType,useMedianFlag);
end