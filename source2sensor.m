function sensorData = source2sensor(sourceData,leadfield,goodSourcePos)

if ~exist('goodSourcePos','var');         goodSourcePos = [];           end

if isempty(goodSourcePos)
    goodSourcePos = find(sourceData.inside); % Take all dipoles
end

numTrials = sourceData.df;
sensorData = cell(1,numTrials);

for i=1:numTrials
    mom = sourceData.trial(i).mom;

    nonEmptySourcePos = find(~cellfun(@isempty, mom));
    goodSourcePosThisTrial = intersect(goodSourcePos,nonEmptySourcePos);

    tmpData = [];
    for j=1:length(goodSourcePosThisTrial)
        lf = leadfield.leadfield{goodSourcePosThisTrial(j)};
        tmpMom = mom{goodSourcePosThisTrial(j)};
        if isempty(tmpData)
            tmpData = lf * tmpMom;
        else
            tmpData = tmpData + (lf * tmpMom);
        end
    end
    sensorData{i} = tmpData;
end
end