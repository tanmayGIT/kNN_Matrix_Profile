
% The objective of this function is to show the 1NN v/s 2NN and 1NN v/s 4NN
% matrix profile difference.

% Please run the code until 58th line to see the matching for self-join
% cases on "Protien Data"


% Tanmoy Mondal, Reza Akbarinia,, and Florent Masseglia, "Matrix Profile 
% Based kNN Search over Large Time Series," submitted to: 
% "Elsevier Pattern Recognition Journal", 2020.
% https://sites.google.com/view/knnmatrixprofile/home


clear
clc
close all

subSeqLen = 40;  % The length of the sub-sequence
whichDimToConsider = 1;
kNN_Uwant = 10;   % How many kNN do you want, just tell me


load ('ProteinData.mat');
load ('ProteinLabel.mat');

queryIndex = 2781; % 1165; % randsample(size(realData,1),1);

queryComplete = realData(queryIndex, 1:end);
queryLabel = realDataLabel{queryIndex, 1};
    

% ##############  This is to perfrom kNN matches of all the sub-sequences of query time series. The query time series
% is included in the concatenated target time series. This can also be imagined as self join but applied on only certain
% portion of concatenated time series    ##############

keepTargetData = zeros((size(realData,1)), length(queryComplete));
keepTargetLabels = cell((size(realData,1)), 1);

keepTargetData(1:queryIndex-1,:) = realData(1:queryIndex-1, 1:end);
keepTargetData(queryIndex:end,:) = realData(queryIndex:end, 1:end);

for ik = 1:1:(queryIndex-1)
    keepTargetLabels{ik,1} = realDataLabel{ik, 1};
end
for ik = queryIndex:1:size(realData,1)
    keepTargetLabels{ik,1} = realDataLabel{ik, 1};
end

for iLoop = 1:1:30
    queryIndex = 1461; % kNNSameLabels(iLoop);

    queryComplete = realData(queryIndex, 1:end);
    queryLabel = realDataLabel{queryIndex, 1};
    
    Match_Subseq_N_Query_Par_Journal_1_1(keepTargetData, queryComplete, ...
        subSeqLen, kNN_Uwant, whichDimToConsider, keepTargetLabels, queryLabel, queryIndex);
end






% ##############  This is to perfrom kNN matches of all the sub-sequences of query time series. Hence the query time series
% is excluded from the concatenated target time series. This can also be imagined as independent join     ##############

keepTargetData = zeros((size(realData,1)-1), length(queryComplete));
keepTargetLabels = cell((size(realData,1)-1), 1);

if(queryIndex == 1)
    keepTargetData(:,:) = realData(2:end, 1:end-1);
    for ik = 2:1:size(realData,1)
        keepTargetLabels{ik,1} = queryLabel{ik, 1};
    end
elseif(queryIndex == size(realData,1))
    keepTargetData(:,:) = realData(1:end-1, 1:end-1);
    for ik = 1:1:(size(realData,1)-1)
        keepTargetLabels{ik,1} = queryLabel{ik, 1};
    end
else
    keepTargetData(1:queryIndex-1,:) = realData(1:queryIndex-1, 1:end);
    keepTargetData(queryIndex:end,:) = realData(queryIndex+1:end, 1:end);
    
    for ik = 1:1:(queryIndex-1)
        keepTargetLabels{ik,1} = realDataLabel{ik, 1};
    end
    for ik = queryIndex+1:1:size(realData,1)
        keepTargetLabels{ik-1,1} = realDataLabel{ik, 1};
    end
end

Match_Subseq_N_Query_Par_Journal_1(keepTargetData, queryComplete, ...
                subSeqLen, kNN_Uwant, whichDimToConsider);   % working for the 2nd column here

