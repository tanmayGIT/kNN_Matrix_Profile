
% The objective of this function is to show the 1NN v/s 2NN and 1NN v/s 4NN
% matrix profile difference.

% Please run the code until 52th line to see the matching for self-join
% cases on "Sheep Data"


% Tanmoy Mondal, Reza Akbarinia,, and Florent Masseglia, "Matrix Profile 
% Based kNN Search over Large Time Series," submitted to: 
% "Elsevier Pattern Recognition Journal", 2020.
% https://sites.google.com/view/knnmatrixprofile/home


clear
clc
close all

subSeqLen = 50;  % The length of the sub-sequence
whichDimToConsider = 1;
kNN_Uwant = 10;   % How many kNN do you want, just tell me

% Loading the sheep data
load('mat_files/SheepData.mat')
% queryIndex = 256; % for label 5
queryIndex = 307; % for label 6
queryLabel = keepAllData(queryIndex, end);


% ##############  This is to perfrom kNN matches of all the sub-sequences of query time series. The query time series
% is included in the concatenated target time series. This can also be imagined as self join but applied on only certain
% portion of concatenated time series    ##############

%for iLoop = 1:1:30
queryIndex = 7882; % kNNSameLabels(iLoop);
queryComplete = keepAllData(queryIndex, 1:end-1);

noOfZerosInQuery = sum(queryComplete(:)==0);
if(noOfZerosInQuery > (length(queryComplete) /2)) % if the number of zeros is more than half the length of query
    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
end

keepTargetData = zeros((size(keepAllData,1)), length(queryComplete));
keepTargetLabels = zeros((size(keepAllData,1)), 1);


keepTargetData(1:queryIndex,:) = keepAllData(1:queryIndex, 1:end-1);
keepTargetData(queryIndex+1:end,:) = keepAllData(queryIndex+1:end, 1:end-1);

keepTargetLabels(1:queryIndex,1) = keepAllData(1:queryIndex, end);
keepTargetLabels(queryIndex+1:end,1) = keepAllData(queryIndex+1:end, end);

Match_Subseq_N_Query_Par_Journal_1_1(keepTargetData, queryComplete, ...
    subSeqLen, kNN_Uwant, whichDimToConsider, keepTargetLabels, queryLabel, queryIndex);
%  #############################                    #####################################








% ##############  This is to perfrom kNN matches of all the sub-sequences of query time series. Hence the query time series
% is excluded from the concatenated target time series. This can also be imagined as independent join     ##############

%for iLoop = 1:1:30
queryIndex = 7882; % kNNSameLabels(iLoop);
queryComplete = keepAllData(queryIndex, 1:end-1);

noOfZerosInQuery = sum(queryComplete(:)==0);
if(noOfZerosInQuery > (length(queryComplete) /2)) % if the number of zeros is more than half the length of query
    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
end

keepTargetData = zeros((size(keepAllData,1)-1), length(queryComplete));
keepTargetLabels = zeros((size(keepAllData,1)-1), 1);

if(queryIndex == 1)
    keepTargetData(:,:) = keepAllData(2:end, 1:end-1);
    keepTargetLabels(:,1) = keepAllData(2:end, end);
elseif(queryIndex == size(keepAllData,1))
    keepTargetData(:,:) = keepAllData(1:end-1, 1:end-1);
    keepTargetLabels(:,1) = keepAllData(1:end-1, end);
else
    keepTargetData(1:queryIndex-1,:) = keepAllData(1:queryIndex-1, 1:end-1);
    keepTargetData(queryIndex:end,:) = keepAllData(queryIndex+1:end, 1:end-1);
    
    keepTargetLabels(1:queryIndex-1,1) = keepAllData(1:queryIndex-1, end);
    keepTargetLabels(queryIndex:end,1) = keepAllData(queryIndex+1:end, end);
end

Match_Subseq_N_Query_Par_Journal_1_1(keepTargetData, queryComplete, ...
    subSeqLen, kNN_Uwant, whichDimToConsider, keepTargetLabels, queryLabel, queryIndex);   % working for the 2nd column here


