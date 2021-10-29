clear;
clc;
close all;

dataFolderPath = '/home/tmondal/UCR_Matrix_Profile_Exp/UCRArchive_2018/';
subFolders = GetSubDirsFirstLevelOnly(dataFolderPath);

getMax2NN = 0;
getMax4NN = 0;
subSeqLen = 64;  % The length of the sub-sequence
keepInfoAllData = cell(length(subFolders), 1);


for iDir = 1:1:length(subFolders)
    datasetName = subFolders{1, iDir}; % 'OSULeaf';
    
    datasetPath_Train = strcat(dataFolderPath, datasetName, '/', datasetName, '_TRAIN.tsv');
    datasetPath_Test = strcat(dataFolderPath, datasetName, '/', datasetName, '_TEST.tsv');
    disp(['The name of the dataset is : ', datasetName ])    
	TRAIN = load(datasetPath_Train); % Only these two lines need to be changed to test a different dataset %
	TEST  = load(datasetPath_Test); % Only these two lines need to be changed to test a different dataset %

	if (subSeqLen < (size(TRAIN,2) / 2))
	    [keepAllQueryInfo, queryIndex, getMax2NN, getMax4NN, noGoodQueryFlag] = Load_UCR_Data(TRAIN, TEST, getMax2NN, getMax4NN, subSeqLen);
    	if(noGoodQueryFlag)
	    	keepInfoAllData{iDir,1}.keepAllQueryInfo = keepAllQueryInfo;
	    	keepInfoAllData{iDir,1}.AllQueryIndex = queryIndex;
	    	keepInfoAllData{iDir,1}.getMax2NN = getMax2NN;
	    	keepInfoAllData{iDir,1}.getMax4NN = getMax4NN;
            keepInfoAllData{iDir,1}.datasetName = datasetName;
		end		
	end
end
save('UCR_Stat_MP.mat', 'keepInfoAllData');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following function is taken from : https://www.cs.ucr.edu/%7Eeamonn/time_series_data_2018/BriefingDocument2018.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [keepAllQueryInfo, queryIndex, getMax2NN, getMax4NN, noGoodQueryFlag] = Load_UCR_Data(TRAIN, TEST, getMax2NN, getMax4NN, subSeqLen)

whichDimToConsider = 1;
kNN_Uwant = 10;   % How many kNN do you want, just tell me

TRAIN_class_labels = TRAIN(:,1);    % Pull out the class labels.
TRAIN(:,1) = [];                    % Remove class labels from training set.
TEST_class_labels = TEST(:,1);      % Pull out the class labels.
TEST(:,1) = [];                     % Remove class labels from testing set.

disp(['The dataset you are working with has ', int2str(length(unique(TRAIN_class_labels))), ' classes'])
disp(['The training set is of size ', int2str(size(TRAIN,1)),', and the test set is of size ',int2str(size(TEST,1)),'.'])
disp(['The time series are of length ', int2str(size(TRAIN,2))])
fprintf('\n\n');

% Joining/concatenating the TRAIN and TEST data together
TRAIN_TEST_Data = zeros(size(TRAIN, 1)+size(TEST, 1), size(TRAIN, 2));
TRAIN_TEST_Label = zeros(size(TRAIN_class_labels, 1)+size(TEST_class_labels, 1), size(TRAIN_class_labels, 2));

TRAIN_TEST_Data(1:size(TRAIN, 1), :) = TRAIN(:,:);
TRAIN_TEST_Data( (size(TRAIN, 1)+1) : end, :) = TEST(:,:);

% Putting the labels together also
TRAIN_TEST_Label(1:size(TRAIN_class_labels, 1), :) = TRAIN_class_labels(:,:);
TRAIN_TEST_Label( (size(TRAIN_class_labels, 1)+1) : end, :) = TEST_class_labels(:,:);

nEntry = size(TRAIN_TEST_Data,1);
queryIndex = randsample(nEntry, 30); % randsample(size(TRAIN_TEST_Data,1),1); % randomly choose an index from train data

keepAllQueryInfo = cell(1,1);
noGoodQueryFlag = false;
cntEnter = 1;
for ii = 1:1:length(queryIndex)
    queryComplete = TRAIN_TEST_Data(queryIndex(ii), 1:end);
    queryLabel = TRAIN_TEST_Label(queryIndex(ii), 1);
    
    noOfZerosInQuery = sum(queryComplete(:)== 0);
    if(noOfZerosInQuery > (length(queryComplete) /2)) % if the number of zeros is more than half the length of query
        error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
    end
    
    [keepTargetData, keepTargetLabels] = OrganizeTargetData(TRAIN_TEST_Data, TRAIN_TEST_Label, queryComplete, queryIndex(ii), queryLabel);

    queryComplete = CheckData(queryComplete);
    badFlag = CheckQuery(queryComplete);
    
    if(~badFlag) % if this flag remains false
        [keepFileDistEntryOnly, goodQuery] = Match_Subseq_N_Query_Par_UCR(keepTargetData, queryComplete, subSeqLen, kNN_Uwant, whichDimToConsider, queryIndex(ii));   % working for the 2nd column here
        
        if(goodQuery) % if it is true 
            % Here in this part of the code, I am checking whether there are some
            % interesting thing in 1NN, 2NN and 4NN matrix profile; If there is a
            % high difference of values between 1NN v/s 2NN and 1NN v/s 4NN then
            % there is a chance that these points are interesting

            diff1NN2NN = abs(keepFileDistEntryOnly(2,:) - keepFileDistEntryOnly(1,:)); % for strong motifs
            [sortedVal1NN2NN, sortedIndx1NN2NN] = sort(diff1NN2NN, 'descend'); % sort these values

            diff1NN4NN = abs(keepFileDistEntryOnly(4,:) - keepFileDistEntryOnly(1,:)); % for strong motifs
            [sortedVal1NN4NN, sortedIndx1NN4NN] = sort(diff1NN4NN, 'descend'); % sort these values

            [~,indxCol1NN_2NN,~] = find(keepFileDistEntryOnly(2,:)  < keepFileDistEntryOnly(1,:) ); % see when 2NN MP is less than 1NN MP; genrally they should be opposite
            [~,indxCol1NN_4NN,~] = find(keepFileDistEntryOnly(4,:)  < keepFileDistEntryOnly(1,:) ); % see when 4NN MP is less than 11NN MP; genrally they should be opposite

            keepAllQueryInfo{cntEnter,1}.SortedVal_1NN_2NN = sortedVal1NN2NN(1, 1:end);
            keepAllQueryInfo{cntEnter,1}.SortedIndx_1NN_2NN = sortedIndx1NN2NN(1, 1:end);
            keepAllQueryInfo{cntEnter,1}.Actual_1NN_MP = keepFileDistEntryOnly(1,:);
            keepAllQueryInfo{cntEnter,1}.Actual_2NN_MP = keepFileDistEntryOnly(2,:);

            keepAllQueryInfo{cntEnter,1}.SortedVal_1NN_4NN = sortedVal1NN4NN(1, 1:end);
            keepAllQueryInfo{cntEnter,1}.SortedIndx_1NN_4NN = sortedIndx1NN4NN(1, 1:end);
            keepAllQueryInfo{cntEnter,1}.Actual_4NN_MP = keepFileDistEntryOnly(4,:);

            keepAllQueryInfo{cntEnter,1}.Strange_1NN_2NN = indxCol1NN_2NN;
            keepAllQueryInfo{cntEnter,1}.Strange_1NN_4NN = indxCol1NN_4NN;

            if(max(sortedVal1NN2NN(1, 1:20)) > getMax2NN)
                getMax2NN = max(sortedVal1NN2NN(1, 1:20));
            end
            if(max(sortedVal1NN4NN(1, 1:20)) > getMax4NN)
                getMax4NN = max(sortedVal1NN4NN(1, 1:20));
            end
            
            noGoodQueryFlag = true;
            cntEnter = cntEnter+1;
        end
    end
end
return;
end


function badFlag = CheckQuery(data)
badFlag = false;
nanIndexesTar = find(isnan(data(:) ));
if(~isempty(nanIndexesTar) )
    badFlag = true;
    return
end
infIndexesTar = find(isinf(data(:) ));
if(~isempty(infIndexesTar) )
    badFlag = true;
    return
end
return;
end


function data = CheckData(data)
nanIndexesTar = (isnan(data(:) ));
infIndexesTar = (isinf(data(:) ));
rowNonInf = (~isfinite(data(:) ));

allBadCell = bitor((bitor(nanIndexesTar,infIndexesTar)),rowNonInf);
tAllIndex = 1:numel(data(:));
data(allBadCell) = interp1(tAllIndex(~allBadCell), data(~allBadCell), tAllIndex(allBadCell));

return;
end


function [keepTargetData, keepTargetLabels] = ...
    OrganizeTargetData(realData, realDataLabel, queryComplete, queryIndex, queryLabel)

keepTargetData = zeros((size(realData,1)-1), length(queryComplete));
keepTargetLabels = zeros((size(realDataLabel,1)-1), 1);

if(queryIndex == 1)
    keepTargetData(:,:) = realData(2:end, :);
    keepTargetLabels(:,1) = realDataLabel(2:end, 1);
    
elseif(queryIndex == size(realData,1))
    keepTargetData(:,:) = realData(1:end-1, :);
    keepTargetLabels(:,1) = realDataLabel(1:end-1, 1);
else
    keepTargetData(1:queryIndex-1,:) = realData(1:queryIndex-1, :);
    keepTargetData(queryIndex:end,:) = realData(queryIndex+1:end, :);
    
    keepTargetLabels(1:queryIndex-1,1) = realDataLabel(1:queryIndex-1, 1);
    keepTargetLabels(queryIndex:end,1) = realDataLabel(queryIndex+1:end, 1);
end

return;
end

function [subDirsNames] = GetSubDirsFirstLevelOnly(parentDir)
    % Get a list of all files and folders in this folder.
    files = dir(parentDir);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subDirs = files(dirFlags);
    subDirsNames = cell(1, numel(subDirs) - 2);
    for i=3:numel(subDirs)
        subDirsNames{i-2} = subDirs(i).name;
    end
end
