clear;
clc;
close all;

load('UCR_Stat_MP.mat')
keepAllNorm = zeros(1,5);

ptCnt = 1;
for ii = 1:1:length(keepInfoAllData)
    if(~isempty(keepInfoAllData{ii}))
        
        getCell = keepInfoAllData{ii};
        allQueryIndex = getCell.AllQueryIndex;
        datasetName = getCell.datasetName; % the dataset name can also be mentioned by the index ii as we know that at the iith index
        % which dataset is saved
        
        getMax2NN = getCell.getMax2NN;
        getMax4NN = getCell.getMax4NN;
        
        for jj = 1:1:length(getCell.keepAllQueryInfo)
            innerCell = getCell.keepAllQueryInfo{jj};
            if(~isempty(innerCell))
                getDiff1NN2NN = innerCell.SortedVal_1NN_2NN;
                SortedIndx_1NN_2NN = innerCell.SortedIndx_1NN_2NN;
                
                getDiff1NN2NNNorm = getDiff1NN2NN/getMax2NN;
                
                getDiff1NN4NN = innerCell.SortedVal_1NN_4NN;
                SortedIndx_1NN_4NN = innerCell.SortedIndx_1NN_4NN;
                
                getDiff1NN4NNNorm = getDiff1NN4NN/getMax4NN;
                
                
                keepAllNorm(ptCnt:ptCnt+length(getDiff1NN2NNNorm)-1,1) = getDiff1NN2NNNorm(:);
                keepAllNorm(ptCnt:ptCnt+length(getDiff1NN2NNNorm)-1,2) = SortedIndx_1NN_2NN(:);
                
                keepAllNorm(ptCnt:ptCnt+length(getDiff1NN4NNNorm)-1,3) = getDiff1NN4NNNorm(:);
                keepAllNorm(ptCnt:ptCnt+length(getDiff1NN4NNNorm)-1,4) = SortedIndx_1NN_4NN(:);
                
                keepAllNorm(ptCnt:ptCnt+length(getDiff1NN4NNNorm)-1,5) = ii; % dataset name
                
                ptCnt = ptCnt+length(getDiff1NN4NNNorm);
                
            end
        end
    end
end
[sortedVal_1NN_2NN, sortedIndx_1NN_2NN]= sort(keepAllNorm(:,1), 'descend');
modifiedKeepAllNorm_1NN_2NN = zeros(size(keepAllNorm,1),3);
modifiedKeepAllNorm_1NN_2NN(:,1) = sortedVal_1NN_2NN(:);
modifiedKeepAllNorm_1NN_2NN(:,2) = keepAllNorm(sortedIndx_1NN_2NN,2);
modifiedKeepAllNorm_1NN_2NN(:,3) = keepAllNorm(sortedIndx_1NN_2NN,5);

[sortedVal_1NN_4NN, sortedIndx_1NN_4NN]= sort(keepAllNorm(:,3), 'descend');
modifiedKeepAllNorm_1NN_4NN = zeros(size(keepAllNorm,1),3);
modifiedKeepAllNorm_1NN_4NN(:,1) = sortedVal_1NN_4NN(:);
modifiedKeepAllNorm_1NN_4NN(:,2) = keepAllNorm(sortedIndx_1NN_4NN,4);
modifiedKeepAllNorm_1NN_4NN(:,3) = keepAllNorm(sortedIndx_1NN_4NN,5);

hist_1NN_2NN = zeros(123,2);
hist_1NN_4NN = zeros(123,2);

for ij = 1:1:size(modifiedKeepAllNorm_1NN_2NN,1)
    getVal1 = modifiedKeepAllNorm_1NN_2NN(ij,3);
    getVal2 = modifiedKeepAllNorm_1NN_4NN(ij,3);
    
    hist_1NN_2NN(getVal1,1) = hist_1NN_2NN(getVal1,1)+1;
    hist_1NN_4NN(getVal2,1) = hist_1NN_4NN(getVal2,1)+1;
    
    % get the avarage indexes of these values
    hist_1NN_2NN(getVal1,2) = hist_1NN_2NN(getVal1,2) + ij;
    hist_1NN_4NN(getVal2,2) = hist_1NN_4NN(getVal2,2) + ij;
end

hist_1NN_2NN(:,2) = round(hist_1NN_2NN(:,2)./hist_1NN_2NN(:,1));
hist_1NN_4NN(:,2) = round(hist_1NN_4NN(:,2)./hist_1NN_4NN(:,1));

% now plot the MP
for ik = 1:1:size(hist_1NN_2NN,1)
    if(hist_1NN_2NN(ik,1) > 0)
        getDatasetindex = ik;
        RecreateMatrixProfile(getDatasetindex, keepInfoAllData);
    end
end


function RecreateMatrixProfile(datasetIndex, keepInfoAllData)
dataFolderPath = '/home/tmondal/Videos/Dataset/Time_Series/UCRArchive_2018/';

subSeqLen = 64;  % The length of the sub-sequence

getCell = keepInfoAllData{datasetIndex};
allQueryIndex = getCell.AllQueryIndex;
datasetName = getCell.datasetName;
disp(['The dataset name is ', datasetName ])
yourFolder = strcat('./UCR_Matrix_Profile_Results/', datasetName, '/');
if ~exist(yourFolder, 'dir')
    mkdir(yourFolder)
    
    
    datasetPath_Train = strcat(dataFolderPath, datasetName, '/', datasetName, '_TRAIN.tsv');
    datasetPath_Test = strcat(dataFolderPath, datasetName, '/', datasetName, '_TEST.tsv');
    
    TRAIN = load(datasetPath_Train); % Only these two lines need to be changed to test a different dataset %
    TEST  = load(datasetPath_Test); % Only these two lines need to be changed to test a different dataset %
    
    Load_UCR_Data(TRAIN, TEST, subSeqLen, allQueryIndex, yourFolder);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following function is taken from : https://www.cs.ucr.edu/%7Eeamonn/time_series_data_2018/BriefingDocument2018.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Load_UCR_Data(TRAIN, TEST, subSeqLen, queryIndex, yourFolder)

whichDimToConsider = 1;
kNN_Uwant = 10;   % How many kNN do you want, just tell me

TRAIN_class_labels = TRAIN(:,1);    % Pull out the class labels.
TRAIN(:,1) = [];                    % Remove class labels from training set.
TEST_class_labels = TEST(:,1);      % Pull out the class labels.
TEST(:,1) = [];                     % Remove class labels from testing set.

disp(['The dataset you tested has ', int2str(length(unique(TRAIN_class_labels))), ' classes'])
disp(['The training set is of size ', int2str(size(TRAIN,1)),', and the test set is of size ',int2str(size(TEST,1)),'.'])
disp(['The time series are of length ', int2str(size(TRAIN,2))])

% Joining/concatenating the TRAIN and TEST data together
TRAIN_TEST_Data = zeros(size(TRAIN, 1)+size(TEST, 1), size(TRAIN, 2));
TRAIN_TEST_Label = zeros(size(TRAIN_class_labels, 1)+size(TEST_class_labels, 1), size(TRAIN_class_labels, 2));

TRAIN_TEST_Data(1:size(TRAIN, 1), :) = TRAIN(:,:);
TRAIN_TEST_Data( (size(TRAIN, 1)+1) : end, :) = TEST(:,:);

% Putting the labels together also
TRAIN_TEST_Label(1:size(TRAIN_class_labels, 1), :) = TRAIN_class_labels(:,:);
TRAIN_TEST_Label( (size(TRAIN_class_labels, 1)+1) : end, :) = TEST_class_labels(:,:);

for ii = 1:1:length(queryIndex)
    queryComplete = TRAIN_TEST_Data(queryIndex(ii), 1:end);
    queryLabel = TRAIN_TEST_Label(queryIndex(ii), 1);
    
    noOfZerosInQuery = sum(queryComplete(:)== 0);
    if(noOfZerosInQuery > (length(queryComplete) /2)) % if the number of zeros is more than half the length of query
        error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
    end
    
    [keepTargetData, ~] = OrganizeTargetData(TRAIN_TEST_Data, TRAIN_TEST_Label, queryComplete, queryIndex(ii), queryLabel);
    badFlag = CheckQuery(queryComplete);
    
    if(~badFlag) % if this flag remains false
        [keepFileDistEntryOnly, ~] = Match_Subseq_N_Query_Par_UCR(keepTargetData, queryComplete, subSeqLen, kNN_Uwant, whichDimToConsider, queryIndex(ii));   % working for the 2nd column here
        
        plotTheGraphOnly2NN(keepFileDistEntryOnly, 1, 2, 1, 2, queryIndex(ii), yourFolder); % comparison between 1NN and 2NN
        plotTheGraphOnly2NN(keepFileDistEntryOnly, 1, 4, 1, 5, queryIndex(ii), yourFolder); % comparison between 1NN and 4NN
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


function plotTheGraphOnly2NN(keepFileDistEntryOnly, firstNN, secondNN, firstColor, secondColor, queryIndex, yourFolder)

colorString = {'b','m','r','y','k','c','g',[.5 .6 .7],[.8 .2 .6], [0.72 0.52 0.04], [0.6 0.8 0.2], [0.619 1 0.4 ],...
    [0.9 0.9 0.9], [0.4 0.4 0.4], [0.8 0.8 0.8], [0.6 0.6 0.6], [0.2 0.2 0.2], [0 0.75 0.75], [0 0.5 0],...
    [0.4 0.58 0.9], [0.75 0 0.75], [0.8 0.7 0.6], [0.6 0.5 0.4 ], [0.8 0.6 1 ], [0 1 1], [1 0.6 0.8]};
hFig1 = figure();
set(hFig1,'WindowState','fullscreen')
hold on;
plot(1:length(keepFileDistEntryOnly(firstNN, :)), keepFileDistEntryOnly(firstNN, :),'color', colorString{firstColor}, 'LineWidth',1);
plot(1:length(keepFileDistEntryOnly(secondNN, :)), keepFileDistEntryOnly(secondNN, :),'color', colorString{secondColor}, 'LineWidth',1);

imageNam = strcat(yourFolder,num2str(queryIndex), '_', num2str(firstNN), '_', num2str(secondNN),'.png');
saveas(hFig1,imageNam);
hold off;
close(hFig1);

end
