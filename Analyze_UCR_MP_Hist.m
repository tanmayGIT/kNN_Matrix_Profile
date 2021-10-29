clear;
clc;
close all;

load('UCR_Stat_MP.mat')
dataFolderPath = '/home/tmondal/Videos/Dataset/Time_Series/UCRArchive_2018/';
ptCnt = 1;

keepDatasetInfoThrsh = cell(1,1);
goodCellCnt = 1;

for ii = 1:1:length(keepInfoAllData)
    if(~isempty(keepInfoAllData{ii}))
        
        getCell = keepInfoAllData{ii};
        datasetName = getCell.datasetName;
        
        datasetPath_Train = strcat(dataFolderPath, datasetName, '/', datasetName, '_TRAIN.tsv');
        datasetPath_Test = strcat(dataFolderPath, datasetName, '/', datasetName, '_TEST.tsv');
        
        TRAIN = load(datasetPath_Train); % Only these two lines need to be changed to test a different dataset %
        TEST  = load(datasetPath_Test); % Only these two lines need to be changed to test a different dataset %
               
           
        getMax2NN = getCell.getMax2NN;
        getMax4NN = getCell.getMax4NN;
        keepAllNormEachQuery = zeros(1,5);
        
        for jj = 1:1:length(getCell.keepAllQueryInfo)
            innerCell = getCell.keepAllQueryInfo{jj};
            if(~isempty(innerCell))
                getDiff1NN2NN = innerCell.SortedVal_1NN_2NN;
                SortedIndx_1NN_2NN = innerCell.SortedIndx_1NN_2NN;
                
                getDiff1NN2NNNorm = getDiff1NN2NN/getMax2NN;
                
                getDiff1NN4NN = innerCell.SortedVal_1NN_4NN;
                SortedIndx_1NN_4NN = innerCell.SortedIndx_1NN_4NN;
                
                getDiff1NN4NNNorm = getDiff1NN4NN/getMax4NN;
                
                
                keepAllNormEachQuery(ptCnt:ptCnt+length(getDiff1NN2NNNorm)-1,1) = getDiff1NN2NNNorm(:);
                keepAllNormEachQuery(ptCnt:ptCnt+length(getDiff1NN2NNNorm)-1,2) = SortedIndx_1NN_2NN(:);
                
                keepAllNormEachQuery(ptCnt:ptCnt+length(getDiff1NN4NNNorm)-1,3) = getDiff1NN4NNNorm(:);
                keepAllNormEachQuery(ptCnt:ptCnt+length(getDiff1NN4NNNorm)-1,4) = SortedIndx_1NN_4NN(:);
                
                keepAllNormEachQuery(ptCnt:ptCnt+length(getDiff1NN4NNNorm)-1,5) = ii; % dataset name
                
                ptCnt = ptCnt+length(getDiff1NN4NNNorm);
                
            end
        end
        % Now check that how many are passing the threshold of 2NN and 4NN
        thresh_2NN_1 = (getMax2NN*2)/100;
        thresh_2NN_2 = (getMax2NN*5)/100;
        thresh_2NN_3 = (getMax2NN*10)/100;
        
        thresh_4NN_1 = (getMax4NN*2)/100;
        thresh_4NN_2 = (getMax4NN*5)/100;
        thresh_4NN_3 = (getMax4NN*10)/100;
        
        valPass2nnThresh_1 = length(find(keepAllNormEachQuery(:,1) > thresh_2NN_1));
        valPass2nnThresh_2 = length(find(keepAllNormEachQuery(:,1) > thresh_2NN_2));
        valPass2nnThresh_3 = length(find(keepAllNormEachQuery(:,1) > thresh_2NN_3));
        
        valPass4nnThresh_1 = length(find(keepAllNormEachQuery(:,3) > thresh_4NN_1));
        valPass4nnThresh_2 = length(find(keepAllNormEachQuery(:,3) > thresh_4NN_2));
        valPass4nnThresh_3 = length(find(keepAllNormEachQuery(:,3) > thresh_4NN_3));
        
        keepDatasetInfoThrsh{goodCellCnt, 1}.DatasetName = datasetName;
        keepDatasetInfoThrsh{goodCellCnt, 1}.valPass2nnThresh_1 = valPass2nnThresh_1;
        keepDatasetInfoThrsh{goodCellCnt, 1}.valPass2nnThresh_2 = valPass2nnThresh_2;
        keepDatasetInfoThrsh{goodCellCnt, 1}.valPass2nnThresh_3 = valPass2nnThresh_3;
        
        keepDatasetInfoThrsh{goodCellCnt, 1}.valPass4nnThresh_1 = valPass4nnThresh_1;
        keepDatasetInfoThrsh{goodCellCnt, 1}.valPass4nnThresh_2 = valPass4nnThresh_2;
        keepDatasetInfoThrsh{goodCellCnt, 1}.valPass4nnThresh_3 = valPass4nnThresh_3;
        
        fprintf('%d & %s & %d & %d & %d & %d & %d & %d & %d & %d & %d \\\\\n', goodCellCnt, datasetName, size(TRAIN,1), ...
        size(TEST,1), size(TEST,2), valPass2nnThresh_1, valPass2nnThresh_2, valPass2nnThresh_3, valPass4nnThresh_1, valPass4nnThresh_2,valPass4nnThresh_3); 
        
        goodCellCnt = goodCellCnt+1;
    end
end

% Perform the bar plot
all_2NN_thresh = zeros(length(keepDatasetInfoThrsh), 3);
for jj = 1:1:length(keepDatasetInfoThrsh)
    getCell = keepDatasetInfoThrsh{jj, 1};
    
    all_2NN_thresh(jj, 1) = getCell.valPass4nnThresh_1;
    all_2NN_thresh(jj, 2) = getCell.valPass4nnThresh_2;
    all_2NN_thresh(jj, 3) = getCell.valPass4nnThresh_3;
end
x = 1:length(keepDatasetInfoThrsh);
xStr = string(x);
figure
bar(x,all_2NN_thresh(:,1), 0.6,'FaceColor',[0.6350 0.0780 0.1840])
hold on
bar(x, all_2NN_thresh(:,2), 0.4, 'FaceColor',[0 0.7 0.7])
hold on
bar(x, all_2NN_thresh(:,3), 0.2, 'FaceColor',	[1 0 1])


xgoLabel = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';...
    '17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';...
    '31';'32';'33';'34';'35';'36';'37';'38';'39';'40';'41';'42';'43';'44';...
    '45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';...
    '59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71';'72';...
    '73';'74';'75';'76';'77';'78';'79';'80';'81';'82';'83';'84';'85';'86';...
    '87'};
set(gca, 'XTick', 1:length(xgoLabel),'XTickLabel',xgoLabel, 'XTickLabelRotation', 50);
grid on
ylabel('No. of extra outliers or inliers')
legend({'Threshold : 2% ','Threshold : 5%', 'Threshold : 10%'},'Location','northwest')
disp('hello')

