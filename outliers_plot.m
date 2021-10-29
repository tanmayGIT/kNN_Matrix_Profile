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


putSameLabelIndxs = zeros(1,1);
labelCnt = 1;
for iLabel = 1:1:size(keepAllData,1)
    getLabel = keepAllData(iLabel,end);
    if( getLabel == queryLabel )
        putSameLabelIndxs(labelCnt, 1) = iLabel; % get all the time series indexes whose labels are same as the label of 307th time series
        labelCnt = labelCnt +1;
    end
end

YSameLabels = randsample(labelCnt-1, 30); % generate 30 random numbers in the range of 0 to (labelCnt-1)
kNNSameLabels = putSameLabelIndxs(YSameLabels); % get the time series indexes which are stored at these indexes i.e. YSameLabels




% ##############  This is to fabricate the data to show the multiple occurance of the outliers ##############

colorString = {'y','m','r','b','k','c','g',[.5 .6 .7],[.8 .2 .6], [0.72 0.52 0.04], [0.6 0.8 0.2], [0.619 1 0.4 ],...
    [0.9 0.9 0.9], [0.4 0.4 0.4], [0.8 0.8 0.8], [0.6 0.6 0.6], [0.2 0.2 0.2], [0 0.75 0.75], [0 0.5 0],...
    [0.4 0.58 0.9], [0.75 0 0.75], [0.8 0.7 0.6], [0.6 0.5 0.4 ], [0.8 0.6 1 ], [0 1 1], [1 0.6 0.8]};
hFig1 = figure();
getIndxLabel = 92;
dataToPlot = keepAllData(getIndxLabel, 1:end-1);
[minVal, minIdx] = min(dataToPlot);

rnadomLoc = [9, 102, 201, 291, 365, 446]; % randi([1 500],[1,12]);
for icpy = 1:1:6
    getLoc = rnadomLoc(icpy);
    calSTD = std(dataToPlot(33:42));
    dataToPlot(getLoc:getLoc+9) = dataToPlot(33:42)+(calSTD*icpy);
end

dataToPlot(120:129) = dataToPlot(minIdx-5 : minIdx +4); % copy at the first place
dataToPlot(250:256) = dataToPlot(minIdx-3 : minIdx +3) ; % copy at the second place

% n = 20;
% R = [-4 4];
% z = rand(n,1)*range(R)+min(R);

z = [
    2.5541
    1.9569
    3.8357
    3.3894
    1.2296
    3.4609
    2.6919
    3.3688
    2.3573
    0.6192];

% dataToPlot(161:170) = z(:);

plot(1:length(dataToPlot), dataToPlot(:),'color', colorString{5}, 'LineWidth',1);
hold on

plot(120:129, dataToPlot(120:129),'Color','r', 'LineWidth',1);
plot(minIdx-5:minIdx +4, dataToPlot(minIdx-5:minIdx +4),'Color','r', 'LineWidth',1);
plot(250:259, dataToPlot(250:259),'Color','b', 'LineWidth',1);

for icpy = 1:1:6
    getLoc = rnadomLoc(icpy);
    
    plot(getLoc:getLoc+9, dataToPlot(getLoc:getLoc+9),'Color','m', 'LineWidth',1);
end
plot(432:441, dataToPlot(432:441),'Color','g', 'LineWidth',2.5);

% xticks([0:10:500])
hold off
[sorted_CrossingDisval] = Match_Subseq_N_Query_Par_SelfJoin(dataToPlot', 20, 3);

plotTheQueryOnly(kNNSameLabels, queryLabel, keepAllData, 10);


function plotTheQueryOnly(kNNSameLabels, queryLabel, realData, kNN_Uwant)


colorString = {'y','m','r','b','k','c','g',[.5 .6 .7],[.8 .2 .6], [0.72 0.52 0.04], [0.6 0.8 0.2], [0.619 1 0.4 ],...
    [0.9 0.9 0.9], [0.4 0.4 0.4], [0.8 0.8 0.8], [0.6 0.6 0.6], [0.2 0.2 0.2], [0 0.75 0.75], [0 0.5 0],...
    [0.4 0.58 0.9], [0.75 0 0.75], [0.8 0.7 0.6], [0.6 0.5 0.4 ], [0.8 0.6 1 ], [0 1 1], [1 0.6 0.8]};
hFig1 = figure();

hold on;
for iOut = 1:1:kNN_Uwant
    getIndxLabel = kNNSameLabels(iOut);
    if( realData(getIndxLabel, end) == queryLabel )  % the same label or not
        plot(1:length(realData(getIndxLabel, 1:end-1)), realData(getIndxLabel, 1:end-1),'color', colorString{iOut}, 'LineWidth',1);
    else
        error('The labels should be same');
    end
end
hold off;
% hleg1 = legend(plt1,'plot');
% set(hleg1,'Location','Best')
% set(hleg1,'FontSize',12)

end
