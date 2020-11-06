
% The following code is to plot the graph of time v/s kNN 
% or time v/s subSequence Length or time v/s time series length or 
% time v/s no. of cores. To plot each graph, you need to commnet and uncommnet 
% following block of code (seperated by several line space)

% For example, the following block of code 55-80 are uncommented. Next, you
% need to comment the lines 55-80 and uncomment the lines 91-110 to see the
% next plot.

% Tanmoy Mondal, Reza Akbarinia,, and Florent Masseglia, "Matrix Profile 
% Based kNN Search over Large Time Series," submitted to: 
% "Elsevier Pattern Recognition Journal", 2020.
% https://sites.google.com/view/knnmatrixprofile/home


clc;
%close all;
clear 


% Self Join Only

textFilePathSubSeq3 = 'Part_3/799randomWalk_kNN_HeapBased_Self_Join.txt';
textFilePathSubSeq4 = 'Part_3/799_seismic_kNN_HeapBased_Self_Join.txt';

textFilePathSubSeq3_1 = 'Part_3/799randomWalk_kNN_SortBased_Self_Join.txt';
textFilePathSubSeq4_1 = 'Part_3/799_seismic_kNN_SortBased_Self_Join.txt';

textFilePathSubSeq3_2 = 'Part_3/799randomWalk_kNN_MaxBased_Self_Join.txt';
textFilePathSubSeq4_2 = 'Part_3/799_seismic_kNN_MaxBased_Self_Join.txt';

keepDataRandomWalk_kNNHeapBased = ReadFile_2(textFilePathSubSeq3);
keepDataSiesmic_kNNHeapBased = ReadFile_2(textFilePathSubSeq4);

keepDataRandomWalk_kNNSortBased = ReadFile_2(textFilePathSubSeq3_1);
keepDataSiesmic_kNNSortBased = ReadFile_2(textFilePathSubSeq4_1);

keepDataRandomWalk_kNNMaxBased = ReadFile_2(textFilePathSubSeq3_2);
keepDataSiesmic_kNNMaxBased = ReadFile_2(textFilePathSubSeq4_2);


figure();
h(1) = subplot(1,1,1);

colorspec = {[0.9 0.9 0.9]; [0.8 0.8 0.8]; [0.6 0.6 0.6]; ...
  [0.4 0.4 0.4]; [0.2 0.2 0.2];[0 0.75 0.75];[0 0.5 0];[0.75 0.75 0];...
  [1 0.50 0.25];[0.75 0 0.75];[0.7 0.7 0.7];[0.8 0.7 0.6];[0.6 0.5 0.4 ]};








q = plot(keepDataRandomWalk_kNNHeapBased(:,1),keepDataRandomWalk_kNNHeapBased(:,2),'r-*', ...
        keepDataRandomWalk_kNNSortBased(:,1),keepDataRandomWalk_kNNSortBased(:,2),'r--d', ...
        keepDataRandomWalk_kNNMaxBased(:,1),keepDataRandomWalk_kNNMaxBased(:,2),'r:s',...
        keepDataSiesmic_kNNHeapBased(:,1),keepDataSiesmic_kNNHeapBased(:,2),'b-*',...
        keepDataSiesmic_kNNSortBased(:,1),keepDataSiesmic_kNNSortBased(:,2),'b--d',...
        keepDataSiesmic_kNNMaxBased(:,1),keepDataSiesmic_kNNMaxBased(:,2),'b:s');
    
hold on;
hleg1 = legend('kNN Heap-Max-Based on RandomWalk Data(50000 * 256)', ...
                'kNN Sorting-Based on RandomWalk Data(50000 * 256)',... 
                'kNN Max-Based on RandomWalk Data(50000 * 256)',...
                'kNN Heap-Max-Based on Seismic Data(50000 * 200)',...
                'kNN Sorting-Based on Seismic Data(50000 * 200)',...
                'kNN Max-Based on Seismic Data(50000 * 200)');
            
set(hleg1,'Location','NorthWest')
set(hleg1,'FontSize',16)
set(gca,'XTick', keepDataRandomWalk_kNNHeapBased(:,1));
grid on;
set(gca,'FontSize',10);
xl = xlabel('k Nearest Neighbors');
yl = ylabel('Time Needed (sec.)');
set(xl,'FontSize',14,'FontWeight','bold','FontName','Courier');
set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
hold off;
disp('see me');









% q = plot(keepDataRandomWalk15_1(:,1),keepDataRandomWalk15_1(:,2),'r-*', ...
%         keepDataRandomWalk15_2(:,1),keepDataRandomWalk15_2(:,2),'r--d', ...
%         keepDataSeismic16_1(:,1),keepDataSeismic16_1(:,2),'b-*',...
%         keepDataSeismic16_2(:,1),keepDataSeismic16_2(:,2),'b--d');
%     
% hold on;
% hleg1 = legend('Proposed Self Join on RandomWalk Data (no. of time series * 256)', 'Yeh et al. Self Join on  RandomWalk Data (no. of time series * 256)', ...
%                 ' Proposed Self Join on Seismic Data (no. of time series * 200)','Yeh et al. Self Join on Seismic Data (no. of time series * 200)');
%             
% set(hleg1,'Location','NorthWest')
% set(hleg1,'FontSize',14)
% set(gca,'XTick', keepDataRandomWalk15_1(:,1));
% grid on;
% set(gca,'FontSize',10);
% xl = xlabel('No. of Time Series');
% yl = ylabel('Time Needed (sec.)');
% set(xl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% hold off;
% disp('see me');














function keepData = ReadFile_1(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    newStr = split(tline,' ');
    subSeqlen = str2double(newStr{5});
    timeNeed = str2double(newStr{11}); 
    
    keepData(lnCnt, 1) = subSeqlen;
    keepData(lnCnt, 2) = timeNeed;
    
    lnCnt = lnCnt +1;
    tline = fgetl(fid);
end
fclose(fid);
return;
end

function keepData = ReadFile_2(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    newStr = split(tline,' ');
    subSeqlen = str2double(newStr{6});
    timeNeed = str2double(newStr{12}); 
    
    keepData(lnCnt, 1) = subSeqlen;
    keepData(lnCnt, 2) = timeNeed;
    
    lnCnt = lnCnt +1;
    tline = fgetl(fid);
end
fclose(fid);
return;
end

function keepData = ReadFile_3(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    newStr = split(tline,' ');
    subSeqlen = str2double(newStr{6});
    timeNeed = str2double(newStr{12}); 
    
    keepData(lnCnt, 1) = subSeqlen;
    keepData(lnCnt, 2) = timeNeed;
    
    lnCnt = lnCnt +1;
    tline = fgetl(fid);
end
fclose(fid);
return;
end

function keepData = ReadFile_4(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    newStr = split(tline,' ');
    subSeqlen = str2double(newStr{6});
    timeNeed = str2double(newStr{12}); 
    
    keepData(lnCnt, 1) = subSeqlen;
    keepData(lnCnt, 2) = timeNeed;
    
    lnCnt = lnCnt +1;
    tline = fgetl(fid);
end
fclose(fid);
return;
end


function keepData = ReadFile_5(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    newStr = split(tline,' ');
    subSeqlen = str2double(newStr{5});
    timeNeed = str2double(newStr{11}); 
    
    keepData(lnCnt, 1) = subSeqlen;
    keepData(lnCnt, 2) = timeNeed;
    
    lnCnt = lnCnt +1;
    tline = fgetl(fid);
end
fclose(fid);
return;
end


function keepData = ReadFile_6(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    newStr = split(tline,' ');
    subSeqlen = str2double(newStr{6});
    timeNeed = str2double(newStr{12}); 
    
    keepData(lnCnt, 1) = subSeqlen;
    keepData(lnCnt, 2) = timeNeed;
    
    lnCnt = lnCnt +1;
    tline = fgetl(fid);
end
fclose(fid);
return;
end

function keepData = ReadFile_7(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    newStr = split(tline,' ');
    subSeqlen = str2double(newStr{6});
    timeNeed = str2double(newStr{12}); 
    
    keepData(lnCnt, 1) = subSeqlen;
    keepData(lnCnt, 2) = timeNeed;
    
    lnCnt = lnCnt +1;
    tline = fgetl(fid);
end
fclose(fid);
return;
end