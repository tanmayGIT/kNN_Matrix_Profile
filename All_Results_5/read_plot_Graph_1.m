
% The following code is to plot the graph of time v/s kNN 
% or time v/s subSequence Length or time v/s time series length or 
% time v/s no. of cores. To plot each graph, you need to commnet and uncommnet 
% following block of code (seperated by several line space)

% For example, the following block of code 183-202 are uncommented. Next, you
% need to comment the lines 183-202 and uncomment the lines 77-90 to see the
% other plots.

% do the same for the other block of codes i.e. 102-127 or 137-150 etc.

% Tanmoy Mondal, Reza Akbarinia,, and Florent Masseglia, "Matrix Profile 
% Based kNN Search over Large Time Series," submitted to: 
% "Elsevier Pattern Recognition Journal", 2020.
% https://sites.google.com/view/knnmatrixprofile/home

clc;
%close all;
clear 


% Self Join Only

textFilePathSubSeq9_1 = '212_protien_subSeq_SelfJoin.txt';
textFilePathSubSeq9_2 = '2452protien_subSeq_SelfJoin_keogh.txt';

textFilePathSubSeq10_1 = '212_sheep_subSeq_SelfJoin.txt';
textFilePathSubSeq10_2 = '2452sheep_subSeq_SelfJoin_keogh.txt';

keepDataProtienSubSeq9_1 = ReadFile_5(textFilePathSubSeq9_1);
keepDataProtienSubSeq9_2 = ReadFile_5(textFilePathSubSeq9_2);

keepDataSheepSubSeq10_1 = ReadFile_5(textFilePathSubSeq10_1);
keepDataSheepSubSeq10_2 = ReadFile_5(textFilePathSubSeq10_2);


textFilePathSubSeq11_1 = '212randomWalk_kNN_SelfJoin.txt';
textFilePathSubSeq11_2 = '2452randomWalk_kNN_SelfJoin_keogh.txt';

textFilePathSubSeq12_1 = '212seismic_kNN_SelfJoin.txt';
textFilePathSubSeq12_2 = '2452seismic_kNN_SelfJoin_keogh.txt';

keepDataRandomWalk11_1 = ReadFile_6(textFilePathSubSeq11_1);
keepDataRandomWalk11_2 = ReadFile_6(textFilePathSubSeq11_2);

keepDataSeismic12_1 = ReadFile_6(textFilePathSubSeq12_1);
keepDataSeismic12_2 = ReadFile_6(textFilePathSubSeq12_2);



textFilePathSubSeq13_1 = '212randomWalk_nFiles_SelfJoin.txt';
textFilePathSubSeq13_2 = '2452randomWalk_nFiles_SelfJoin_keogh.txt';

textFilePathSubSeq14_1 = '212seismic_nFiles_SelfJoin.txt';
textFilePathSubSeq14_2 = '2452seismic_nFiles_SelfJoin_keogh.txt';

keepDataRandomWalk13_1 = ReadFile_7(textFilePathSubSeq13_1);
keepDataRandomWalk13_2 = ReadFile_7(textFilePathSubSeq13_2);

keepDataSeismic14_1 = ReadFile_7(textFilePathSubSeq14_1);
keepDataSeismic14_2 = ReadFile_7(textFilePathSubSeq14_2);


figure();
h(1) = subplot(1,1,1);

colorspec = {[0.9 0.9 0.9]; [0.8 0.8 0.8]; [0.6 0.6 0.6]; ...
  [0.4 0.4 0.4]; [0.2 0.2 0.2];[0 0.75 0.75];[0 0.5 0];[0.75 0.75 0];...
  [1 0.50 0.25];[0.75 0 0.75];[0.7 0.7 0.7];[0.8 0.7 0.6];[0.6 0.5 0.4 ]};








% q = plot(keepDataProtienSubSeq(:,1),keepDataProtienSubSeq(:,2),'r-*',keepDataSheepSubSeq(:,1),keepDataSheepSubSeq(:,2),'b-d');
% hold on;
% hleg1 = legend('Protien Data(4074 * 680)','Sheep Data(16880 * 500)');
% set(hleg1,'Location','NorthEast')
% set(hleg1,'FontSize',10)
% set(gca,'XTick', keepDataProtienSubSeq(:,1));
% grid on;
% set(gca,'FontSize',10);
% xl = xlabel('Sub-sequence Length');
% yl = ylabel('Time Needed (sec.)');
% set(xl,'FontSize',12,'FontWeight','bold','FontName','Courier');
% set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% hold off;
% disp('see me');









% q = plot(keepDataRandomWalk_kNN(:,1),keepDataRandomWalk_kNN(:,2),'r-*', ...
%         keepDataRandomWalk_kNNSortBased(:,1),keepDataRandomWalk_kNNSortBased(:,2),'r--d', ...
%         keepDataRandomWalk_kNNMaxBased(:,1),keepDataRandomWalk_kNNMaxBased(:,2),'r:s',...
%         keepDataSiesmic_kNN(:,1),keepDataSiesmic_kNN(:,2),'b-*',...
%         keepDataSiesmic_kNNSortBased(:,1),keepDataSiesmic_kNNSortBased(:,2),'b--d',...
%         keepDataSiesmic_kNNMaxBased(:,1),keepDataSiesmic_kNNMaxBased(:,2),'b:s');
%     
% hold on;
% hleg1 = legend('kNN Heap-Max-Based on RandomWalk Data(50000 * 256)', ...
%                 'kNN Sorting-Based on RandomWalk Data(50000 * 256)',... 
%                 'kNN Max-Based on RandomWalk Data(50000 * 256)',...
%                 'kNN Heap-Max-Based on Seismic Data(50000 * 200)',...
%                 'kNN Sorting-Based on Seismic Data(50000 * 200)',...
%                 'kNN Max-Based on Seismic Data(50000 * 200)');
%             
% set(hleg1,'Location','NorthWest')
% set(hleg1,'FontSize',16)
% set(gca,'XTick', keepDataRandomWalk_kNN(:,1));
% grid on;
% set(gca,'FontSize',10);
% xl = xlabel('k Nearest Neighbors');
% yl = ylabel('Time Needed (sec.)');
% set(xl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% hold off;
% disp('see me');









% q = plot(keepDataRandomWalk_NumRws(:,1),keepDataRandomWalk_NumRws(:,2),'r-*',keepDataSiesmic_NumRws(:,1),keepDataSiesmic_NumRws(:,2),'b-d');
% hold on;
% hleg1 = legend('RandomWalk Data (50000 * 256)','Seismic Data (50000 * 200)');
% set(hleg1,'Location','NorthEast')
% set(hleg1,'FontSize',14)
% set(gca,'XTick', keepDataRandomWalk_NumRws(:,1));
% grid on;
% set(gca,'FontSize',10);
% xl = xlabel('Number of Time Series');
% yl = ylabel('Time Needed (sec.)');
% set(xl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% hold off;
% disp('see me');









% q = plot(keepDataRandomWalk_Cores(2:end,1),keepDataRandomWalk_Cores(2:end,2),...
%             'r-*',keepDataSiesmic_Cores(2:end,1),keepDataSiesmic_Cores(2:end,2),'b-d');
% hold on;
% hleg1 = legend('RandomWalk Data (50000 * 256)','Seismic Data (50000 * 200)');
% set(hleg1,'Location','NorthEast')
% set(hleg1,'FontSize',16)
% set(gca,'XTick', keepDataRandomWalk_Cores(2:end,1));
% grid on;
% set(gca,'FontSize',10);
% xl = xlabel('Number of Used Cores');
% yl = ylabel('Time Needed (sec.)');
% set(xl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% hold off;
% disp('see me');








% Plot for the protien and sheep data

q = plot(keepDataProtienSubSeq9_1(:,1),keepDataProtienSubSeq9_1(:,2),'r-*', ...
        keepDataProtienSubSeq9_2(:,1),keepDataProtienSubSeq9_2(:,2),'r--d', ...
        keepDataSheepSubSeq10_1(:,1),keepDataSheepSubSeq10_1(:,2),'b-*',...
        keepDataSheepSubSeq10_2(:,1),keepDataSheepSubSeq10_2(:,2),'b--d');
    
hold on;
hleg1 = legend('Proposed Self Join on Protien Data (200 * 680 = 2.77 M)', 'Yeh et al. Self Join on  Protien Data (200 * 680 = 2.77 M)', ...
                ' Proposed Self Join on Sheep Data (200 * 500 = 8.44 M)','Yeh et al. Self Join on Sheep Data (200 * 500 = 8.44 M)');
            
set(hleg1,'Location','NorthWest')
set(hleg1,'FontSize',14)
set(gca,'XTick', keepDataProtienSubSeq9_1(:,1));
grid on;
set(gca,'FontSize',10);
xl = xlabel('Sub-sequence length');
yl = ylabel('Time Needed (sec.)');
set(xl,'FontSize',14,'FontWeight','bold','FontName','Courier');
set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
hold off;
disp('see me');







% q = plot(keepDataRandomWalk11_1(:,1),keepDataRandomWalk11_1(:,2),'r-*', ...
%         keepDataRandomWalk11_2(:,1),keepDataRandomWalk11_2(:,2),'r--d', ...
%         keepDataSeismic12_1(:,1),keepDataSeismic12_1(:,2),'b-*',...
%         keepDataSeismic12_2(:,1),keepDataSeismic12_2(:,2),'b--d');
%     
% hold on;
% hleg1 = legend('Proposed Self Join on RandomWalk Data (1000 * 256 = 256,000)', 'Yeh et al. Self Join on  RandomWalk Data (1000 * 256 = 256,000)', ...
%                 ' Proposed Self Join on Seismic Data (1000 * 200 = 200,000)','Yeh et al. Self Join on Seismic Data (1000 * 200 = 200,000)');
%             
% set(hleg1,'Location','NorthWest')
% set(hleg1,'FontSize',14)
% set(gca,'XTick', keepDataRandomWalk11_1(:,1));
% grid on;
% set(gca,'FontSize',10);
% xl = xlabel('k Nearest Neighbors');
% yl = ylabel('Time Needed (sec.)');
% set(xl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
% hold off;
% disp('see me');








% q = plot(keepDataRandomWalk13_1(:,1),keepDataRandomWalk13_1(:,2),'r-*', ...
%         keepDataRandomWalk13_2(:,1),keepDataRandomWalk13_2(:,2),'r--d', ...
%         keepDataSeismic14_1(:,1),keepDataSeismic14_1(:,2),'b-*',...
%         keepDataSeismic14_2(:,1),keepDataSeismic14_2(:,2),'b--d');
%     
% hold on;
% hleg1 = legend('Proposed Self Join on RandomWalk Data (no. of time series * 256)', 'Yeh et al. Self Join on  RandomWalk Data (no. of time series * 256)', ...
%                 ' Proposed Self Join on Seismic Data (no. of time series * 200)','Yeh et al. Self Join on Seismic Data (no. of time series * 200)');
%             
% set(hleg1,'Location','NorthWest')
% set(hleg1,'FontSize',14)
% set(gca,'XTick', keepDataRandomWalk13_1(:,1));
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