function wholetime = Match_SelfJoin_VariableCore(subFoldersTarget, subSeqLen, kNN_Uwant, nCores)

% myCluster = parcluster('local');
n_work = nCores; % myCluster.NumWorkers;
fprintf('The Number of Worker present is : %d \n',n_work);

% creating an array to keep the block of files together
totNumFiles = size(subFoldersTarget,1);
remainderFiles = rem(totNumFiles,n_work);
keepDivideClusInfo = zeros(n_work,2);

if(remainderFiles == 0) % the perfect condition
    startFiles = 1;
    %   endFiles = size(subFoldersTarget,1);
    gapFiles = totNumFiles / n_work;
    for iClus = 1:1:n_work
        keepDivideClusInfo(iClus,1) = startFiles;
        keepDivideClusInfo(iClus,2) = (startFiles + gapFiles)-1;
        startFiles = startFiles + gapFiles;
    end
elseif (remainderFiles >0)
    startFiles = 1;
    gapFiles = (totNumFiles - remainderFiles) / n_work;
    %     endFiles = size(subFoldersTarget,1) -(gapFiles+remainderFiles);

    for iClus = 1:1:(n_work-1) % last cell is kept for remainder filese
        keepDivideClusInfo(iClus,1) = startFiles;
        keepDivideClusInfo(iClus,2) = (startFiles + gapFiles)-1;
        startFiles = startFiles + gapFiles;
    end
    % putting all the files in the last cell
    keepDivideClusInfo(end,1) = startFiles;
    keepDivideClusInfo(end,2) = size(subFoldersTarget,1);
end

allStartFiles = keepDivideClusInfo(:,1);
allEndFiles = keepDivideClusInfo(:,2);

% myCluster = parcluster('local');
% getAllCrores = myCluster.NumWorkers;

pool = gcp('nocreate'); % Get current parallel pool
if isempty(gcp('nocreate')) % if the current parallel pool is empty
    parpool(n_work);
    fprintf('The parallel pool exist from before and the opened pool are : %d \n', n_work);
elseif (pool.NumWorkers ~= n_work)  % if the current parallel pool is not equal to available number of workers I wanted
    delete(gcp('nocreate'));
    parpool(n_work);
    fprintf('The parallel pool exist from before but it is less than available. So, the no. of newly opened pool are : %d \n', n_work);
end

tic
[keepAllTargetTogether, keepDataFileInfo] = ConcatenateAllSeries(subFoldersTarget, subSeqLen); % concatenate all the time series
% fprintf('The data is concatenated, the size of the data is : %d %d', size(keepAllTargetTogether, 1), size(keepAllTargetTogether, 2));
keep_KNN_Cell = cell(n_work,1);

parfor divFilesInCluster = 1:1:n_work
    startFiles = allStartFiles(divFilesInCluster);
    endFiles = allEndFiles(divFilesInCluster);

    [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = ...
        ProcessTheFilesInParallel(keepAllTargetTogether, keepDataFileInfo, subSeqLen, kNN_Uwant, startFiles, endFiles);

    keep_KNN_Cell{divFilesInCluster}.DistArray =  keep_KNN_Dist_ToQuery;
    keep_KNN_Cell{divFilesInCluster}.IndexArray = keep_KNN_Index_ToQuery;
end

numOfpossibleQuery = (size(keepAllTargetTogether,1) - subSeqLen +1);

% Merge all the distances and their corresponding indexes
keep_KNN_Dist_Merge = zeros((kNN_Uwant * n_work),numOfpossibleQuery);
keep_KNN_Index_Merge = zeros((kNN_Uwant * n_work),numOfpossibleQuery);
keep_Index_LoopUp = zeros((kNN_Uwant * n_work),numOfpossibleQuery); % to keep the information about which core has been used for the operation

indexMe = 1;
% cntCells = 1;
for mergeMe = 1:1:n_work
    keep_KNN_Dist_Merge(indexMe:((indexMe+kNN_Uwant)-1),:) = keep_KNN_Cell{mergeMe}.DistArray(:,:);
    keep_KNN_Index_Merge(indexMe:((indexMe+kNN_Uwant)-1),:) = keep_KNN_Cell{mergeMe}.IndexArray(:,:);
    keep_Index_LoopUp(indexMe:((indexMe+kNN_Uwant)-1),:) = mergeMe; % which core has given this output

    indexMe = indexMe + kNN_Uwant;
end

[sorted_Disval, sorted_idx] = sort(keep_KNN_Dist_Merge(:, :), 1); % columnwise sorting

wholetime = toc;

return;
end






function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = ...
                ProcessTheFilesInParallel(keepAllTargetTogether, keepDataFileInfoConcat, subSeqLen, kNN_Uwant, startMe, endMe)

 startFiles = keepDataFileInfoConcat{startMe}.DataStart;
 endFiles = keepDataFileInfoConcat{endMe}.DataEnd;

[keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = SingleSeries_ParMaxHeapBased(keepAllTargetTogether, subSeqLen, kNN_Uwant, startFiles, endFiles);

keep_KNN_Dist_ToQuery = keep_KNN_Dist_ToQuery(:, :);
keep_KNN_Index_ToQuery = keep_KNN_Index_ToQuery(:, :);
end


function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = SingleSeries_ParMaxHeapBased(targetSeqComplete, subSeqLen, kNN_Uwant, startFiles, endFiles)

exc_zone = round(subSeqLen / 2);
dataLenTarget = size(targetSeqComplete, 1);
proLen = dataLenTarget - subSeqLen + 1;


keep_KNN_Dist_ToQuery = Inf(kNN_Uwant, proLen);
keep_KNN_Index_ToQuery = Inf(kNN_Uwant, proLen);

sorted_Disval = Inf(kNN_Uwant+1, proLen);
sorted_idx = Inf(kNN_Uwant+1, proLen);

if subSeqLen > dataLenTarget / 2
    error(['Error: Time series is too short relative ', ...
        'to desired subsequence length']);
end
if subSeqLen < 4
    error('Error: Subsequence length must be at least 4');
end

target_data_freq = zeros((subSeqLen + dataLenTarget), 1);
target_data_mu = zeros(proLen, 1);
target_data_sig = zeros(proLen, 1);
drop_val = zeros(1, 1);

dist_pro = zeros(proLen, 1);
last_prod = zeros(proLen, 1);
first_prod = zeros(proLen, 1);

% tic
[target_data_freq(:, 1), target_data_mu(:, 1), target_data_sig(:, 1), targetDataCleaned] = mass_pre(targetSeqComplete(:, 1), dataLenTarget, subSeqLen);


cutTarget = targetDataCleaned(startFiles:(startFiles+subSeqLen-1), 1);
[~, first_prod(:, 1)] = mass(target_data_freq(:, 1), cutTarget, dataLenTarget, ...
    subSeqLen, target_data_mu(:, 1), target_data_sig(:, 1), target_data_mu(1, 1), target_data_sig(1, 1), 226);
loopCnt = 1;

% fprintf('I am executing the group which is starting from %d and endign at %d',startFiles, (endFiles-subSeqLen+1) );


for ii = startFiles:1:(endFiles-subSeqLen+1)
    cutQuery = targetDataCleaned(ii:ii+subSeqLen-1, 1);


    if(ii == 1)
        [dist_pro(:, 1), last_prod(:, 1)] = mass(target_data_freq(:, 1), cutQuery, dataLenTarget, ...
            subSeqLen, target_data_mu(:, 1), target_data_sig(:, 1), target_data_mu(ii, 1), target_data_sig(ii, 1), 234); % calculate distance profile with the very first query sequence
    else
        last_prod(2:dataLenTarget - subSeqLen + 1, :) = last_prod(1:dataLenTarget - subSeqLen, :) - targetDataCleaned(1:dataLenTarget - subSeqLen, :) ...
            .* repmat(drop_val, proLen - 1, 1) + targetDataCleaned(subSeqLen + 1:dataLenTarget, :) .* repmat(cutQuery(subSeqLen, :), proLen - 1, 1);

        last_prod(1, :) = first_prod(ii, :);
        dist_pro = 2 * (subSeqLen - (last_prod ...
            - subSeqLen * target_data_mu .* repmat(target_data_mu(ii, :), proLen, 1)) ...
            ./ (target_data_sig .* repmat(target_data_sig(ii, :), proLen, 1)));
    end
    drop_val(:) = cutQuery(1, :);
    if( (target_data_sig(ii, :) ~= 0) )
        dist_pro = abs(dist_pro);
        dist_pro = real(dist_pro);

        % apply exclusion zone
        exc_st = max(1, ii - exc_zone);
        exc_ed = min(proLen, ii+exc_zone);
        dist_pro(exc_st:exc_ed, :) = inf;
        dist_pro(target_data_sig < eps) = inf;

%         dist_pro(isnan(dist_pro(:,1)),1) = intmax('int32');
        dist_pro(isinf(dist_pro(:,1)), 1) = intmax('int32');
%         dist_pro(~isreal(dist_pro(:,1)), 1) = intmax('int32');
%         dist_pro(~isfinite(dist_pro(:,1)), 1) = intmax('int32');
        dist_pro((dist_pro(:,1)< 0), 1) = intmax('int32');

        if (loopCnt <= kNN_Uwant)
            sorted_Disval(loopCnt,:) = sqrt(dist_pro(:));
            sorted_idx(loopCnt,:) = ii;
        else

            % Get the top KNN from this sub-set of distances after sorting.
            % Now when you have filled-up the kNN_Uwant elements then you do columnwise sorting
            sorted_Disval((kNN_Uwant+1),:) = sqrt(dist_pro(:));
            sorted_idx((kNN_Uwant+1),:) = ii;

            if(loopCnt == (kNN_Uwant+1)) % do it for once only
                for ioT = 1:1:size(sorted_Disval,2)
                    [sorted_Disval(1:kNN_Uwant, ioT), heapSortedIndexes]  = ...
                        HeapBasedPriorityQueue.buildMaxBasedHeapOnArray(sorted_Disval(1:kNN_Uwant, ioT), kNN_Uwant);
                    sorted_idx(1:kNN_Uwant, ioT) = sorted_idx (heapSortedIndexes(1:kNN_Uwant), ioT);
                end
            end
            % [maxDisvalForeachCol, maxRwidxForEachCol] = max(sorted_Disval(:, :),[], 1); % columnwise sorting
            [~, idxMax] = find(sorted_Disval(1,:) > sorted_Disval((kNN_Uwant+1),:));

            if(~isempty(idxMax))
                idxMaxActual = idxMax + (1-1);

                % replacing the maximum value
                %for iKO = 1:1:length(idxMaxActual)
                sorted_Disval(1, idxMaxActual) = sorted_Disval((kNN_Uwant+1),idxMaxActual);  % replacing the 11th row
                sorted_idx(1, idxMaxActual) = sorted_idx((kNN_Uwant+1),idxMaxActual);
                %end

                % now only run heap based maximum on those column which are modified to obtain again the maximum
                for ioT = 1:1:length(idxMaxActual)
                    [sorted_Disval(1:kNN_Uwant, idxMaxActual(ioT)), heapSortedIndexes]  = ...
                        HeapBasedPriorityQueue.buildMaxBasedHeapOnArray(sorted_Disval(1:kNN_Uwant, idxMaxActual(ioT)), kNN_Uwant);
                    sorted_idx(1:kNN_Uwant, idxMaxActual(ioT)) = sorted_idx(heapSortedIndexes(1:kNN_Uwant), idxMaxActual(ioT));
                end

            end
        end

        loopCnt = loopCnt+1;
    end
end

keep_KNN_Dist_ToQuery(1:kNN_Uwant, :) = sorted_Disval(1:kNN_Uwant,:); % putting the final sorted distance
keep_KNN_Index_ToQuery(1:kNN_Uwant, :) = sorted_idx(1:kNN_Uwant,:); % putting the final sorted index
% groupTime = toc;
% fprintf('This group is calculated and the indexes which are calculated are from %d to %d and the time needed %f \n', startFiles, (endFiles-subSeqLen+1), groupTime);
return ;
end

function [keepAllTargetTogether, keepDataFileInfo] = ConcatenateAllSeries(subFoldersTarget, subSeqLen )
% Get all the target sequences together and merged
keepAllTargetTogether = zeros(1,1);
fullPtCnt = 1;
keepDataFileInfo = cell(1,1);

getGoodFileCnt = 1;
for lTarget = 1:1:size(subFoldersTarget,1) % get the target files
    C1Target = subFoldersTarget(lTarget,:);
    getLengthTarget = length(C1Target);

    C1TargetArr = zeros(getLengthTarget,1);

    C1TargetArr(:,1) = C1Target(1, :);

    clearvars C1Target

    if (getGoodFileCnt == 1)
        keepAllTargetTogether(1:(getLengthTarget),1) = C1TargetArr(1:end);

        keepDataFileInfo{getGoodFileCnt,1}.FileNum = lTarget;
        keepDataFileInfo{getGoodFileCnt,1}.DataStart = 1;
        keepDataFileInfo{getGoodFileCnt,1}.DataEnd = (getLengthTarget);
        fullPtCnt = fullPtCnt + (getLengthTarget);
    else
        keepAllTargetBackup = keepAllTargetTogether;
        keepAllTargetTogether = zeros ((size(keepAllTargetBackup,1)+size(C1TargetArr,1)),1);
        keepAllTargetTogether(1:size(keepAllTargetBackup,1),:) = keepAllTargetBackup(:,:);

        clearvars keepAllTargetBackup
        keepAllTargetTogether(fullPtCnt:(fullPtCnt+(getLengthTarget-1)),1) = C1TargetArr(:,1);

        keepDataFileInfo{getGoodFileCnt,1}.FileNum = lTarget;
        keepDataFileInfo{getGoodFileCnt,1}.DataStart = fullPtCnt;
        keepDataFileInfo{getGoodFileCnt,1}.DataEnd = (fullPtCnt+(getLengthTarget-1));

        if(( keepDataFileInfo{getGoodFileCnt,1}.DataEnd - ...
                keepDataFileInfo{getGoodFileCnt,1}.DataStart) > getLengthTarget)
            error('need to check it, there could be some problem here');
        end
        fullPtCnt = fullPtCnt + (getLengthTarget);
    end
    getGoodFileCnt = getGoodFileCnt + 1;
end

% no sub-sequecne should be used
% validSubSeq = ones(length(keepAllTargetTogether),1);
% for it = 1:1:size(keepDataFileInfo,1)
%     getDataEnd =  keepDataFileInfo{it,1}.DataEnd;
%     lastSubSeqStart = getDataEnd -subSeqLen +1;
%     blankStart = lastSubSeqStart+1;
%     blankEnd = getDataEnd;
%
%     validSubSeq(blankStart:blankEnd,1) = 0;
% end
return
end

% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [data_freq, data_mu, data_sig,data] = mass_pre(data, data_len, sub_len)

% nanIndexesTar = find(isnan(data(:) ));
% if(~isempty(nanIndexesTar) )
%     data(nanIndexesTar,1) = 0;
% end
% infIndexesTar = find(isinf(data(:) ));
% if(~isempty(infIndexesTar) )
%     data(infIndexesTar,1) = 0;
% end
% [rowNonInf,~,~] = find(~isfinite(data));
% data(rowNonInf) = 0; % making NAN and INF to zero to avaoid the problems

nanIndexesTar = (isnan(data(:) ));
infIndexesTar = (isinf(data(:) ));
rowNonInf = (~isfinite(data));

allBadCell = bitor((bitor(nanIndexesTar,infIndexesTar)),rowNonInf);
tAllIndex = 1:numel(data(:));
data(allBadCell) = interp1(tAllIndex(~allBadCell), data(~allBadCell), tAllIndex(allBadCell));

data(data_len+1:(sub_len+data_len)) = 0;
data_freq = fft(data);
data_cum = cumsum(data);
data2_cum =  cumsum(data.^2);
data2_sum = data2_cum(sub_len:data_len) - ...
    [0; data2_cum(1:data_len-sub_len)];
data_sum = data_cum(sub_len:data_len) - ...
    [0; data_cum(1:data_len-sub_len)];
data_mu = data_sum./sub_len;
data_sig2 = (data2_sum./sub_len)-(data_mu.^2);
data_sig2 = real(data_sig2);
data_sig2 = max(data_sig2, 0);
data_sig = sqrt(data_sig2);
end


function [dist_pro, last_prod] = mass(data_freq, query, ...
    data_len, sub_len, data_mu, data_sig, query_mu, query_sig, lnNum)
% pre-process query for fft
query = query(end:-1:1);
query(sub_len+1:(sub_len+data_len)) = 0;

% compute the product
query_freq = fft(query);
product_freq = data_freq.*query_freq;
product = ifft(product_freq);

% compute the distance profile
dist_pro = 2 * (sub_len - ...
    (product(sub_len:data_len) - sub_len*data_mu*query_mu)./...
    (data_sig * query_sig));
if(~isfinite(dist_pro))
    fprintf('The bad line number %d \n', lnNum);
     %fprintf('The data_freq is \n');
     %disp(data_freq');

    % fprintf('The query is \n');
     %disp(query');

    % fprintf('The data_mu %d data_sig %d query_mu %d query_sig %d \n', data_mu, data_sig, query_mu, query_sig);
    disp('need checking');
end
last_prod = real(product(sub_len:data_len));
end
