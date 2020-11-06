function [wholetime] = Match_Big_Data(subFoldersTarget, querySeqComplete, subSeqLen,...
    kNN_Uwant, whichDimToConsider)


myCluster = parcluster('local');
n_work = myCluster.NumWorkers;
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



pool = gcp('nocreate'); % Get current parallel pool
if isempty(gcp('nocreate')) % if the current parallel pool is empty
    parpool(n_work);
    fprintf('The parallel pool exist from before and the opened pool are : %d \n', n_work);
elseif (pool.NumWorkers ~= n_work)  % if the current parallel pool is not equal to available number of workers
    delete(gcp('nocreate'));
    parpool(n_work);
    fprintf('The parallel pool exist from before but it is less than available. So, the no. of newly opened pool are : %d \n', n_work);
end



keep_KNN_Cell = cell(n_work,1);
tic

parfor divFilesInCluster = 1:1:n_work
    startFiles = allStartFiles(divFilesInCluster);
    endFiles = allEndFiles(divFilesInCluster);
    
    [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery, keepDataFileInfo] = ...
        ProcessTheFilesInParallel(subFoldersTarget,querySeqComplete, ...
        subSeqLen, kNN_Uwant, whichDimToConsider, startFiles, endFiles);
    
    keep_KNN_Cell{divFilesInCluster}.DistArray =  keep_KNN_Dist_ToQuery;
    keep_KNN_Cell{divFilesInCluster}.IndexArray = keep_KNN_Index_ToQuery;
    keep_KNN_Cell{divFilesInCluster}.FileInfo = keepDataFileInfo;
end

numOfpossibleQuery = (size(querySeqComplete',1) - subSeqLen +1);

% Merge all the distances and their corresponding indexes
keep_KNN_Dist_Merge = zeros((kNN_Uwant * n_work),numOfpossibleQuery);
keep_KNN_Index_Merge = zeros((kNN_Uwant * n_work),numOfpossibleQuery);
keep_Index_LoopUp = zeros((kNN_Uwant * n_work),numOfpossibleQuery); % to keep the information about which core has been used for the operation
keep_KNN_DataFileInfo_Merge = cell((kNN_Uwant * n_work),numOfpossibleQuery);

indexMe = 1;
% cntCells = 1;
for mergeMe = 1:1:n_work
    keep_KNN_Dist_Merge(indexMe:((indexMe+kNN_Uwant)-1),:) = keep_KNN_Cell{mergeMe}.DistArray(:,:);
    keep_KNN_Index_Merge(indexMe:((indexMe+kNN_Uwant)-1),:) = keep_KNN_Cell{mergeMe}.IndexArray(:,:);
    keep_Index_LoopUp(indexMe:((indexMe+kNN_Uwant)-1),:) = mergeMe; % which core has given this output
    
    getGroupFileInfo = keep_KNN_Cell{mergeMe}.FileInfo;
    cntMe = 1;
    for io = indexMe:1:(indexMe+kNN_Uwant-1)
        for jo = 1:1:numOfpossibleQuery
            keep_KNN_DataFileInfo_Merge{io, jo} = getGroupFileInfo{cntMe,jo};
        end
        cntMe = cntMe +1;
    end
    indexMe = indexMe + kNN_Uwant;
end

[sorted_Disval, sorted_idx] = sort(keep_KNN_Dist_Merge(:, :), 1); % columnwise sorting
wholetime = toc;
% fprintf('The total time spend is : %f \n',wholetime);

keepAllNearestNNInfo = cell(kNN_Uwant,numOfpossibleQuery);
keepFileDistEntryOnly = zeros(kNN_Uwant,numOfpossibleQuery);


for goEachQuery = 1:1:numOfpossibleQuery
    for iNN = 1:1:kNN_Uwant
        getSortedIdx = sorted_idx(iNN,goEachQuery);
        if (keep_Index_LoopUp(getSortedIdx,goEachQuery) > 0) % that means it is a normal distance obtained from normal grouping operation
            pickIndex = keep_KNN_Index_Merge(getSortedIdx,goEachQuery);
            
            if(isfinite(sorted_Disval(iNN,goEachQuery))) % if the distance at that index is finite
                
                keepFileDistEntryOnly(iNN,goEachQuery) = sorted_Disval(iNN,goEachQuery);
                getAllDataFileInfo = keep_KNN_DataFileInfo_Merge{getSortedIdx,goEachQuery} ;
                
                fileNam = getAllDataFileInfo.FileNum;
                seriesStart = getAllDataFileInfo.DataStart;
                seriesEnd = getAllDataFileInfo.DataEnd;
                findFlag = false;
                if((seriesStart <= pickIndex) && (seriesEnd >= pickIndex) )
                    fileFullName = fileNam;
                    if(seriesStart == 1) % then it is the very first file treated in that specific group
                        forward = pickIndex + 1; % then it is the very first file
                    else
                        pickPrevEnd = seriesStart-1; % get the previous end
                        forward = pickIndex - pickPrevEnd;
                    end
                    keepAllNearestNNInfo{iNN, goEachQuery}.Indicator = 0; % means there is no existence of the second file
                    keepAllNearestNNInfo{iNN, goEachQuery}.NNIndex = iNN;
                    
                    keepAllNearestNNInfo{iNN, goEachQuery}.FirstFileName = fileFullName;
                    keepAllNearestNNInfo{iNN, goEachQuery}.IntialFileStart = forward;
                    keepAllNearestNNInfo{iNN, goEachQuery}.InitialFileEnd = ((forward)+subSeqLen)-1;
                    
                    findFlag = true;
                end
                if(~findFlag) % if find flag remains false then there is some problem
                    error('This flag should not remain false that means the file is not found')
                end
            end
        end
    end
end


end



function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery, keepAllFileInfo] = ...
    ProcessTheFilesInParallel(subFoldersTarget,querySeqComplete,...
    subSeqLen, kNN_Uwant, whichDimToConsider, startFiles, endFiles)


% Get all the target sequences together and merged
keepAllTargetTogether = zeros(1,1);
fullPtCnt = 1;
keepDataFileInfo = cell(1,1);

getGoodFileCnt = 1;
for lTarget = startFiles:1:endFiles % get the target files
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
        fullPtCnt = fullPtCnt + (getLengthTarget-1);
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
        fullPtCnt = fullPtCnt + (getLengthTarget-1);
    end
    getGoodFileCnt = getGoodFileCnt + 1;
end

% no sub-sequecne should be used
validSubSeq = ones(length(keepAllTargetTogether),1);
for it = 1:1:size(keepDataFileInfo,1)
    getDataEnd =  keepDataFileInfo{it,1}.DataEnd;
    lastSubSeqStart = getDataEnd -subSeqLen +1;
    blankStart = lastSubSeqStart+1;
    blankEnd = getDataEnd;
    
    validSubSeq(blankStart:blankEnd,1) = 0;
end


[keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = CalculateMatrixProfile_ParMaxHeapBased(keepAllTargetTogether, querySeqComplete', subSeqLen, ...
    kNN_Uwant,whichDimToConsider, validSubSeq);

% Now verify each index and put in the right cell order
keepAllFileInfo = cell(size(keep_KNN_Index_ToQuery,1), size(keep_KNN_Index_ToQuery,2));
for io = 1:1:size(keep_KNN_Index_ToQuery,1)
    for jo = 1:1:size(keep_KNN_Index_ToQuery,2)
        getIndex = keep_KNN_Index_ToQuery(io, jo);
        enterFlag = false;
        for pickFiles = 1:1:size(keepDataFileInfo,1)
            getFullFileInfo = keepDataFileInfo{pickFiles,1};
            
            if( (getFullFileInfo.DataStart <= getIndex) && ...
                    (getIndex <= getFullFileInfo.DataEnd) )
                keepAllFileInfo{io, jo} = getFullFileInfo;
                enterFlag = true;
                break;
            end
        end
        if(~enterFlag)
            error('I am trying to index inside the calling function and if the enter flag remains false then there is some problem');
        end
    end
end


end



% change this actually
function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = CalculateMatrixProfile_ParMaxHeapBased(targetSeqComplete, querySeqComplete, ...
    subSeqLen, kNN_Uwant, whichDimToConsider, validSubSeq)

querySeqExchanged = targetSeqComplete;
targetSeqExchanged = querySeqComplete;
% get various length
dataLenQuery = size(querySeqExchanged, 1); % it is taken as query sequence
dataLenTarget = size(targetSeqExchanged, 1); % it is taken as target sequence
proLen = dataLenTarget - subSeqLen + 1; % the distance profile will be computed using this length
nPossibleQuery = dataLenQuery - subSeqLen + 1;

keep_KNN_Dist_ToQuery = Inf(kNN_Uwant, proLen); % the row-wise I keep the nearest neighbors and columnwise for all the queries
keep_KNN_Index_ToQuery = Inf(kNN_Uwant, proLen);

sorted_Disval = Inf(kNN_Uwant+1, proLen);
sorted_idx = Inf(kNN_Uwant+1, proLen);
refinedSortIndex = Inf(kNN_Uwant+1, proLen);

% check input
%if subSeqLen > dataLenTarget / 2
 %   error(['Error: Time series is too short relative ', ...
  %      'to desired subsequence length']);
%end
if subSeqLen < 4
    error('Error: Subsequence length must be at least 4');
end

% initialization
query_data_freq = zeros((subSeqLen + dataLenQuery), 1);
query_data_mu = zeros(nPossibleQuery, 1);
query_data_sig = zeros(nPossibleQuery, 1);

target_data_freq = zeros((subSeqLen + dataLenTarget), 1);
target_data_mu = zeros(proLen, 1);
target_data_sig = zeros(proLen, 1);
drop_val = zeros(1, 1);

dist_pro = zeros(proLen, 1);
last_prod = zeros(proLen, 1);
first_prod = zeros(nPossibleQuery, 1);

% the reason of getting back the sequence vectors from the below functions i.e. "querySeqExchanged" and "targetSeqExchanged" is that after passing
% from this function, the cells having NAN or Inf values are replaced by 0.
[query_data_freq(:, 1), query_data_mu(:, 1), query_data_sig(:, 1),querySeqExchanged ] = mass_pre(querySeqExchanged(:, whichDimToConsider), dataLenQuery, subSeqLen);
[target_data_freq(:, 1), target_data_mu(:, 1), target_data_sig(:, 1), targetSeqExchanged] = mass_pre(targetSeqExchanged(:, 1), dataLenTarget, subSeqLen);

% here you need to calculate the distance between first sub-sequence of target with all the possible sub-sequence of query
cutTarget = targetSeqExchanged(1:1+subSeqLen-1, 1);
[~, first_prod(:, 1)] = mass(query_data_freq(:, 1), cutTarget, dataLenQuery, ...
    subSeqLen, query_data_mu(:, 1), query_data_sig(:, 1), target_data_mu(1, 1), target_data_sig(1, 1)); % calculate distance profile with the very first query sequence

for ii = 1:1:nPossibleQuery % total number of sub-sequnce possible in the query, remember that query is bigger in length
    
    % fprintf('I am good query subsequence at index : %d \n',ii);
    cutQuery = querySeqExchanged(ii:ii+subSeqLen-1, 1);
    
    if( (query_data_sig(ii, :) ~= 0) && (validSubSeq(ii,1) == 1) ) % otherwise all the denominator will become 0, hence all last_prod elements become -Inf
        if(ii == 1)
            [dist_pro(:, 1), last_prod(:, 1)] = mass(target_data_freq(:, 1), cutQuery, dataLenTarget, ...
                subSeqLen, target_data_mu(:, 1), target_data_sig(:, 1), query_data_mu(ii, 1), query_data_sig(ii, 1)); % calculate distance profile with the very first query sequence
        else
            last_prod(2:dataLenTarget - subSeqLen + 1, :) = last_prod(1:dataLenTarget - subSeqLen, :) - targetSeqExchanged(1:dataLenTarget - subSeqLen, :) ...
                .* repmat(drop_val, proLen - 1, 1) + targetSeqExchanged(subSeqLen + 1:dataLenTarget, :) .* repmat(cutQuery(subSeqLen, :), proLen - 1, 1);
            
            last_prod(1, :) = first_prod(ii, :);
            dist_pro = 2 * (subSeqLen - (last_prod ...
                - subSeqLen * target_data_mu .* repmat(query_data_mu(ii, :), proLen, 1)) ...
                ./ (target_data_sig .* repmat(query_data_sig(ii, :), proLen, 1)));
        end
        dist_pro = real(dist_pro);
        
        % after removing the non finite values (NAN, Inf) in the function
        % "mass_pre", if there are still some non fininte values present in the calcualted distance then
        % I replace them by very high value 100000. I do the same for any
        % negative distance value
        
        dist_pro(isnan(dist_pro(:,1)),1) = intmax('int32');  % if find nan value then replace by high value
        dist_pro(isinf(dist_pro(:,1)), 1) = intmax('int32'); % if find inf value then replace by high value
        dist_pro(~isreal(dist_pro(:,1)), 1) = intmax('int32'); % if find complex vallue then replace by high value
        dist_pro(~isfinite(dist_pro(:,1)), 1) = intmax('int32');
        dist_pro((dist_pro(:,1)< 0), 1) = intmax('int32');   % if the value if negative then replace it by a bigh value
        
        drop_val(:) = cutQuery(1, :);
        if (ii <= kNN_Uwant)
            sorted_Disval(ii,:) = sqrt(dist_pro(:));
            sorted_idx(ii,:) = ii;
        else
            
            % Get the top KNN from this sub-set of distances after sorting.
            % Now when you have filled-up the kNN_Uwant elements then you do columnwise sorting
            sorted_Disval((kNN_Uwant+1),:) = sqrt(dist_pro(:));
            sorted_idx((kNN_Uwant+1),:) = ii;
            
            if(ii == (kNN_Uwant+1)) % do it for once only
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
    end
end

keep_KNN_Dist_ToQuery(1:kNN_Uwant, :) = sorted_Disval(1:kNN_Uwant,:); % putting the final sorted distance
keep_KNN_Index_ToQuery(1:kNN_Uwant, :) = sorted_idx(1:kNN_Uwant,:); % putting the final sorted index
return ;
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
    data_len, sub_len, data_mu, data_sig, query_mu, query_sig)
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
    disp('need checking');
end
last_prod = real(product(sub_len:data_len));
end
