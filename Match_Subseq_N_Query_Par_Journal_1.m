function Match_Subseq_N_Query_Par_Journal_1(subFoldersTarget, querySeqComplete, subSeqLen,...
    kNN_Uwant, whichDimToConsider)

%  THIS IS SAME VERSION AS "Match_1_Subseq_N_Query_Par_Ver1", THE ONLY
%  DIFFERENCE IS HERE i MAKE THE INDEX OF ALL COMPONENT INFO IN THE
%  FUNCTION "ProcessTheFilesInParallel" AND GET THE SOTER COMPINFO FROM ALL
%  THE GROUP OF FILES

% In this function we are considering that when we are grouping several
% files together then concatenating them together, it may happen that the
% best sub-sequence which starts from at the end of one file and ends at
% the beginning of another file. So, we are handling this kind of
% situations here

% Moreover, here we also have considered that there can be situation where
% the best match has occurred between two different groups. That means the last
% file of a group and first file of another group and the time series
% built by this can contain a best match.

myCluster = parcluster('local');
n_work = myCluster.NumWorkers;

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


% ## setup pool
if isempty(which('parpool'))
    if matlabpool('size') <= 0 %#ok<*DPOOL>
        matlabpool(n_work);
    elseif matlabpool('size')~= n_work
        matlabpool('close');
        matlabpool(n_work);
    end
else
    pool = gcp('nocreate');
    if isempty(gcp('nocreate'))
        parpool(n_work);
    elseif pool.NumWorkers ~= n_work
        delete(gcp('nocreate'));
        parpool(n_work);
    end
end

keep_KNN_Cell = cell(n_work,1);
tic
parfor divFilesInCluster = 1:1:n_work
    startFiles = allStartFiles(divFilesInCluster);
    endFiles = allEndFiles(divFilesInCluster);
    
    [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery, keepDataFileInfo] = ...
        ProcessTheFilesInParallel(subFoldersTarget,startFiles, endFiles, querySeqComplete, ...
        subSeqLen, kNN_Uwant, whichDimToConsider);
    
    keep_KNN_Cell{divFilesInCluster}.DistArray =  keep_KNN_Dist_ToQuery;
    keep_KNN_Cell{divFilesInCluster}.IndexArray = keep_KNN_Index_ToQuery;
    keep_KNN_Cell{divFilesInCluster}.FileInfo = keepDataFileInfo;

end
querySeqComplete = querySeqComplete';
numOfpossibleQuery = (size(querySeqComplete,1) - subSeqLen +1);

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

% Get the top KNN from this sub-set of distances after sorting
lenNormal = size(keep_KNN_Dist_Merge,1);
mergeCompleteAllDist = zeros(lenNormal, numOfpossibleQuery);
mergeCompleteAllIndexes = zeros(lenNormal, numOfpossibleQuery);
mergeCompleteAllInfo = cell(lenNormal, numOfpossibleQuery);

mergeCompleteAllDist(1:lenNormal,:) = keep_KNN_Dist_Merge(1:lenNormal,:);
mergeCompleteAllIndexes(1:lenNormal,:) = keep_KNN_Index_Merge(1:lenNormal,:);


for io = 1:1:lenNormal
    for getAllQuery = 1:1:numOfpossibleQuery
        mergeCompleteAllInfo{io, getAllQuery} = keep_KNN_DataFileInfo_Merge{io, getAllQuery};
    end
end



% remember this sorted index is extracted from the groups which are containing files. This is not
% for all the files which belongs to a single group
[sorted_Disval, sorted_idx] = sort(mergeCompleteAllDist(:, :), 1); % columnwise sorting
keepAllNearestNNInfo = cell(kNN_Uwant,numOfpossibleQuery);

keepFileIndxEntryOnly = zeros(kNN_Uwant,numOfpossibleQuery);
keepFileDistEntryOnly = zeros(kNN_Uwant,numOfpossibleQuery);

for goEachQuery = 1:1:numOfpossibleQuery
    %  fprintf('We are considering the query number_%d \n\n', goEachQuery);
    for iNN = 1:1:kNN_Uwant
        getSortedIdx = sorted_idx(iNN,goEachQuery);
        pickIndex = mergeCompleteAllIndexes(getSortedIdx,goEachQuery);
        
        if(isfinite(sorted_Disval(iNN,goEachQuery))) % if the distance at that index is finite
            getAllDataFileInfo = mergeCompleteAllInfo{getSortedIdx,goEachQuery} ;
            
            seriesStart = getAllDataFileInfo.DataStart;
            seriesEnd = getAllDataFileInfo.DataEnd;
            findFlag = false;
            if((seriesStart <= pickIndex) && (seriesEnd >= pickIndex) )
                if(seriesStart == 1) % then it is the very first file treated in that specific group
                    forward = pickIndex + 1; % then it is the very first file
                else
                    pickPrevEnd = seriesStart-1; % get the previous end
                    forward = pickIndex - pickPrevEnd;
                end
                
                keepAllNearestNNInfo{iNN, goEachQuery}.NNIndex = iNN;
                keepAllNearestNNInfo{iNN, goEachQuery}.FileEntryIndx = getAllDataFileInfo.FileNum;
                keepAllNearestNNInfo{iNN, goEachQuery}.IntialFileStart = forward;
                keepAllNearestNNInfo{iNN, goEachQuery}.InitialFileEnd = ((forward)+subSeqLen)-1;
                
                keepFileIndxEntryOnly(iNN, goEachQuery) = getAllDataFileInfo.FileNum;
                keepFileDistEntryOnly(iNN, goEachQuery) = sorted_Disval(iNN,goEachQuery);
                
                findFlag = true;
            end
            if(~findFlag) % if find flag remains false then there is some problem
                error('This flag should not remain false that means the file is not found')
            end
        end
    end
end
toc

 PlotTheDataJournal_2(querySeqComplete, keepAllNearestNNInfo,subFoldersTarget, subSeqLen, 324);
 % plotTheGraph(keepFileDistEntryOnly, kNN_Uwant);

end




function plotTheGraph(keepFileDistEntryOnly, kNN)

colorString = {'y','m','r','b','k','c','g',[.5 .6 .7],[.8 .2 .6], [0.72 0.52 0.04], [0.6 0.8 0.2], [0.619 1 0.4 ],...
                        [0.9 0.9 0.9], [0.4 0.4 0.4], [0.8 0.8 0.8], [0.6 0.6 0.6], [0.2 0.2 0.2], [0 0.75 0.75], [0 0.5 0],...
                        [0.4 0.58 0.9], [0.75 0 0.75], [0.8 0.7 0.6], [0.6 0.5 0.4 ], [0.8 0.6 1 ], [0 1 1], [1 0.6 0.8]};
hFig1 = figure();

hold on;
for iOut = 1:1:kNN
    plot(1:length(keepFileDistEntryOnly(iOut, :)), keepFileDistEntryOnly(iOut, :),'color', colorString{iOut}, 'LineWidth',1);  
end
hold off;



end


function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery, keepAllFileInfo, firstFileData, lastFileData] = ...
    ProcessTheFilesInParallel(subFoldersTarget,startFiles, endFiles, querySeqComplete,...
    subSeqLen, kNN_Uwant, whichDimToConsider)


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
if subSeqLen > dataLenTarget / 2
    error(['Error: Time series is too short relative ', ...
        'to desired subsequence length']);
end
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

% change this actually
function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = CalculateMatrixProfile_ParMax(targetSeqComplete, querySeqComplete, ...
    subSeqLen, kNN_Uwant, whichDimToConsider)

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
if subSeqLen > dataLenTarget / 2
    error(['Error: Time series is too short relative ', ...
        'to desired subsequence length']);
end
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
    
    if(query_data_sig(ii, :) ~= 0) % otherwise all the denominator will become 0, hence all last_prod elements become -Inf
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
            [maxDisvalForeachCol, maxRwidxForEachCol] = max(sorted_Disval(:, :),[], 1); % columnwise sorting
            [~, idxMax] = find(maxDisvalForeachCol > sorted_Disval((kNN_Uwant+1),:));
            idxMaxActual = idxMax + (1-1);
            
            % replacing the maximum value
            for iKO = 1:1:length(idxMaxActual)
                sorted_Disval(maxRwidxForEachCol(idxMax(iKO)), idxMaxActual(iKO)) = sorted_Disval((kNN_Uwant+1),idxMaxActual(iKO));  % replacing the 11th row
                sorted_idx(maxRwidxForEachCol(idxMax(iKO)),idxMaxActual(iKO)) = sorted_idx((kNN_Uwant+1),idxMaxActual(iKO));
            end
            
            % now just take top kNN_Uwant for your need
            sorted_Disval(1:kNN_Uwant,:) = sorted_Disval(1:kNN_Uwant,:);  % just keep changing the top 10 elements bcoz the last element will be modified anyhow in the next iteration
            sorted_idx(1:kNN_Uwant,:) = sorted_idx(1:kNN_Uwant,:); % just keep changing the top 10 elements bcoz the last element will be modified anyhow in the next iteration
        end
    end
end

keep_KNN_Dist_ToQuery(1:kNN_Uwant, :) = sorted_Disval(1:kNN_Uwant,:); % putting the final sorted distance
keep_KNN_Index_ToQuery(1:kNN_Uwant, :) = sorted_idx(1:kNN_Uwant,:); % putting the final sorted index
return ;
end


% change this actually
function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = CalculateMatrixProfile_Par(targetSeqComplete, querySeqComplete, ...
    subSeqLen, kNN_Uwant, whichDimToConsider)

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
if subSeqLen > dataLenTarget / 2
    error(['Error: Time series is too short relative ', ...
        'to desired subsequence length']);
end
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
    
    if(query_data_sig(ii, :) ~= 0) % otherwise all the denominator will become 0, hence all last_prod elements become -Inf
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
            [sorted_CrossingDisval, sorted_Crossingidx] = sort(sorted_Disval(:, :), 1); % columnwise sorting
            
            for io = 1:1:proLen
                refinedSortIndex(:,io) = sorted_idx(sorted_Crossingidx(:,io),io);
            end
            
            % now just take top kNN_Uwant for your need
            sorted_Disval(1:kNN_Uwant,:) = sorted_CrossingDisval(1:kNN_Uwant,:);  % just keep changing the top 10 elements bcoz the last element will be modified anyhow in the next iteration
            sorted_idx(1:kNN_Uwant,:) = refinedSortIndex(1:kNN_Uwant,:); % just keep changing the top 10 elements bcoz the last element will be modified anyhow in the next iteration
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