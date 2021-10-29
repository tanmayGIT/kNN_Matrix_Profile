function [sorted_CrossingDisval] = Match_Subseq_N_Query_Par_SelfJoin(subFoldersTarget, subSeqLen, kNN_Uwant)

tic
 [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = ...
        ProcessTheFilesInParallel(subFoldersTarget,subSeqLen, kNN_Uwant);
    
    
% Get the top KNN from this sub-set of distances after sorting
[sorted_CrossingDisval, sorted_Crossingidx] = sort(keep_KNN_Dist_ToQuery(:, :), 1); % columnwise sorting
sorted_CrossingDisval = sorted_CrossingDisval(1:kNN_Uwant, :);
sorted_CrossingidxActual = zeros(size(sorted_Crossingidx));

for jj = 1:1:size(sorted_Crossingidx,2) % for each column
    sorted_CrossingidxActual(:,jj) = keep_KNN_Index_ToQuery(sorted_Crossingidx(1:kNN_Uwant, jj), jj);
end
wholeTime = toc;

return;
end


function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = ProcessTheFilesInParallel(keepAllTargetTogether, subSeqLen, kNN_Uwant)

[keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = CalculateMatrixProfile_ParMaxHeapBased(keepAllTargetTogether, subSeqLen, kNN_Uwant);

end


function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = CalculateMatrixProfile_ParMaxHeapBased(targetSeqComplete, subSeqLen, kNN_Uwant)

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


[target_data_freq(:, 1), target_data_mu(:, 1), target_data_sig(:, 1), targetDataCleaned] = mass_pre(targetSeqComplete(:, 1), dataLenTarget, subSeqLen);


cutTarget = targetDataCleaned(1:(1+subSeqLen-1), 1);
[~, first_prod(:, 1)] = mass(target_data_freq(:, 1), cutTarget, dataLenTarget, ...
    subSeqLen, target_data_mu(:, 1), target_data_sig(:, 1), target_data_mu(1, 1), target_data_sig(1, 1)); 

for ii = 1:1:proLen 
    cutQuery = targetDataCleaned(ii:ii+subSeqLen-1, 1);
    
    if( (target_data_sig(ii, :) ~= 0) )
        if(ii == 1)
            [dist_pro(:, 1), last_prod(:, 1)] = mass(target_data_freq(:, 1), cutQuery, dataLenTarget, ...
                subSeqLen, target_data_mu(:, 1), target_data_sig(:, 1), target_data_mu(ii, 1), target_data_sig(ii, 1)); % calculate distance profile with the very first query sequence
        else
            last_prod(2:dataLenTarget - subSeqLen + 1, :) = last_prod(1:dataLenTarget - subSeqLen, :) - targetDataCleaned(1:dataLenTarget - subSeqLen, :) ...
                .* repmat(drop_val, proLen - 1, 1) + targetDataCleaned(subSeqLen + 1:dataLenTarget, :) .* repmat(cutQuery(subSeqLen, :), proLen - 1, 1);
            
            last_prod(1, :) = first_prod(ii, :);
            dist_pro = 2 * (subSeqLen - (last_prod ...
                - subSeqLen * target_data_mu .* repmat(target_data_mu(ii, :), proLen, 1)) ...
                ./ (target_data_sig .* repmat(target_data_sig(ii, :), proLen, 1)));
        end
        dist_pro = real(dist_pro);
        
         % apply exclusion zone
        exc_st = max(1, ii - exc_zone);
        exc_ed = min(proLen, ii+exc_zone);
        dist_pro(exc_st:exc_ed, :) = inf;
        dist_pro(target_data_sig < eps) = inf;
        
        dist_pro(isnan(dist_pro(:,1)),1) = intmax('int32');  
        dist_pro(isinf(dist_pro(:,1)), 1) = intmax('int32'); 
        dist_pro(~isreal(dist_pro(:,1)), 1) = intmax('int32'); 
        dist_pro(~isfinite(dist_pro(:,1)), 1) = intmax('int32');
        dist_pro((dist_pro(:,1)< 0), 1) = intmax('int32');
        
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

