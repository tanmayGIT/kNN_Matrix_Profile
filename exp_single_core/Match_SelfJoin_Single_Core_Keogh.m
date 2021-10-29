function wholetime = Match_SelfJoin_Single_Core_Keogh(subFoldersTarget, subSeqLen, kNN_Uwant)

[keepAllTargetTogether, ~] = ConcatenateAllSeries(subFoldersTarget, subSeqLen); % concatenate all the time series

tic
[~, ~] = ProcessTheFilesInParallel(keepAllTargetTogether, subSeqLen, kNN_Uwant);
wholetime = toc;

return;
end


function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = ...
                ProcessTheFilesInParallel(keepAllTargetTogether, subSeqLen, kNN_Uwant)

[keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = SingleSeries_Par_STOMP_MP(keepAllTargetTogether, subSeqLen, kNN_Uwant, 1, length(keepAllTargetTogether));

end


function [keep_KNN_Dist_ToQuery, keep_KNN_Index_ToQuery] = SingleSeries_Par_STOMP_MP(targetSeqComplete, subSeqLen, kNN_Uwant, startFiles, endFiles)

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

  pro_mul = Inf(proLen, 1);
  pro_idx = zeros(proLen, 1);

  dist_pro = zeros(proLen, 1);
  last_prod = zeros(proLen, 1);
  first_prod = zeros(proLen, 1);


  [target_data_freq(:, 1), target_data_mu(:, 1), target_data_sig(:, 1), targetDataCleaned] = mass_pre(targetSeqComplete(:, 1), dataLenTarget, subSeqLen);

  validSubSeq = ones(length(targetSeqComplete),1);
  for knn = 1:1:kNN_Uwant

      cutTarget = targetDataCleaned(startFiles:(startFiles+subSeqLen-1), 1);
      [~, first_prod(:, 1)] = mass(target_data_freq(:, 1), cutTarget, dataLenTarget, ...
          subSeqLen, target_data_mu(:, 1), target_data_sig(:, 1), target_data_mu(1, 1), target_data_sig(1, 1));

      for ii = startFiles:1:(endFiles-subSeqLen+1)
          cutQuery = targetDataCleaned(ii:ii+subSeqLen-1, 1);

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
          drop_val(:) = cutQuery(1, :);

          if( (target_data_sig(ii, :) ~= 0) && (validSubSeq(ii, :) ~= 0))
              dist_pro = abs(dist_pro);
              dist_pro = real(dist_pro);

              %         dist_pro(isnan(dist_pro(:,1)),1) = intmax('int32');
              dist_pro(isinf(dist_pro(:,1)), 1) = intmax('int32');
              %         dist_pro(~isreal(dist_pro(:,1)), 1) = intmax('int32');
              %         dist_pro(~isfinite(dist_pro(:,1)), 1) = intmax('int32');
              dist_pro((dist_pro(:,1)< 0), 1) = intmax('int32');

              % apply exclusion zone
              exc_st = max(1, ii - exc_zone);
              exc_ed = min(proLen, ii+exc_zone);
              dist_pro(exc_st:exc_ed, :) = inf;
              dist_pro(target_data_sig < eps) = inf;

              dist_pro = sqrt(dist_pro);

              updatePos = dist_pro < pro_mul;
              pro_idx(updatePos) = ii;
              pro_mul(updatePos) = dist_pro(updatePos);

  %             [min_val, min_idx] = min(dist_pro);
  %             pro_mul(ii, 1) = min_val;
  %             pro_idx(ii, 1) = min_idx;
          end
      end
      sorted_Disval(knn,:) = pro_mul(:);
      sorted_idx(knn,:) = pro_idx(:);

      for iGet = 1:1:size(pro_idx,1)
          if(pro_idx(iGet) > 0 )
              validSubSeq(pro_idx(iGet),1) = 0;
          end
      end
  end


  keep_KNN_Dist_ToQuery(1:kNN_Uwant, :) = sqrt(sorted_Disval(1:kNN_Uwant,:)); % putting the final sorted distance
  keep_KNN_Index_ToQuery(1:kNN_Uwant, :) = sorted_idx(1:kNN_Uwant,:); % putting the final sorted index
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

return
end


% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [data_freq, data_mu, data_sig,data] = mass_pre(data, data_len, sub_len)

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
    disp('need checking');
end
last_prod = real(product(sub_len:data_len));
end
