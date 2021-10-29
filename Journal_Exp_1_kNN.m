function Journal_Exp_1_kNN(dataPath, textfileToSavePath, whichAlgoToConsider)
	% This experiment is for varying number of cores and it is for independent join 

	subSeqLen = 50;  % The length of the sub-sequence
	whichDimToConsider = 1;
	kNN_Uwant = 5;   % How many kNN do you want, just tell me



	% Loading the sheep data
	load(dataPath); % load('seismic_50000.mat');

	keepAllData = seeMe1;
	clear seeMe1;
	disp(size(keepAllData));

	queryIndex = 307;
	queryComplete = keepAllData(queryIndex, 1:end-1);

	noOfZerosInQuery = sum(queryComplete(:)==0);
	if(noOfZerosInQuery > (length(queryComplete) /2)) % if the number of zeros is more than half the length of query
	    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
	end

	textFileNam = textfileToSavePath; % 'seismic_kNN_SortBased.txt';

	if (exist(textFileNam, 'file') == 0)
	    disp('File does not exist, creating file.')
	    fid = fopen( textFileNam, 'w' );
	else % if the file exist
	    disp('File exists.');
	    fid = fopen(textFileNam, 'wt' ); % create a new file each time 
	end


	keepTargetData = zeros((size(keepAllData,1)-1), length(queryComplete));
	keepTargetLabels = zeros((size(keepAllData,1)-1), 1);

	if(queryIndex == 1)
	    keepTargetData(:,:) = keepAllData(2:end, 1:end-1);
	    keepTargetLabels(:,1) = keepAllData(2:end, end);
	elseif(queryIndex == size(keepAllData,1))
	    keepTargetData(:,:) = keepAllData(1:end-1, 1:end-1);
	    keepTargetLabels(:,1) = keepAllData(1:end-1, end);
	else
	    keepTargetData(1:queryIndex-1,:) = keepAllData(1:queryIndex-1, 1:end-1);
	    keepTargetData(queryIndex:end,:) = keepAllData(queryIndex+1:end, 1:end-1);
	    
	    keepTargetLabels(1:queryIndex-1,1) = keepAllData(1:queryIndex-1, end);
	    keepTargetLabels(queryIndex:end,1) = keepAllData(queryIndex+1:end, end);
	end

	while (kNN_Uwant < 101)
	    wholetime = Match_Big_Data_Exp1(keepTargetData, queryComplete, ...
		subSeqLen, kNN_Uwant, whichDimToConsider, whichAlgoToConsider);   % working for the 2nd column here
	    
	    fprintf(fid, 'The Number of kNN is %d and time needed is : %d \n',kNN_Uwant, wholetime);
	    kNN_Uwant = kNN_Uwant + 5;
	end

	fclose(fid);
end

