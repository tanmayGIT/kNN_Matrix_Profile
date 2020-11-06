function Journal_Exp_2_Subseq_Protien(dataPath, textfileToSavePath)

	subSeqLen = 30;  % The length of the sub-sequence
	whichDimToConsider = 1;
	kNN_Uwant = 10;   % How many kNN do you want, just tell me

	load(dataPath) % load ('ProteinData.mat');

	queryIndex = 1165;      % randsample(size(realData,1),1);

	queryComplete = realData(queryIndex, 1:end);

	textFileNam = textfileToSavePath; % 'Protien_VariableSubSeqLen.txt';

	if (exist(textFileNam, 'file') == 0)
	    disp('File does not exist, creating file.')
	    fid = fopen( textFileNam, 'w' );
	else % if the file exist
	    disp('File exists.');
	    fid = fopen(textFileNam, 'wt' ); % create a new file each time 
	end

	noOfZerosInQuery = sum(queryComplete(:)==0);
	if(noOfZerosInQuery > (length(queryComplete) /2)) % if the number of zeros is more than half the length of query
	    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
	end

	keepTargetData = zeros((size(realData,1)-1), length(queryComplete));

	if(queryIndex == 1)
	    keepTargetData(:,:) = realData(2:end, 1:end-1);
	elseif(queryIndex == size(realData,1))
	    keepTargetData(:,:) = realData(1:end-1, 1:end-1);
	else
	    keepTargetData(1:queryIndex-1,:) = realData(1:queryIndex-1, 1:end);
	    keepTargetData(queryIndex:end,:) = realData(queryIndex+1:end, 1:end);
	end



	while (subSeqLen < ( round((length(queryComplete)*90)/100) ))
	    wholetime = Match_Big_Data(keepTargetData, queryComplete, ...
		subSeqLen, kNN_Uwant, whichDimToConsider);   % working for the 2nd column here
	    
	    fprintf(fid, 'The sub-sequ length is %d and time needed is : %d \n',subSeqLen, wholetime);
	    subSeqLen = subSeqLen + 20;
	end
	fclose(fid);
end
