function Journal_Exp_2_Subseq_Sheep(dataPath, textfileToSavePath)
	% This experiment is for varying kNN

	% Saving the sheep data
	% txtFilePath = '/home/mondal/Documents/Dataset/Time_Series/Donnees/allBrebisDataA.txt';
	% fid = fopen(txtFilePath);
	% tline = fgetl(fid);
	% keepAllData = zeros(1,1);
	% 
	% % HERE I AM NOT TAKING ALL THE FILES AS I AM IGNORING SOME SPECIFIC LABELS WHICH ARE NOT NEEDED FOR MY CASE 
	% cntLn = 1;
	% while ischar(tline)
	%     splitName = strsplit(tline,' '); % split the line
	% 
	%     for ii = 1:1:(length(splitName)-1)
	%         keepAllData(cntLn,ii) = str2double(splitName{ii});
	%     end
	%     cntLn = cntLn +1;
	%     
	%     tline = fgetl(fid);
	% end
	% save('SheepDataFull.mat','keepAllData');



	% Loading the sheep data
	load(dataPath) % load('SheepDataFull.mat')

	queryIndex = 307;
	queryComplete = keepAllData(queryIndex, 1:end-1);

	noOfZerosInQuery = sum(queryComplete(:)==0);
	if(noOfZerosInQuery > (length(queryComplete) /2)) % if the number of zeros is more than half the length of query
	    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
	end


	textFileNam = textfileToSavePath; % 'SheepDataFull_VariableSubSeqLen.txt';

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


	subSeqLen = 30;  % The length of the sub-sequence
	whichDimToConsider = 1;
	kNN_Uwant = 10;   % How many kNN do you want, just tell me

	while (subSeqLen < ( round((length(queryComplete)*90)/100) ))
	    wholetime = Match_Big_Data(keepTargetData, queryComplete, ...
		subSeqLen, kNN_Uwant, whichDimToConsider);   % working for the 2nd column here
	    
	    fprintf(fid, 'The sub-sequ length is %d and time needed is : %d \n',subSeqLen, wholetime);
	    subSeqLen = subSeqLen + 20;
	end
	fclose(fid);
end
