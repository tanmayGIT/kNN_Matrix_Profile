function Journal_Server_2_Func_Protien_subSeq(datasetName, fileSavingPath)
  subSeqLen = 256;  % The length of the sub-sequence
  kNN_Uwant = 10;   % How many kNN do you want, just tell me

  load(datasetName);
  % fprintf('The data is correctly loaded and the size is : %d  %d \n', size(realData, 1), size(realData,2) )
  textFileNam = fileSavingPath;

  if (exist(textFileNam, 'file') == 0)
      disp('File does not exist, creating file.')
      fid = fopen( textFileNam, 'w' );
  else % if the file exist
      disp('File exists.');
      fid = fopen(textFileNam, 'wt' ); % create a new file each time
  end

  keepTargetData = realData(1:1000,:);
  % fprintf('The transformed data size is : %d  %d \n', size(keepTargetData, 1), size(keepTargetData,2) )
  clear realData;

  while (subSeqLen < 1001)
      wholetime = Match_Subseq_N_Query_SelfJoin_MultiCoreShort(keepTargetData, subSeqLen, kNN_Uwant);   % working for the 2nd column here

      fprintf(fid, 'The sub-sequ length is %d and time needed is : %d \n',subSeqLen, wholetime);
      subSeqLen = subSeqLen + 100;
  end
  fclose(fid);
end
