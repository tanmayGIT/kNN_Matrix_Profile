function Journal_Server_2_Func_Sheep_subSeq(datasetName, fileSavingPath)

  % Loading the sheep data
  load(datasetName)
  textFileNam = fileSavingPath;

  if (exist(textFileNam, 'file') == 0)
      disp('File does not exist, creating file.')
      fid = fopen( textFileNam, 'w' );
  else % if the file exist
      disp('File exists.');
      fid = fopen(textFileNam, 'wt' ); % create a new file each time
  end


  keepTargetData = keepAllData(1:1000,:);
  clear keepAllData;
  disp(size(keepTargetData));	
  subSeqLen = 256;  % The length of the sub-sequence
  whichDimToConsider = 1;
  kNN_Uwant = 10;   % How many kNN do you want, just tell me

  while (subSeqLen < 1001)
      wholetime = Match_Subseq_N_Query_SelfJoin_MultiCoreShort(keepTargetData, subSeqLen, kNN_Uwant);   % working for the 2nd column here

      fprintf(fid, 'The sub-sequ length is %d and time needed is : %d \n',subSeqLen, wholetime);
      subSeqLen = subSeqLen + 100;
  end
  fclose(fid);
end
