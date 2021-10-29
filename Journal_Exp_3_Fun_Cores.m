function Journal_Exp_3_Fun_Cores(datasetName, savingPathFileName)

 % this is for self join where the number of cores will be varying 

  subSeqLen = 256;  % The length of the sub-sequence
  whichDimToConsider = 1;
  kNN_Uwant = 10;   % How many kNN do you want, just tell me

  % Loading the sheep data
  load(datasetName);
  keepAllData = seeMe1(1:1000,:); % take only 10000 time series
  clear seeMe1;

  textFileNam = savingPathFileName;

  if (exist(textFileNam, 'file') == 0)
      disp('File does not exist, creating file.')
      fid = fopen( textFileNam, 'w' );
  else % if the file exist
      disp('File exists.');
      fid = fopen(textFileNam, 'wt' ); % create a new file each time
  end

  keepTargetData = keepAllData;
  clear keepAllData;
  nCores = 2;	
  while (nCores < 37)
      wholetime = Match_SelfJoin_VariableCore(keepTargetData, subSeqLen, kNN_Uwant, nCores); % working for the 2nd column here

      fprintf(fid, 'The Number of cores is %d and time needed is : %d \n',nCores, wholetime);
      nCores = nCores + 5;
  end

  fclose(fid);
end
