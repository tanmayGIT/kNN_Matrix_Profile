function Journal_Exp_1_Fun_kNN_keogh(datasetName, savingPathFileName)

  subSeqLen = 256;  % The length of the sub-sequence
  whichDimToConsider = 1;
  kNN_Uwant = 10;   % How many kNN do you want, just tell me

  % Loading the sheep data
  load(datasetName);
  keepAllData = seeMe1(1:100,:); % take only 10000 time series
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

  while (kNN_Uwant < 101)
      wholetime = Match_SelfJoin_Single_Core_Keogh(keepTargetData, subSeqLen, kNN_Uwant); % working for the 2nd column here

      fprintf(fid, 'The Number of kNN is %d and time needed is : %d \n',kNN_Uwant, wholetime);
      kNN_Uwant = kNN_Uwant + 20;
  end

  fclose(fid);
end
