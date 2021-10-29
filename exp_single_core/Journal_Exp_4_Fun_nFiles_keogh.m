% This experiment is for varying kNN
function Journal_Exp_4_Fun_nFiles_keogh(datasetName, textFileSavingName)

  subSeqLen = 256;  % The length of the sub-sequence
  whichDimToConsider = 1;
  kNN_Uwant = 10;   % How many kNN do you want, just tell me

  % Loading the sheep data
  load(datasetName);
  keepAllData = seeMe1(1:100,:); % take only 16000 time series ;
  clear seeMe1;


  textFileNam = textFileSavingName;

  if (exist(textFileNam, 'file') == 0)
      disp('File does not exist, creating file.')
      fid = fopen( textFileNam, 'w' );
  else % if the file exist
      disp('File exists.');
      fid = fopen(textFileNam, 'wt' ); % create a new file each time
  end

  keepTargetData = keepAllData;
  clear keepAllData;


  takenRws = 10;
  while (takenRws < 101)

      % fprintf('The number of rows to be considered %d and the number of rows present in the matrix %d', takenRws, size(keepTargetData,1));
      if(takenRws > size(keepTargetData,1))
         takenRws = size(keepTargetData,1);
      end
      keepTargetDataCut = keepTargetData(1:takenRws,:);
      wholetime = Match_SelfJoin_Single_Core_Keogh(keepTargetDataCut, subSeqLen, kNN_Uwant);   % working for the 2nd column here

      fprintf(fid, 'The Number of rows is %d and time needed is : %d \n',takenRws, wholetime);
      takenRws = takenRws + 20;
      clear keepTargetDataCut
  end
  fclose(fid);

end
