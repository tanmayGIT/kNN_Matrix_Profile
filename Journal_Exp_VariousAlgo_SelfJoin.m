
% The variation of computational time with the increasing number of kNN similarity search. 
% The computational time analysis for three different proposed approaches 
% i.e. sort based (refer to Section 4.2.1),  max based (refer to Section
% 4.2.2 ) and heap-max based (refer to Section 4.2.3 ) for kNN similarity search.
%
% Journal_Exp_VariousAlgo_SelfJoin(datasetName, savingPathFileName, whichAlgoToConsider)
%
% Input:
%     datasetName: input dataset (matrix: .mat file)
%     savingPathFileName: the text file path where the computational time
%     are saved
%     whichAlgoToConsider: the name of the algorithm to run : "sort",
%     "max", "heap"
%
% Tanmoy Mondal, Reza Akbarinia,, and Florent Masseglia, "Matrix Profile 
% Based kNN Search over Large Time Series," submitted to: 
% "Elsevier Pattern Recognition Journal", 2020.
% https://sites.google.com/view/knnmatrixprofile/home

function Journal_Exp_VariousAlgo_SelfJoin(datasetName, savingPathFileName, whichAlgoToConsider)

  subSeqLen = 256;  % The length of the sub-sequence
  kNN_Uwant = 5;   % How many kNN do you want, just tell me

  % Loading the sheep data
  load(datasetName);
  keepAllData = seeMe1(1:100,:); % take only 100 time series
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

  while (kNN_Uwant < 65)
      wholetime = Match_SelfJoin_VariousAlgo(keepTargetData, subSeqLen, kNN_Uwant, whichAlgoToConsider, 36); 

      fprintf(fid, 'The Number of kNN is %d and time needed is : %d \n',kNN_Uwant, wholetime);
      kNN_Uwant = kNN_Uwant + 5;
  end

  fclose(fid);
end
