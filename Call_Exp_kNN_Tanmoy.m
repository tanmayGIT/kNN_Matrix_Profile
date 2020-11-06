% This function is used for the following :

% i) The variation of computational time with the increasing number of kNN similarity search
%     --> Journal_Exp_1_Fun_kNN();
% ii) The variation of computational time with the increasing sub-sequence length 
%     --> Journal_Server_2_Func_Sheep_subSeq();
% iii) The variation of computational time with the increasing time series length 
%     --> Journal_Exp_4_Fun_nFiles();
% iv) The variation of computational time with the increasing number of cores 
%     --> Journal_Exp_3_Fun_Cores();

% The objective of this function is to compute the computational time by
% using the proposed techique, mentioend in :

% Tanmoy Mondal, Reza Akbarinia,, and Florent Masseglia, "Matrix Profile 
% Based kNN Search over Large Time Series," submitted to: 
% "Elsevier Pattern Recognition Journal", 2020.
% https://sites.google.com/view/knnmatrixprofile/home


close
clc
clear

% rand_num = randi(1000);
Call_All_Code_In_Batch(799);
Call_All_Code_In_Batch(2492);


function Call_All_Code_In_Batch(rand_num)
	
	Journal_Server_2_Func_Sheep_subSeq('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/SheepDataFull.mat', "/local/tmp/tmondal/outputSelf/" + num2str(rand_num) + "_sheep_subSeq_SelfJoin.txt");
	fprintf('The sheep data processing is done to calculate variable sub-sequence exp for selfjoin case \n');
	Journal_Server_2_Func_Protien_subSeq('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/ProteinData.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "_protien_subSeq_SelfJoin.txt");
	fprintf('The protien data processing is done to calculate variable sub-sequence exp for selfjoin case \n');

	Journal_Exp_1_Fun_kNN('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "seismic_kNN_SelfJoin.txt");
	fprintf('The seismic_5000 data processing is done to calculate kNN for selfjoin case \n');
	Journal_Exp_1_Fun_kNN('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_kNN_SelfJoin.txt");
	fprintf('The randomWalk_5000 data processing is done to calculate kNN for selfjoin case \n');


	Journal_Exp_4_Fun_nFiles('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) +  "seismic_nFiles_SelfJoin.txt");
	fprintf('The seismic_5000 data processing is done to calculate nFiles exp for selfjoin case \n');
	Journal_Exp_4_Fun_nFiles('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_nFiles_SelfJoin.txt");
	fprintf('The randomWalk_5000 data processing is done to calculate nFiles exp for selfjoin case \n');


	Journal_Exp_3_Fun_Cores('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) +  "seismic_nCores_SelfJoin.txt");
	fprintf('The seismic_5000 data processing is done to calculate number of cores exp for selfjoin case \n');
	Journal_Exp_3_Fun_Cores('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_nCores_SelfJoin.txt");
	fprintf('The randomWalk_5000 data processing is done to calculate number of cores exp for selfjoin case \n');
end
