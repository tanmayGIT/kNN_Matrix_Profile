close
clc
clear

% rand_num = randi(1000);
Call_All_Code_In_Batch(799);
Call_All_Code_In_Batch(2492);


function Call_All_Code_In_Batch(rand_num)
	
	Journal_Server_2_Func_Sheep_subSeq('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/SheepDataFull.mat', "/local/tmp/tmondal/outputSelf/" + num2str(rand_num) + "_sheep_subSeq_SelfJoin_single_core.txt");
	fprintf('The sheep data processing is done to calculate variable sub-sequence exp for selfjoin case\n');
	Journal_Server_2_Func_Protien_subSeq('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/ProteinData.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "_protien_subSeq_SelfJoin_single_core.txt");
	fprintf('The protien data processing is done to calculate variable sub-sequence exp for selfjoin case \n');

	Journal_Exp_1_Fun_kNN('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "seismic_kNN_SelfJoin_single_core.txt");
	fprintf('The seismic_5000 data processing is done to calculate kNN for selfjoin case \n');
	Journal_Exp_1_Fun_kNN('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_kNN_SelfJoin_single_core.txt");
	fprintf('The randomWalk_5000 data processing is done to calculate kNN for selfjoin case \n');


	Journal_Exp_4_Fun_nFiles('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) +  "seismic_nFiles_SelfJoin_single_core.txt");
	fprintf('The seismic_5000 data processing is done to calculate nFiles exp for selfjoin case \n');
	Journal_Exp_4_Fun_nFiles('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_nFiles_SelfJoin_single_core.txt");
	fprintf('The randomWalk_5000 data processing is done to calculate nFiles exp for selfjoin case \n');


% 	Journal_Exp_3_Fun_Cores('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) +  "seismic_nCores_SelfJoin.txt");
% 	fprintf('The seismic_5000 data processing is done to calculate number of cores exp for selfjoin case \n');
% 	Journal_Exp_3_Fun_Cores('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_nCores_SelfJoin.txt");
% 	fprintf('The randomWalk_5000 data processing is done to calculate number of cores exp for selfjoin case \n');
end