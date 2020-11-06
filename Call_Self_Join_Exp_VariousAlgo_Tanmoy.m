
% The variation of computational time with the increasing number of kNN similarity search. 
% The computational time analysis for three different proposed approaches 
% i.e. sort based (refer to Section 4.2.1),  max based (refer to Section
% 4.2.2 ) and heap-max based (refer to Section 4.2.3 ) for kNN similarity search.


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
	

	Journal_Exp_VariousAlgo_SelfJoin('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) +  "_seismic_kNN_SortBased_Self_Join.txt", "sortBased");
	fprintf('The seismic_5000 data processing is done to calculate kNN processing experiment for self join case using sort based algorithm  \n');

	Journal_Exp_VariousAlgo_SelfJoin('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) +  "_seismic_kNN_MaxBased_Self_Join.txt", "maxBased");
	fprintf('The seismic_5000 data processing is done to calculate kNN processing experiment for self join case using max based algorithm  \n');

	Journal_Exp_VariousAlgo_SelfJoin('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) +  "_seismic_kNN_HeapBased_Self_Join.txt", "heapBased");
	fprintf('The seismic_5000 data processing is done to calculate kNN processing experiment for self join case using sort based algorithm  \n');



	Journal_Exp_VariousAlgo_SelfJoin('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_kNN_SortBased_Self_Join.txt", "sortBased");
	fprintf('The randomWalk_5000 data processing is done to calculate kNN processing experiment for self join case using sort based algorithm  \n');

	Journal_Exp_VariousAlgo_SelfJoin('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_kNN_MaxBased_Self_Join.txt", "maxBased");
	fprintf('The randomWalk_5000 data processing is done to calculate kNN processing experiment for self join case using sort based algorithm  \n');

	Journal_Exp_VariousAlgo_SelfJoin('/local/tmp/tmondal/inputSelf/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/outputSelf/" +  num2str(rand_num) + "randomWalk_kNN_HeapBased_Self_Join.txt", "heapBased");
	fprintf('The randomWalk_5000 data processing is done to calculate kNN processing experiment for self join case using sort based algorithm  \n');
end
