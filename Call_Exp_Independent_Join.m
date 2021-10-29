close
clc
clear

Call_All_Code_In_Batch(12);
Call_All_Code_In_Batch(22);



function Call_All_Code_In_Batch(rand_num)
    Journal_Exp_1_kNN('/local/tmp/tmondal/input/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/output/" + num2str(rand_num) + "_seismic_kNN_SortBased_Independent_Join.txt", "sortBased");
    fprintf('The seismic data processing is done to calculate kNN processing experiment for independent join case using sort based algorithm \n');

    Journal_Exp_1_kNN('/local/tmp/tmondal/input/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/output/" + num2str(rand_num) + "_seismic_kNN_MaxBased_Independent_Join.txt", "maxBased");
    fprintf('The seismic data processing is done to calculate kNN processing experiment for independent join case using max based algorithm \n');

    Journal_Exp_1_kNN('/local/tmp/tmondal/input/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/output/" + num2str(rand_num) + "_seismic_kNN_HeapBased_Independent_Join.txt", "heapBased");
    fprintf('The seismic data processing is done to calculate kNN processing experiment for independent join case using heap based algorithm \n');




    Journal_Exp_1_kNN('/local/tmp/tmondal/input/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/output/" + num2str(rand_num) + "_randomWalk_kNN_SortBased_Independent_Join.txt", "sortBased");
    fprintf('The random walk data processing is done to calculate kNN processing experiment for independent join case using sort based algorithm \n');

    Journal_Exp_1_kNN('/local/tmp/tmondal/input/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/output/" + num2str(rand_num) + "_randomWalk_kNN_MaxBased_Independent_Join.txt", "maxBased");
    fprintf('The random walk data processing is done to calculate kNN processing experiment for independent join case using max based algorithm \n');

    Journal_Exp_1_kNN('/local/tmp/tmondal/input/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/output/" + num2str(rand_num) + "_randomWalk_kNN_HeapBased_Independent_Join.txt", "heapBased");
    fprintf('The random walk data processing is done to calculate kNN processing experiment for independent join case using heap based algorithm \n');




    Journal_Exp_2_Subseq_Protien('/local/tmp/tmondal/input/Matrix_Profile_Server/ProteinData.mat', "/local/tmp/tmondal/output/" +  num2str(rand_num) + "_protien_subSeq_Independent_Join.txt");
    fprintf('The protien data processing is done to calculate variable sub-sequence exp for independent join case \n');

    Journal_Exp_2_Subseq_Sheep('/local/tmp/tmondal/input/Matrix_Profile_Server/SheepDataFull.mat', "/local/tmp/tmondal/output/" +  num2str(rand_num) + "_sheep_subSeq_Independent_Join.txt");
    fprintf('The sheep data processing is done to calculate variable sub-sequence exp for independent join case \n');




    Journal_Exp_3_Cores('/local/tmp/tmondal/input/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/output/" +  num2str(rand_num) + "seismic_cores_Independent_Join.txt");
    fprintf('The seismic_5000 data processing is done to calculate kNN for independent join case \n');

    Journal_Exp_3_Cores('/local/tmp/tmondal/input/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/output/" +  num2str(rand_num) + "randomWalk_cores_Independent_Join.txt");
    fprintf('The randomWalk_5000 data processing is done to calculate kNN for independent join case \n');




    Journal_Exp_4_NumOfFiles('/local/tmp/tmondal/input/Matrix_Profile_Server/seismic_50000.mat', "/local/tmp/tmondal/output/" +  num2str(rand_num) +  "seismic_nFiles_Independent_Join.txt");
    fprintf('The seismic_5000 data processing is done to calculate nFiles exp for independent join case \n');

    Journal_Exp_4_NumOfFiles('/local/tmp/tmondal/input/Matrix_Profile_Server/randomWalk_50000.mat', "/local/tmp/tmondal/output/" +  num2str(rand_num) + "randomWalk_nFiles_Independent_Join.txt");
    fprintf('The randomWalk_5000 data processing is done to calculate nFiles exp for independent join case \n');
end
