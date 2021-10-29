clear;
clc;
close all;

dataFolderPath = '/home/tmondal/Documents/All_Safrans/Safran_Project/Matrix_Profile_Server/UCR_Matrix_Profile_Results/';
newSavePath = '/home/tmondal/Documents/All_Safrans/Safran_Project/Matrix_Profile_Server/UCR_Selected/';
subFolders = GetSubDirsFirstLevelOnly(dataFolderPath);



%% This part of the code is for renaming the files according to the name of the folder
% for iDir = 1:1:length(subFolders)
%     datasetName = subFolders{1, iDir};
%     datasetPath_Full = strcat(dataFolderPath, datasetName, '/');
%     
%     dirData = dir(datasetPath_Full);      % Get the data for the current directory
%     dirIndex = [dirData.isdir];  % Find the index for directories
%     fileList = {dirData(~dirIndex).name}';
%     
%     % Loop through each
%     for id = 1:length(fileList)
%         % Get the file name (minus the extension)
%         [~, f, ext] = fileparts(fileList{id});
%         new_file_name = strcat(datasetName, '_', f, ext);
%         
%         movefile( fullfile(datasetPath_Full, fileList{id}), fullfile(datasetPath_Full, new_file_name) );
%     end
% end


%% This part of the code is for picking selected images from each dataset folder and to put it inside a single folder
for iDir = 1:1:length(subFolders)
    datasetName = subFolders{1, iDir};
    datasetPath_Full = strcat(dataFolderPath, datasetName, '/');
    dirData = dir(datasetPath_Full);      % Get the data for the current directory
    dirIndex = [dirData.isdir];  % Find the index for directories
    fileList = {dirData(~dirIndex).name}';

    switch datasetName 
        case 'ACSF1'
            goodIndexes = [26, 32, 35, 79, 132, 136, 192];
        case 'Adiac'
            goodIndexes = [43, 60, 109, 117, 154, 199, 202, 274, 275, 415, 428, 458, 589, 636, 696, 717, 726];
        case 'ArrowHead'
            goodIndexes = [32, 79, 85, 88, 92, 175, 208];
        case 'BME'
            goodIndexes = [5, 20, 30, 48, 57, 62, 70, 78, 87, 105, 112, 117, 148, 149, 177];
        case 'Beef'
            goodIndexes = [2, 6 ,11, 13, 19, 21, 25, 26, 28, 36, 40, 46, 50];
        case 'BeetleFly'
            goodIndexes = [3, 10, 11, 13, 16, 20, 21, 23, 24, 29];
        case 'BirdChicken'
            goodIndexes = [3, 8, 9, 12, 13, 14, 15, 18, 27, 30];
        case 'CBF'
            goodIndexes = [186, 921];
        case 'Car'
            goodIndexes = [];
        case 'ChlorineConcentration'
            goodIndexes = 3045;
        case 'CinCECGTorso'
            goodIndexes = [480, 1108, 1173];
        case 'Coffee'
            goodIndexes = 5;
        case 'CricketX'
            goodIndexes = [58, 77, 85, 99, 154, 166, 202, 247, 248, 436, 549];
        case 'CricketY'
            goodIndexes = [120, 128 ,549, 643];
        case 'CricketZ'
            goodIndexes = [20, 36, 138, 175, 347, 472, 500, 509, 567, 742, 747];
        case 'DiatomSizeReduction'
            goodIndexes = [14, 35, 109, 188, 199, 247, 290];
        case 'DodgerLoopDay'
            goodIndexes = [16, 29, 31];
        case 'DodgerLoopWeekend'
            goodIndexes = [4, 39, 45, 56, 72, 78, 93, 94, 99, 154, 136, 133, 116, 107, 106, 105];
        case 'ECG5000'
            goodIndexes = 1346;
        case 'ECGFiveDays'
            goodIndexes = [109, 132, 244, 588, 596, 675, 764, 781, 793, 808];
        case 'EOGHorizontalSignal'
            goodIndexes = [];
        case 'EOGVerticalSignal'
            goodIndexes = [];
        case 'Earthquakes'
            goodIndexes = [41, 55, 64, 65, 88, 102, 113, 119, 167, 180, 186, 187, 198, 277, 326, 362];
        case 'EthanolLevel'
            goodIndexes = [];
        case 'FaceAll'
            goodIndexes = [129, 130, 262, 641, 842, 880];
        case 'FaceFour'
            goodIndexes = [2, 16, 17, 26, 27, 31, 35, 40, 51, 63, 66, 90, 110];
        case 'FacesUCR'
            goodIndexes = [439, 1803, 1974, 439];
        case 'FiftyWords'
            goodIndexes = [45, 83, 162, 200, 219, 297, 396, 418, 443, 888, 873, 831];
        case 'Fish'
            goodIndexes = 13;
        case 'FordA'
            goodIndexes = [2, 2016, 2040];
        case 'FordB'
            goodIndexes = [906, 994, 1292, 1389, 1790, 1913, 2599, 3005, 3394, 3560, 3690, 3833];
        case 'FreezerRegularTrain'
            goodIndexes = [2705, 2956, 109];
        case 'FreezerSmallTrain'
            goodIndexes = [576, 1324, 1459, 1479, 1727];
        case 'Fungi'
            goodIndexes = [12, 43, 66, 91, 168, 187];
        case 'GunPoint'
            goodIndexes = 130;
        case 'GunPointAgeSpan'
            goodIndexes = [80, 215, 323];
        case 'GunPointMaleVersusFemale'
            goodIndexes = [60, 387];
        case 'GunPointOldVersusYoung'
            goodIndexes = [56, 64, 185, 257, 406];
        case 'Ham'
            goodIndexes = [12, 42, 73, 81, 98, 102,  110,  114, 117, 119, 121, 126, 211, 212];
        case 'HandOutlines'
            goodIndexes = [];
        case 'Haptics'
            goodIndexes = [];
        case 'Herring'
            goodIndexes = [4, 7, 27, 48, 57, 72];
        case 'HouseTwenty'
            goodIndexes = [];
        case 'InlineSkate'
            goodIndexes = [];
        case 'InsectEPGRegularTrain'
            goodIndexes = [62, 101, 241, 295, 304];
        case 'InsectEPGSmallTrain'
            goodIndexes = [];
        case 'InsectWingbeatSound'
            goodIndexes = [];
        case 'LargeKitchenAppliances'
            goodIndexes = [];
        case 'Lightning7'
            goodIndexes = [5, 12, 14, 82, 102];
        case 'Mallat'
            goodIndexes = [117, 475, 517, 1562, 1875, 1977, 2282, 2385];
        case 'Meat'
            goodIndexes = [9, 25, 31, 37, 42, 44, 47, 49, 56, 70, 85, 107];
        case 'MixedShapesRegularTrain'
            goodIndexes = 103;
        case 'MixedShapesSmallTrain'
            goodIndexes = [598, 935, 968, 1172];
        case 'NonInvasiveFetalECGThorax1'
            goodIndexes = [828, 1098, 1334];
        case 'NonInvasiveFetalECGThorax2'
            goodIndexes = [];
        case 'OSULeaf'
            goodIndexes = [8, 67, 80, 81, 150, 211, 233, 235, 281, 291, 382, 424];
        case 'OliveOil'
            goodIndexes = [6, 17];
        case 'Phoneme'
            goodIndexes = [];
        case 'PigAirwayPressure'
            goodIndexes = [];
        case 'PigArtPressure'
            goodIndexes = 83;
        case 'PigCVP'
            goodIndexes = 147;
        case 'Plane'
            goodIndexes = [12, 17, 36, 61, 69, 73, 84, 88, 122, 143, 154, 187, 200];
        case 'PowerCons'
            goodIndexes = [40, 118, 206, 244, 259, 323, 338];
        case 'Rock'
            goodIndexes = [7, 11, 53, 55];
        case 'SemgHandGenderCh2'
            goodIndexes = [];
        case 'SemgHandMovementCh2'
            goodIndexes = [142, 212, 352];
        case 'SemgHandSubjectCh2'
            goodIndexes = [399, 613];
        case 'ShapeletSim'
            goodIndexes = [9, 14, 18, 20, 22, 23, 24, 37, 42, 50, 54, 85, 91, 109, 130, 141, 147, 149, 163, 164, 169, 191];
        case 'ShapesAll'
            goodIndexes = [5, 93, 185, 188, 221, 376, 439, 767, 781, 785, 884, 981];
        case 'StarLightCurves'
            goodIndexes = [5267, 7589, 8508];
        case 'Strawberry'
            goodIndexes = [7, 67, 132, 201, 223, 249, 387, 439, 511, 516, 627, 734, 793, 809, 931, 938, 951];
        case 'SwedishLeaf'
            goodIndexes = [8, 58, 123, 185, 204, 235, 292, 403, 609, 616, 742, 789, 979, 986, 1024];
        case 'Symbols'
            goodIndexes = [289, 294, 403 460, 675, 806, 823];
        case 'ToeSegmentation1'
            goodIndexes = [31, 62, 88, 114, 138, 207, 211, 267];
        case 'ToeSegmentation2'
            goodIndexes = [];
        case 'Trace'
            goodIndexes = [];
        case 'TwoPatterns'
            goodIndexes = [304, 338, 1583, 1725, 2388, 3202, 3332, 3370 ,3536, 4069, 4357, 4399, 4798];
        case 'UMD'
            goodIndexes = [4, 12, 14, 25, 29, 30, 36, 44, 51, 88, 97, 98, 119, 123];
        case 'UWaveGestureLibraryAll'
            goodIndexes = [217, 1456];
        case 'UWaveGestureLibraryX'
            goodIndexes = [];
        case 'UWaveGestureLibraryZ'
            goodIndexes = [137, 154, 252, 759, 865, 1293, 1323, 1564, 1569, 1640, 1721, 1790, 1948, 2079, 2488, 2499, 2658, 4443];
        case 'Wafer'
            goodIndexes = 133;
        case 'Wine'
            goodIndexes = [7, 24, 27, 29, 34, 37, 39, 43, 48, 52, 58];
        case 'WordSynonyms'
            goodIndexes = [21, 36, 57, 134, 184, 224, 254, 293, 378, 420, 542, 707, 787, 814];
        case 'Worms'
            goodIndexes = [];
        case 'WormsTwoClass'
            goodIndexes = 73;
        case 'Yoga'
            goodIndexes = [122, 201, 561, 599, 1777, 1994, 2490];
    end
    
    if(~isempty(goodIndexes))
        for id = 1:length(goodIndexes)
            formFileName_1 = strcat(datasetName, '_',num2str(goodIndexes(id)), '_' , '1', '_', '2', '.png');
            formFileName_2 = strcat(datasetName, '_',num2str(goodIndexes(id)), '_' , '1', '_', '4','.png');

            fullPathFormFileName_1 = strcat(datasetPath_Full, formFileName_1);
            fullPathFormFileName_2 = strcat(datasetPath_Full, formFileName_2);

            fullNewFilePath_1 = strcat(newSavePath, formFileName_1);
            fullNewFilePath_2 = strcat(newSavePath, formFileName_2);

            copyfile(fullPathFormFileName_1,fullNewFilePath_1);
            copyfile(fullPathFormFileName_2,fullNewFilePath_2);

        end
    end
    
end


function [subDirsNames] = GetSubDirsFirstLevelOnly(parentDir)
% Get a list of all files and folders in this folder.
files = dir(parentDir);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subDirs = files(dirFlags);
subDirsNames = cell(1, numel(subDirs) - 2);
for i=3:numel(subDirs)
    subDirsNames{i-2} = subDirs(i).name;
end
end