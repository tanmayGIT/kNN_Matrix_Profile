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
            goodIndexes = [136];
        case 'Adiac'
            goodIndexes = [43, 109, 199, 428, 589, 726];
        case 'BME'
            goodIndexes = [30, 48, 87, 105, 117, 148, 149, 177];
        case 'Beef'
            goodIndexes = [6 , 28, 36, 40, 46, 50];
        case 'BeetleFly'
            goodIndexes = [21];
        case 'BirdChicken'
            goodIndexes = [9, 12, 18];
        case 'CinCECGTorso'
            goodIndexes = [480];
        case 'CricketX'
            goodIndexes = [436];
        case 'CricketZ'
            goodIndexes = [509];
        case 'DiatomSizeReduction'
            goodIndexes = [247];
        case 'DodgerLoopDay'
            goodIndexes = [31];
        case 'DodgerLoopWeekend'
            goodIndexes = [4, 56, 78, 93, 94, 99, 133, 107];
        case 'Earthquakes'
            goodIndexes = [41, 88, 102, 113, 119, 186];
        case 'FaceAll'
            goodIndexes = [842, 880];
        case 'FaceFour'
            goodIndexes = [16];
        case 'FacesUCR'
            goodIndexes = [439, 1974];
        case 'FiftyWords'
            goodIndexes = [888];
        case 'FordA'
            goodIndexes = [2016, 2040];
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