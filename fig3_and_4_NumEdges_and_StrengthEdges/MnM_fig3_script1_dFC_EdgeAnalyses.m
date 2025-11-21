%% Set-up
close all; clear; clc

xcpddir = '/Users/josephchen/data/MnM/MnM-XCPD-NoGSRNoCensor';
outputdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/fig3_and_4_NumEdges_and_StrengthEdges';
datadir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data';
addpath('/Users/josephchen/Library/CloudStorage/GoogleDrive-nzjosephchen@gmail.com/My Drive/nzjosephchen/software/2019_03_03_BCT')

% Read personality data
PersonalityFile = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/MnM_PersonalityInfo.csv';
PersonalityData = readtable(PersonalityFile);

%% Read Atlas data
AtlasFile = readtable('/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/atlas-4S156Parcels_dseg.tsv', ...
    'FileType', 'text', 'Delimiter', '\t');
NumROIs = size(AtlasFile,1);

% Reordering the ROIs into their 7 networks
[~, ROIold2new] = sort(AtlasFile.network_label);
AtlasFile2 = AtlasFile(ROIold2new,:);
ROInew2old = AtlasFile2.index;

CEN_ROIs = AtlasFile.index(strcmp(AtlasFile.network_label, 'Cont'));
DMN_ROIs = AtlasFile.index(strcmp(AtlasFile.network_label, 'SalVentAttn'));
SN_ROIs = AtlasFile.index(strcmp(AtlasFile.network_label, 'Default'));
CEN_ROIs_num = length(CEN_ROIs);
DMN_ROIs_num = length(DMN_ROIs);
SN_ROIs_num = length(SN_ROIs);
NotCEN_ROIs = setdiff(1:NumROIs, CEN_ROIs);
NotDMN_ROIs = setdiff(1:NumROIs, DMN_ROIs);
NotSN_ROIs = setdiff(1:NumROIs, SN_ROIs);

%% Analyze FC to DMN/CEN/SN

AllSubjDir = dir([datadir '/MnM_dFC/r_0.5/k_1.6/sub-*']);
AllSubj = {AllSubjDir.name};
AllSubj(find(strcmp(AllSubj, 'sub-2004'))) = []; %remove sub-2004 because no extraversion rating
AllSubj(find(strcmp(AllSubj, 'sub-2003'))) = []; %remove sub-2003
AllSubj(find(strcmp(AllSubj, 'sub-2020'))) = []; %remove sub-2020
NumSubj = length(AllSubj);
AllTasks = {'BreathNVHA' 'BreathNVLA' 'Meditation'};

OutTable = [];
OutTable.ID = strings(1,4000);
OutTable.Task = strings(1,4000);
OutTable.Extraversion = strings(1,4000);
OutTable.Extravert = strings(1,4000);
OutTable.AllConn_CEN = strings(1,4000);
OutTable.AllConn_SN = strings(1,4000);
OutTable.AllConn_DMN = strings(1,4000);
OutTable.Intra_CEN = strings(1,4000);
OutTable.Intra_SN = strings(1,4000);
OutTable.Intra_DMN = strings(1,4000);
OutTable.Extra_CEN = strings(1,4000);
OutTable.Extra_SN = strings(1,4000);
OutTable.Extra_DMN = strings(1,4000);
OutTable.CEN_SN = strings(1,4000);
OutTable.SN_DMN = strings(1,4000);
OutTable.DMN_CEN = strings(1,4000);
OutTable.Intra_CEN_z = strings(1,4000);
OutTable.Intra_SN_z = strings(1,4000);
OutTable.Intra_DMN_z = strings(1,4000);
OutTable.Extra_CEN_z = strings(1,4000);
OutTable.Extra_SN_z = strings(1,4000);
OutTable.Extra_DMN_z = strings(1,4000);
OutTable.CEN_SN_z = strings(1,4000);
OutTable.SN_DMN_z = strings(1,4000);
OutTable.DMN_CEN_z = strings(1,4000);

count = 1;
for iSubj = 1:NumSubj
    tempSubj = AllSubj{iSubj};

    % get extraversion + extravert scores
    [~,loc] = ismember(tempSubj, PersonalityData.ID);
    tempExtraversionRaw = PersonalityData.Personality_ExtraversionRaw(loc);
    if tempExtraversionRaw < 24; tempExtrovert = 0; end
    if tempExtraversionRaw > 24; tempExtrovert = 1; end
    if tempExtraversionRaw == 24; tempExtrovert = nan; end

    for iTask = 1:3
        tempTask = AllTasks{iTask};

        % read in the original XCPD output to check for any NAN columns
        warning off
        tempOrigTimeSeries = readtable([xcpddir '/' tempSubj '/func/' tempSubj '_task-' tempTask ...
            '_space-MNI152NLin2009cAsym_atlas-4S156Parcels_timeseries.tsv'], ...
            'FileType', 'text', 'Delimiter', '\t');
        warning on
        NANChecker = false(1,156);
        if size(tempOrigTimeSeries,2) ~= 156; error('Not 156 ROIs. exiting...'); end
        for iColumn = 1:156
            tempValue = tempOrigTimeSeries{1,iColumn};
            if ~isnumeric(tempValue)
                NANChecker(iColumn) = true;
            end
        end
        NANindexes = find(NANChecker);
        
        % now read each FC matrix as delineated by dynamic FC analyses
        tempFCMatrixDir = dir([datadir '/MnM_dFC/r_0.5/k_1.6/' tempSubj '/' tempTask '/FC_*']);
        tempFCMatrixFiles = {tempFCMatrixDir.name};
        tempNumFCMatrix = length(tempFCMatrixFiles);

        for iMatrix = 1:tempNumFCMatrix
            tempMatrixFile = [datadir '/MnM_dFC/r_0.5/k_1.6/' tempSubj '/' tempTask '/' tempFCMatrixFiles{iMatrix}];
            tempMatrixTable = readtable(tempMatrixFile);
            tempMatrixArray = table2array(tempMatrixTable);
            %tempMatrixArray_TR = triu(tempMatrixArray);
            %tempMatrixArray_TR(tempMatrixArray_TR == 0) = nan;

            % now apply threshold to connectivity matrix
            tempSortValues = sort(abs(tempMatrixArray(:)), 'descend');
            tempSortValues = tempSortValues(~isnan(tempSortValues));
            tempThreshold = tempSortValues(floor(0.2*length(tempSortValues)));
            tempAdjacencyMatrix = abs(tempMatrixArray) >= tempThreshold;

            % now use the threshold on the original matrix
            tempMatrixArray = tempMatrixArray .* tempAdjacencyMatrix;

            % now make bottom left of AdjacencyMatrix zero, so you don't double count
            %tempAdjacencyMatrix_TR = triu(tempAdjacencyMatrix);

            % check if matrix is 156x156. If not, add in 0 columns in correct location
            if ~isequal(size(tempAdjacencyMatrix), [156, 156])
                
                %sanity check: if NANindexes is EMPTY, then we've got mismatch
                if sum(NANChecker) == 0; error('mismatch. exiting...'); end
                
                %add column of zeros in correct location
                for iNANCol = 1:size(NANindexes,2)
                    tempNAN = NANindexes(iNANCol);
                    tempSize = size(tempAdjacencyMatrix,1);
                    % add column of zeros at tempNAN's column
                    tempAdjacencyMatrix = [tempAdjacencyMatrix(:,1:tempNAN-1), zeros(tempSize,1), tempAdjacencyMatrix(:,tempNAN:end)];
                    tempMatrixArray = [tempMatrixArray(:,1:tempNAN-1), nan(tempSize,1), tempMatrixArray(:,tempNAN:end)];
                    % add row of zeros at tempNAN's row
                    tempAdjacencyMatrix = [tempAdjacencyMatrix(1:tempNAN-1,:); zeros(1,tempSize+1); tempAdjacencyMatrix(tempNAN:end,:)];
                    tempMatrixArray = [tempMatrixArray(1:tempNAN-1,:); nan(1,tempSize+1); tempMatrixArray(tempNAN:end,:)];
                end
            end

            % get the Column-wise sums of CEN, SN, DMN
            OutTable.AllConn_CEN(count) = sum(triu(tempAdjacencyMatrix(:,CEN_ROIs)), 'all');
            OutTable.AllConn_DMN(count) = sum(triu(tempAdjacencyMatrix(:,DMN_ROIs)), 'all');
            OutTable.AllConn_SN(count) = sum(triu(tempAdjacencyMatrix(:,SN_ROIs)), 'all');

            % get the intra-network sums of CEN, SN, DMN
            OutTable.Intra_CEN(count) = sum(triu(tempAdjacencyMatrix(CEN_ROIs,CEN_ROIs)), 'all');
            OutTable.Intra_DMN(count) = sum(triu(tempAdjacencyMatrix(DMN_ROIs,DMN_ROIs)), 'all');
            OutTable.Intra_SN(count) = sum(triu(tempAdjacencyMatrix(SN_ROIs,SN_ROIs)), 'all');

            % get the extra-network sums of CEN, SN, DMN
            OutTable.Extra_CEN(count) = sum(triu(tempAdjacencyMatrix(NotCEN_ROIs,CEN_ROIs)), 'all');
            OutTable.Extra_DMN(count) = sum(triu(tempAdjacencyMatrix(NotDMN_ROIs,DMN_ROIs)), 'all');
            OutTable.Extra_SN(count) = sum(triu(tempAdjacencyMatrix(NotSN_ROIs,SN_ROIs)), 'all');

            % get the inter-network sums of CEN, SN, DMN
            OutTable.CEN_SN(count) = sum(triu(tempAdjacencyMatrix(CEN_ROIs,SN_ROIs)), 'all');
            OutTable.SN_DMN(count) = sum(triu(tempAdjacencyMatrix(SN_ROIs,DMN_ROIs)), 'all');
            OutTable.DMN_CEN(count) = sum(triu(tempAdjacencyMatrix(DMN_ROIs,CEN_ROIs)), 'all');

            % get the intra-network STRENGTHS of CEN, SN, DMN
            OutTable.Intra_CEN_z(count) = mean(nonzeros(triu(tempMatrixArray(CEN_ROIs, CEN_ROIs),1)), 'all', 'omitnan');
            OutTable.Intra_DMN_z(count) = mean(nonzeros(triu(tempMatrixArray(DMN_ROIs, DMN_ROIs),1)), 'all', 'omitnan');
            OutTable.Intra_SN_z(count) = mean(nonzeros(triu(tempMatrixArray(SN_ROIs, SN_ROIs),1)), 'all', 'omitnan');

            % get the extra-network STRENGTHS of CEN, SN, DMN
            OutTable.Extra_CEN_z(count) = mean(nonzeros(triu(tempMatrixArray(NotCEN_ROIs, CEN_ROIs),1)), 'all', 'omitnan');
            OutTable.Extra_DMN_z(count) = mean(nonzeros(triu(tempMatrixArray(NotDMN_ROIs, DMN_ROIs),1)), 'all', 'omitnan');
            OutTable.Extra_SN_z(count) = mean(nonzeros(triu(tempMatrixArray(NotSN_ROIs, SN_ROIs),1)), 'all', 'omitnan');

            % get the extra-network STRENGTHS of CEN, SN, DMN
            OutTable.CEN_SN_z(count) = mean(nonzeros(triu(tempMatrixArray(CEN_ROIs, SN_ROIs),1)), 'all', 'omitnan');
            OutTable.SN_DMN_z(count) = mean(nonzeros(triu(tempMatrixArray(SN_ROIs, DMN_ROIs),1)), 'all', 'omitnan');
            OutTable.DMN_CEN_z(count) = mean(nonzeros(triu(tempMatrixArray(DMN_ROIs, CEN_ROIs),1)), 'all', 'omitnan');

            % populate any other detail
            OutTable.ID(count) = tempSubj;
            OutTable.Task(count) = tempTask;
            OutTable.Extraversion(count) = tempExtraversionRaw;
            OutTable.Extravert(count) = tempExtrovert;

            % add to the counter; to move to the next row
            count = count + 1;
        end
    end
    disp(['Completed: ' tempSubj ' of ' num2str(NumSubj)])
end


%% save the output file
% only get the non-empty rows
KeepRows = zeros(4000,1);
for iRow = 1:4000
    if OutTable.ID(iRow) ~= ""
        KeepRows(iRow) = 1;
    end
end
KeepRows = logical(KeepRows);
ID = OutTable.ID(KeepRows)';
Task = OutTable.Task(KeepRows)';
Extraversion = OutTable.Extraversion(KeepRows)';
Extravert = OutTable.Extravert(KeepRows)';
AllConn_CEN = OutTable.AllConn_CEN(KeepRows)';
AllConn_DMN = OutTable.AllConn_DMN(KeepRows)';
AllConn_SN = OutTable.AllConn_SN(KeepRows)';
Intra_CEN = OutTable.Intra_CEN(KeepRows)';
Intra_DMN = OutTable.Intra_DMN(KeepRows)';
Intra_SN = OutTable.Intra_SN(KeepRows)';
Extra_CEN = OutTable.Extra_CEN(KeepRows)';
Extra_DMN = OutTable.Extra_DMN(KeepRows)';
Extra_SN = OutTable.Extra_SN(KeepRows)';
CEN_SN = OutTable.CEN_SN(KeepRows)';
SN_DMN = OutTable.SN_DMN(KeepRows)';
DMN_CEN = OutTable.DMN_CEN(KeepRows)';
Intra_CEN_z = OutTable.Intra_CEN_z(KeepRows)';
Intra_DMN_z = OutTable.Intra_DMN_z(KeepRows)';
Intra_SN_z = OutTable.Intra_SN_z(KeepRows)';
Extra_CEN_z = OutTable.Extra_CEN_z(KeepRows)';
Extra_DMN_z = OutTable.Extra_DMN_z(KeepRows)';
Extra_SN_z = OutTable.Extra_SN_z(KeepRows)';
CEN_SN_z = OutTable.CEN_SN_z(KeepRows)';
SN_DMN_z = OutTable.SN_DMN_z(KeepRows)';
DMN_CEN_z = OutTable.DMN_CEN_z(KeepRows)';

OutTable_table = table(ID, Task, Extraversion, Extravert, ...
    AllConn_CEN, AllConn_DMN, AllConn_SN, ...
    Intra_CEN, Intra_DMN, Intra_SN, ...
    Extra_CEN, Extra_DMN, Extra_SN, ...
    CEN_SN, SN_DMN, DMN_CEN, ...
    Intra_CEN_z, Intra_DMN_z, Intra_SN_z, ...
    Extra_CEN_z, Extra_DMN_z, Extra_SN_z, ...
    CEN_SN_z, SN_DMN_z, DMN_CEN_z);
OutTable_name = [outputdir '/MnM_fig34_dFC_EdgeAnalyses.csv'];
writetable(OutTable_table, OutTable_name, 'Delimiter', ',', 'FileType', 'text');

disp('Script is complete.')