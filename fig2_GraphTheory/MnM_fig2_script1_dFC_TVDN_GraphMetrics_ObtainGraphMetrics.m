%% Set-up
close all; clear; clc

inputdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/MnM_dFC/r_0.5/k_1.6';
ratingsdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/events';
scriptdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub';

addpath('/Users/josephchen/Library/CloudStorage/GoogleDrive-nzjosephchen@gmail.com/My Drive/nzjosephchen/software/2019_03_03_BCT')
addpath('/Users/josephchen/Library/CloudStorage/GoogleDrive-nzjosephchen@gmail.com/My Drive/nzjosephchen/software/NetworkVisualizer')

PersonalityData = readtable('/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/MnM_PersonalityInfo.csv');

%% Analyze Graph Metrics

cd(inputdir);
df1 = readtable('df1_modified.csv', 'FileType', 'text', 'delimiter', ',');

% initialise empty matrices for EVERY FC matrix network metric
OutputDF_all = [];
OutputDF_all.ID = strings(1,4000);
OutputDF_all.Task  = strings(1,4000);
OutputDF_all.Extraversion = strings(1,4000);
OutputDF_all.Extravert  = strings(1,4000);
OutputDF_all.GEFF  = strings(1,4000);
OutputDF_all.Modularity = strings(1,4000);

% initialise empty matrices for EVERY FC matrix network metric
OutputDF_summ = [];
OutputDF_summ.ID = strings(1,200);
OutputDF_summ.Task  = strings(1,200);
OutputDF_summ.Extraversion = strings(1,200);
OutputDF_summ.Extravert  = strings(1,200);
OutputDF_summ.Attention = strings(1,200);
OutputDF_summ.GEFF  = strings(1,200);
OutputDF_summ.Modularity = strings(1,200);

count = 1;
summarycount = 1;
AllSubj = dir('sub-*');
AllSubj = {AllSubj(:).name};
NumSubj = length(AllSubj);
for iSubj = 1:length(AllSubj)
    tempSubj = AllSubj{iSubj};

    % get extroversion scores
    [~,loc] = ismember(tempSubj, PersonalityData.ID);
    tempExtraversionRaw = PersonalityData.Personality_ExtraversionRaw(loc);
    if tempExtraversionRaw < 24; tempExtravert = 0; end
    if tempExtraversionRaw > 24; tempExtravert = 1; end
    if tempExtraversionRaw == 24; tempExtravert = nan; end


    AllTasks = {'BreathNVHA' 'BreathNVLA' 'Meditation'};
    for iTask = 1:3
        tempTask = AllTasks{iTask};
        tempFCmatrices = dir([tempSubj '/' tempTask '/FC_*']);
        tempNumFCmatrices = length(tempFCmatrices);

        % get average attention score within this task
        tempbehavdffile = [ratingsdir '/' tempSubj '_task-' tempTask '_events.csv'];
        tempbehavdf = readtable(tempbehavdffile);
        tempresponsetrials = find(contains(tempbehavdf.trial_type, "response"));
        tempAvgAttn = mean(tempbehavdf.rating(tempresponsetrials));

        % create empty summary matrix
        tempAllGEFF = zeros(tempNumFCmatrices,1);
        tempAllModularity = zeros(tempNumFCmatrices,1);

        for iFC = 1:tempNumFCmatrices
            tempFCnum = sprintf('%02d', iFC);
            tempFCfile = [tempSubj '/' tempTask '/FC_matrix-' tempFCnum '.csv'];
            tempFC = readtable(tempFCfile, 'FileType', 'text', ...
                'delimiter', ',', 'ReadVariableNames', false);
            tempFC = table2array(tempFC);

            % now apply threshold to connectivity matrix
            tempSortValues = sort(abs(tempFC(:)), 'descend');
            tempSortValues = tempSortValues(~isnan(tempSortValues));
            tempThreshold = tempSortValues(floor(0.2*length(tempSortValues)));
            tempAdjacencyMatrix = abs(tempFC) >= tempThreshold;

            % network analysis
                     
            % ---- Global measures of integration -----
            % Global efficiency - The average of the inverse shortest path between two points of the network
            geff = efficiency_bin(tempAdjacencyMatrix);
            tempAllGEFF(iFC) = geff;

            % -----Modularity -----
            [~,modularity] = modularity_und(tempAdjacencyMatrix);
            tempAllModularity(iFC) = modularity;

            OutputDF_all.ID(count) = tempSubj;
            OutputDF_all.Task(count)  = tempTask;
            OutputDF_all.Extraversion(count) = num2str(tempExtraversionRaw);
            OutputDF_all.Extravert(count)  = num2str(tempExtravert);
            OutputDF_all.GEFF(count) = num2str(geff);
            OutputDF_all.Modularity(count) = num2str(modularity);
            count = count+1;
        end

        OutputDF_summ.ID(summarycount) = tempSubj;
        OutputDF_summ.Task(summarycount)  = tempTask;
        OutputDF_summ.Extraversion(summarycount) = tempExtraversionRaw;
        OutputDF_summ.Extravert(summarycount)  = tempExtravert;
        OutputDF_summ.Attention(summarycount) = tempAvgAttn;
        OutputDF_summ.GEFF(summarycount)  = mean(tempAllGEFF);
        OutputDF_summ.Modularity(summarycount) = mean(tempAllModularity);
        summarycount = summarycount+1;

        disp(['Completed ' num2str(iSubj) ' of ' num2str(NumSubj) '. Task: ' tempTask]);
    end
end

%% Save outputs
% save outputs in final analysis directories
OutputDFFileALL = 'graphmatrics_ALL.csv';
OutputDFFileSummary = 'graphmatrics_summary.csv';

% save the output_all file
% only get the non-empty rows
KeepRows = zeros(4000,1);
for iRow = 1:4000
    if OutputDF_all.ID(iRow) ~= ""
        KeepRows(iRow) = 1;
    end
end
KeepRows = logical(KeepRows);
ID = OutputDF_all.ID(KeepRows)';
Task = OutputDF_all.Task(KeepRows)';
Extraversion = OutputDF_all.Extraversion(KeepRows)';
Extravert = OutputDF_all.Extravert(KeepRows)';
GEFF = OutputDF_all.GEFF(KeepRows)';
Modularity = OutputDF_all.Modularity(KeepRows)';
OutputDF_all_table = table(ID, Task, Extraversion, Extravert,GEFF, Modularity);
writetable(OutputDF_all_table, OutputDFFileALL, 'Delimiter', ',', 'FileType', 'text');

% save the output_summ file
% only get the non-empty rows
KeepRows = zeros(200,1);
for iRow = 1:200
    if OutputDF_summ.ID(iRow) ~= ""
        KeepRows(iRow) = 1;
    end
end
KeepRows = logical(KeepRows);
ID = OutputDF_summ.ID(KeepRows)';
Task = OutputDF_summ.Task(KeepRows)';
Extraversion = OutputDF_summ.Extraversion(KeepRows)';
Extravert = OutputDF_summ.Extravert(KeepRows)';
Attention = OutputDF_summ.Attention(KeepRows)';
GEFF = OutputDF_summ.GEFF(KeepRows)';
Modularity = OutputDF_summ.Modularity(KeepRows)';
OutputDF_summ_table = table(ID, Task, Extraversion, Extravert, GEFF, Modularity);
writetable(OutputDF_summ_table, OutputDFFileSummary, 'Delimiter', ',', 'FileType', 'text');

