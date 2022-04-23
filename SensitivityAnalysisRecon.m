%Created by Ignacio X. Partarrieu
%% Problem formulation
% We wish to carry out a sensitivity analysis of a large matrix using eta
% squared. This code imports the relevant excel spreadsheets and
% concatenates them for analysis. This script specifically focuses on
% analysing reconstruction

%% Approach
% We load in the data to use, look at the columns representing factors and
% identify the unique levels. Once established, we loop through the table
% to calculate the group sum of squares using the following definition: 
% https://mathworld.wolfram.com/ANOVA.html and https://person.hst.aau.dk/cdahl/BiostatPhD/ANOVA. 
%% Practicalities
% In order to do the ANOVA calculation, we balance the classes by removing
% a single luminal case for the full analysis to have 10 cases of each. We
% also perform the ANOVA on the standardised cases.
clearvars

%find csv files to import
dirpath = "C:\Users\Coding\Desktop\SensitivityAnalysis\CambridgeColabPaperVersions2022\Radiomics";
searchTerm = "*21*.csv";
dirNames = dir(fullfile(dirpath,searchTerm));

%import and concatenate files into a single table for the analysis
T0 = [];
for iii = 1:numel(dirNames)
    %import single table at a time
    T = readtable(fullfile(dirpath,dirNames(iii).name));
    if contains(dirNames(iii).name,'MB')
        [ReconType{1:size(T,1)}] = deal('MB');%Model linear
    else
        [ReconType{1:size(T,1)}] = deal('BP');%Backprojection
    end
    T = addvars(T,ReconType','Before','SliceThickness');
    T0 = [T0;T];
end
%Rename added column
T0 = renamevars(T0,'Var5','ReconType');

%% Now we run the sensitivity analysis

%First, remove randomly selected luminal patient
Tnew = T0;
toDelete = Tnew.PatientName == 246;
Tnew(toDelete,:) = [];

%And declare indices of radiomics features for analysis
radInd = 43:135;

%We want to see the effect of the reconstruction so find indices of each bin for analysis
reconLevels = unique(Tnew.ReconType);
factorNames = {'Model', 'Wavelength','GL Bins'};

%Manual ANOVA implementation since MATLAB changed the code after 2018 to
%include undesired smoothing
for jjj = 1:numel(reconLevels)
    toAnalyze = ismember(Tnew.ReconType, reconLevels{jjj});
    for kk = 1:numel(radInd)
        [intFactorSSE,intFactorNames] = ...
            manualAnova2(Tnew{toAnalyze,radInd(kk)},...
            {Tnew{toAnalyze,3},Tnew{toAnalyze,4},Tnew{toAnalyze,2}},...
            factorNames);
        radSS = intFactorSSE(1:3)/sum(intFactorSSE);%extracts the feature SS values from anova table and converts them to fractions
        radSS(4) = sum(intFactorSSE(4:end))/sum(intFactorSSE);
        if kk == 1
            SS = array2table(radSS,'VariableNames',...
                {Tnew.Properties.VariableNames{radInd(kk)}});
        else
            SS = addvars(SS,radSS,'NewVariableNames',...
                Tnew.Properties.VariableNames{radInd(kk)});
        end
        SSTable{jjj} = SS;
    end
end

%% do the k fold analysis and save results to excel

% find which PatientNames are associated with each model
[~,~,IndModels]=unique(T0(:,3));
basalNum = table2array(unique(T0(IndModels==1,1)));
luminalNum = table2array(unique(T0(IndModels==2,1)));

%Do fivefold analysis, randomly selecting 3 luminal and 2 basal cases each
%time
for ff = 1:5
    Tfold = T0;
    toDeleteFoldVal = [basalNum(randperm(numel(basalNum),2),1);...
        luminalNum(randperm(numel(luminalNum),3),1)];
    toDeleteFold = ismember(Tfold.PatientName,toDeleteFoldVal(:));
    Tfold(toDeleteFold,:) = [];
    %Manual ANOVA implementationsince MATLAB changed the code after 2018 to
    %include undesired smoothing
    for jjj = 1:numel(reconLevels)
        toAnalyze = ismember(Tfold.ReconType, reconLevels{jjj});
        for kk = 1:numel(radInd)
            [intFactorSSE,intFactorNames] = ...
                manualAnova2(Tfold{toAnalyze,radInd(kk)},...
                {Tfold{toAnalyze,3},Tfold{toAnalyze,4},Tfold{toAnalyze,2}},...
                factorNames);
            radSS = intFactorSSE(1:3)/sum(intFactorSSE);%extracts the feature SS values from anova table and converts them to fractions
            radSS(4) = sum(intFactorSSE(4:end))/sum(intFactorSSE);
            if kk == 1
                SS = array2table(radSS,'VariableNames',...
                    {Tfold.Properties.VariableNames{radInd(kk)}});
            else
                SS = addvars(SS,radSS,'NewVariableNames',...
                    Tfold.Properties.VariableNames{radInd(kk)});
            end
            SSTableFold{jjj,ff} = SS;
        end
        writetable(SSTableFold{jjj,ff},['ReconFold',num2str(ff),'.xlsx'],'Sheet',...
            ['ReconValueEq',reconLevels{jjj}]);
    end    
end

%% Now plot the graphs
%We want to seperate the different feature groups using whitespaces so we
%do a bit of pre-processing
close all
for jjj = 1:numel(reconLevels)
    %add whitespace to SSTable plot
    SSTable{jjj} = addvars(SSTable{jjj},[0;0;0;0],'After',...
        'original_firstorder_Variance');
    SSTable{jjj} = addvars(SSTable{jjj},[0;0;0;0],'After',...
        'original_glcm_SumSquares');
    SSTable{jjj} = addvars(SSTable{jjj},[0;0;0;0],'After',...
        'original_gldm_SmallDependenceLowGrayLevelEmphasis');
    SSTable{jjj} = addvars(SSTable{jjj},[0;0;0;0],'After',...
        'original_glrlm_ShortRunLowGrayLevelEmphasis');
    SSTable{jjj} = addvars(SSTable{jjj},[0;0;0;0],'After',...
        'original_glszm_ZoneVariance');
    figure('WindowState','maximized'),
    ba = bar(table2array(SSTable{jjj})','stacked', 'FaceColor','flat')
    ba(1).CData = [34 136 51]/255;
    ba(2).CData = [204 187 68]/255;
    ba(3).CData = [102 204 238]/255;
    ba(4).CData = [0 0 0]/255;
    legend({'Model','Wavelength','GL Bins','Error'},'Location','NorthEastOutside')
    ylim([0 1.05]);
    xticks([10 21 31 41 52 63 73 84 95])
    xticklabels({'10','20','30','40','50','60','70','80','90'})
    labels = {'FOS','GLCM','GLDM','GLRLM','GLSZM','NGTDM'};
    X = [1,20,45,60,77,94];
    Y = [1.07,1.07,1.07,1.07,1.07,1.07];
    text(X,Y,labels,'HorizontalAlignment','left','VerticalAlignment',...
        'bottom','FontSize',20)
    set(gca,'FontSize',20)
    saveas(gcf,reconLevels{jjj},'png')
end