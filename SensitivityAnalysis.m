%Created by Ignacio X. Partarrieu
%% Problem formulation
% We wish to carry out a sensitivity analysis of a large matrix using eta
% squared. This code imports the relevant excel spreadsheets and
% concatenates them for analysis. It then does a sensitivity analysis of
% all effects of interest.

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

%We want to see the effect of the all factors 
factorNames = {'Model', 'GL Bins', 'Wavelength','Reconstruction'};

%Manual ANOVA implementation

for kk = 1:numel(radInd)
    [intFactorSSE,intFactorNames] = ...
        manualAnova2(Tnew{:,radInd(kk)},...
        {Tnew{:,3},Tnew{:,2},Tnew{:,4},...
        Tnew{:,5}},...
        factorNames);
    radSS = intFactorSSE(1:4)/sum(intFactorSSE);%extracts the feature SS values from anova table and converts them to fractions
    radSS(5) = sum(intFactorSSE(5:end))/sum(intFactorSSE);
    if kk == 1
        SS = array2table(radSS,'VariableNames',...
            {Tnew.Properties.VariableNames{radInd(kk)}});
    else
        SS = addvars(SS,radSS,'NewVariableNames',...
            Tnew.Properties.VariableNames{radInd(kk)});
    end
    SSTable = SS;
end

%% Now the k fold analysis
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

        for kk = 1:numel(radInd)
            [intFactorSSE,intFactorNames] = ...
                manualAnova2(Tfold{:,radInd(kk)},...
                {Tfold{:,3},Tfold{:,2},Tfold{:,4},Tfold{:,5}},...
                factorNames);
            radSS = intFactorSSE(1:4)/sum(intFactorSSE);%extracts the feature SS values from anova table and converts them to fractions
            radSS(5) = sum(intFactorSSE(5:end))/sum(intFactorSSE);
            if kk == 1
                SS = array2table(radSS,'VariableNames',...
                    {Tfold.Properties.VariableNames{radInd(kk)}});
            else
                SS = addvars(SS,radSS,'NewVariableNames',...
                    Tfold.Properties.VariableNames{radInd(kk)});
            end
            SSTableFold{ff} = SS;
        end
        writetable(SSTableFold{ff},['MainFold',num2str(ff),'.xlsx'],'Sheet',...
            'AllValue');  
end

%% Now plot the graphs
%We want to seperate the different feature groups using whitespaces so we
%do a bit of pre-processing
close all

    %add whitespace to SSTable plot
    SSTable = addvars(SSTable,[0;0;0;0;0],'After',...
        'original_firstorder_Variance');
    SSTable = addvars(SSTable,[0;0;0;0;0],'After',...
        'original_glcm_SumSquares');
    SSTable = addvars(SSTable,[0;0;0;0;0],'After',...
        'original_gldm_SmallDependenceLowGrayLevelEmphasis');
    SSTable = addvars(SSTable,[0;0;0;0;0],'After',...
        'original_glrlm_ShortRunLowGrayLevelEmphasis');
    SSTable = addvars(SSTable,[0;0;0;0;0],'After',...
        'original_glszm_ZoneVariance');
    figure('WindowState','maximized'),
    ba = bar(table2array(SSTable)','stacked', 'FaceColor','flat')
    ba(1).CData = [34 136 51]/255;
    ba(2).CData = [102 204 238]/255;
    ba(3).CData = [204 187 68]/255;
    ba(4).CData = [170 51 119]/255;
    ba(5).CData = [0 0 0]/255;
    legend({'Model','GL Bins','Wavelength','Reconstruction','Error'},'Location','NorthEastOutside')
    ylim([0 1.05]);
    xticks([10 21 31 41 52 63 73 84 95])
    xticklabels({'10','20','30','40','50','60','70','80','90'})
    labels = {'FOS','GLCM','GLDM','GLRLM','GLSZM','NGTDM'};
    X = [1,20,45,60,77,94];
    Y = [1.07,1.07,1.07,1.07,1.07,1.07];
    text(X,Y,labels,'HorizontalAlignment','left','VerticalAlignment',...
        'bottom','FontSize',20)
    set(gca,'FontSize',20)
    saveas(gcf,['GL',num2str(binLevels(jjj))],'png')
