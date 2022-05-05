function [intFactorSSE,intFactorNames] = manualAnova2(data,factors,factorNames)
%% Name of project and work package go here
% Code written by Ignacio X. Partarrieu 18/12/2020
% ver 1.0 - confidently returns the same values as the MATLAB 2018b and
% before versions
%% Problem formulation
% We wish to carry out ANOVA on a large fully sampled factorial matrix. As
% MATLAB anova function chokes on lack of memory, we implement
% this calculation as a for loop, which will produce a result

%% Approach
% We load in the data to use, look at the columns representing factors and
% identify the unique levels. Once established, we loop through the table
% to calculate the group sum of squares using the following definition: 
% https://mathworld.wolfram.com/ANOVA.html and https://person.hst.aau.dk/cdahl/BiostatPhD/ANOVA. 
%% Practicalities
% In order to do the ANOVA calculation, we use the nchoosek function to
% calculate all possible interations between the terms. We also calculate
% the grand mean. These sets of information are used to calculate the sum
% of squares of the groupings. Interactions are then calculated by looking
% at the within cell sum of squares and substracting the preceeding factors
% and interaction terms. In order to determine which matrix cells contain
% relevant information we index the matrices as we build, making subsequent
% calculations easier
%data: is the array of numerical data to use
%factors: is the list of labels for the data
%factorNames: name of each label
%intFactorSSE: the output sum of squares
%intFactorNames: the associated factor 

%% version details
% -----------------------------------------------------------------------------------------------------
% MATLAB Version: 9.5.0.944444 (R2018b)
% MATLAB License Number: 40761551
% Operating System: Microsoft Windows 10 Enterprise Version 10.0 (Build 17134)
% Java Version: Java 1.8.0_152-b16 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% -----------------------------------------------------------------------------------------------------
% MATLAB                                                Version 9.5         (R2018b)
% Communications Toolbox                                Version 7.0         (R2018b)
% Curve Fitting Toolbox                                 Version 3.5.8       (R2018b)
% DSP System Toolbox                                    Version 9.7         (R2018b)
% Image Processing Toolbox                              Version 10.3        (R2018b)
% Parallel Computing Toolbox                            Version 6.13        (R2018b)
% Signal Processing Toolbox                             Version 8.1         (R2018b)
% Statistics and Machine Learning Toolbox               Version 11.4        (R2018b)
% Symbolic Math Toolbox 
numMeas = numel(data);
numFactors = numel(factors);

%calculate total sum of squares
grandMeandata=mean(data);

%Find and label all unique occurences
ic = zeros(numMeas,numFactors);
for iii = 1:numFactors
    [~,~,ic(:,iii)]=unique(factors{iii});
end

for kk = 1:numFactors%number of interactions
    %list all nth order interactions
    possibleComb = nchoosek(1:numFactors,kk);
    for fff = 1:size(possibleComb,1)
        int = ic(:,possibleComb(fff,:));%goes through all possible interations

        %As we are now dealing with multiple columns, 
        %we write code to condense the indexing into a single one
        %this will allow for an easier search and match
        maxPowerOfTen = size(int,2)-1;
        A1 = sum(int .* 10.^(maxPowerOfTen:-1:0),2);
        %We use the indexes to calculate the means of the data chosen
        [cint,~,~]=unique(A1);
        multiplier = numMeas/numel(cint);
        for jjj = 1:numel(cint)
            intlevelmean{kk,fff}(jjj) = mean(data(A1==...
                cint(jjj)));
        end
        
        %Calculate cell sum of squares
        cellFactorSSE(kk,fff)=multiplier*sum((intlevelmean{kk,fff}-grandMeandata).^2);
        %calculate the sum of squares of factors and their interactions
        %(type 3 SS)
        if kk == 1
            intFactorSSE(kk,fff) = cellFactorSSE(kk,fff);
        else
            SStoRemove{kk,fff} = 0;
            % this loop does the 'magic' of finding the values which are to
            % be subtracted from the cell sum of squares, and adding them
            % up
            for ggg = 1:size(possibleComb,2)-1  
                combSearch = nchoosek(possibleComb(fff,:),ggg);
                logicIndex = ismember(indFactor{ggg},combSearch,'rows');
                SStoRemove{kk,fff} = SStoRemove{kk,fff}+sum(intFactorSSE(ggg,...
                    logicIndex));
            end
            %This is our result!! Importantly, it stores the sum of squares
            %of the various factors and interactions. Interactions are on
            %each row (i.e. row 1: no interactions, row 2, 2 level
            %interactions (i.e. AxB, AxC, AxD,... BxC, BxD,... etc) and so
            %on
            intFactorSSE(kk,fff) = cellFactorSSE(kk,fff) - SStoRemove{kk,fff};
        end
        %this index becomes usefull after the first looping, and gets used
        %by the loop ggg as an index to compare against stating which
        %factors were used
        indFactor{kk}(fff,:) = possibleComb(fff,:);
        clear int
        clear maxPowerOfTen
        clear A1
        clear cint
        clear multiplier

    end
end

%We also create a labelling matrix for the interactions, where
%intFactorNames cells correspond to intFactorSSE cells 
for kk = 1:7 
    possibleComb = nchoosek(1:numFactors,kk);
    for fff = 1:size(possibleComb,1)
        factNames = factorNames{possibleComb(fff,1)};
        for ggg = 2:size(possibleComb,2)
            factNames = join([factNames,factorNames{possibleComb(fff,ggg)}],'x');
        end 
        intFactorNames{kk,fff} = factNames;
        clear factNames;
    end
end
    
intFactorSSE = intFactorSSE';
intFactorSSE = intFactorSSE(:);
intFactorNames = intFactorNames';
intFactorNames = intFactorNames(:);
end

