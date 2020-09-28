clear variables; close all; clc;

%% =========== Internetwork connectivity (inter-FC) analysis =========== %%

% Performs INC among FC networks
% Select IC of interest. Vector 1xn with n as the index of the ICs of
% interest 1:20; for all networks.
% Written by: M.E. Archila-Meléndez and modified by S. Küchenhoff and
% A.L. Ruiz-Rizzo

%% General settings =======================================================

saveVars = 0;          % 0 to test and 1 to save .mat .xls and .pngs   
subjectToExclude = 62; % It has to be 1 more than the actual file name/
                       % number (e.g. Subj003 subjectToExclude = 4);
                       % "62" is a mock
timePoints = 595;      % Number of fMRI time points (rows in files)
dirData  = ['/Users/lmuresearchfellowship/Documents/Adriana'...
    '/LMU_Psychology/Projects/Svenja/Global_Young/']; % where the data...
                       % (dr_stage_1) are stored.
                       % put the full path if this script is not in
                       % the same one where the data are stored
Frequency = 'Global';

%% Settings application ===================================================

%   Networks of interest
    intNets  = [44, 3, 4, 11, 1, 5, 14]; % Vector

    NetNames = {'COn', 'RFPn', 'Vis-39', 'Vis-46',...
        'Vis-59', 'Vis-64', 'Vis-67'}; % Names
                   
%% Load the data set (time courses) of subjects ===========================

cd Global_Young
D = dir(fullfile('*.txt'));
ii = 1;

for i=(1:length(D))
if i==subjectToExclude
 disp (['**EXCLUDED**:dr_stage1_subject000', num2str(i-1),'.txt'])
else
       if(i<11)
            subjCON{1,ii} = load([dirData, filesep,...
                'dr_stage1_subject0000',...
                num2str(i-1),'.txt']);
             disp (['Included: ' 'dr_stage1_subject0000',...
                 num2str(i-1),'.txt'])
       else
           (i>=11 && i<=99);
           subjCON{1,ii} = load([dirData, filesep,...
                'dr_stage1_subject000',...
                num2str(i-1),'.txt']);
            disp (['Included: ' 'dr_stage1_subject0000',...
                num2str(i-1),'.txt'])
       end
        
        ii = ii + 1;
end
end

clear ii
clear i % subjCon is now the complete dataset

%% Extraction of ICs of interest ==========================================

SelectSubData = zeros(timePoints,length(NetNames));
% rows: number of time points or volumes;
% colums: number of ICs selected
for i=1:numel(subjCON)
    currSubjData = subjCON{1,i};
    SelectSubData = currSubjData(:,intNets(1,:));
    subjCON{1,i} = SelectSubData;
end

clear i currSubjData 

%% Compute the correlation between ICs for subjCON subject by subject =====

for i = 1:numel(subjCON)        % subjects
    RCON{i} = corr(subjCON{i}); % correlations within each subject's file
end

% r-to-Z transform with Fisher for all ICs for each data file
for i = 1:numel(subjCON)
    for j = 1:numel(intNets)
        for k = 1:numel(intNets)
            fisherZCON{i}(j,k) = 0.5 * (log((1 + RCON{i}(j,k)) /...
                (1 - RCON{i}(j,k)))); % Fisher Z transformation formula
        end
    end
end

% The next step will concatenate all z-transformation correlation... 
% ...matrices per subject in one file. But we can use the index to... 
% ...create a vector with the Subjects Numbers to later name the per...
% ...subject result table

%% Concatenate all z-transformed correlation matrices =====================

SubjNum = [];
FullFisherZCON = [];
for i = 1:size(fisherZCON,2) % returns the number of columns
    FullFisherZCON = cat (3, FullFisherZCON, fisherZCON{1,i});
    SubjNum = [SubjNum; i];
end

%% Extract Z values =======================================================

% E.g., if you want to use them outside of Matlab,
% select the specific network numbers that you want the Z-value for

% Write manually the heading you want (content of cells)
Names_ExtractZval = {'RFP and CO','RFP and 39', 'RFP and 46',...
    'RFP and 59', 'RFP and 64', 'RFP and 67'}; 
ExtractZval = [];

for i = 1:size(fisherZCON,2)
    tempZval(1,1) = FullFisherZCON(2,1,i);
    tempZval(1,2) = FullFisherZCON(2,3,i);
    tempZval(1,3) = FullFisherZCON(2,4,i);
    tempZval(1,4) = FullFisherZCON(2,5,i);
    tempZval(1,5) = FullFisherZCON(2,6,i);
    tempZval(1,6) = FullFisherZCON(2,7,i);
    tempZval = num2cell(tempZval);
    Names_ExtractZval = [Names_ExtractZval; tempZval];
    ExtractZval = cat (1, ExtractZval, tempZval);
    clear tempZval
end

% Write manually the heading you want (content of cells)
Names_ExtractZval2 = {'CO and 39', 'CO and 46',...
    'CO and 59', 'CO and 64', 'CO and 67'}; 
ExtractZval2 = [];

for i = 1:size(fisherZCON,2)
    tempZval(1,1) = FullFisherZCON(1,3,i);
    tempZval(1,2) = FullFisherZCON(1,4,i);
    tempZval(1,3) = FullFisherZCON(1,5,i);
    tempZval(1,4) = FullFisherZCON(1,6,i);
    tempZval(1,5) = FullFisherZCON(1,7,i);
    tempZval = num2cell(tempZval);
    Names_ExtractZval2 = [Names_ExtractZval2; tempZval];
    ExtractZval2 = cat (1, ExtractZval2, tempZval);
    clear tempZval
end

Names_ExtractZval = [Names_ExtractZval, Names_ExtractZval2];
ExtractZval = [ExtractZval, ExtractZval2];

% Add the subject numbers to the matrix
xm = 0;
SubjNum = [xm; SubjNum];
SubjNum = num2cell (SubjNum);
Names_ExtractZval = [SubjNum, Names_ExtractZval];

clear i k j

%% Significance section ===================================================
    valid        = zeros(numel(intNets));
    
% One-sample t-test for NON-independent samples
    for x = 1:length(FullFisherZCON(:,1,1))
        for y = 1:length(FullFisherZCON(1,:,1))
            [h(x,y),pval(x,y),ci,stats] = ttest(FullFisherZCON(x,y,:));
            
            ciAll(x,y,1)        = ci(:,:,1);
            ciAll(x,y,2)        = ci(:,:,2);
            
            statsAll.tstat(x,y) = stats.tstat;
            statsAll.df(x,y)    = stats.df;
            statsAll.sd(x,y)    = stats.sd;

            if (x < y)
                valid(y,x) = valid(y,x) + 1; % create lower diagonal matrix
            end
        end
    end
    
    Zavg = mean(FullFisherZCON,3);
    
    clear x y 
    
   % Significance output
   ValidVec    = find(valid==1); % only in the lower diagonal for all...
                                  % 3 pVals
    
   % Multiple comparison correction with FDR=false discovery rate for INC 
   q           = mafdr(pval(valid==1),'BHFDR','true');
   pval2       = pval(ValidVec);
   Rfdr        = zeros(numel(intNets));
   Rfdr(ValidVec(pval2 < 0.05 & q < 0.05)) = 1; % extracting the...
   % ...significant averaged INC
   [MarkI,MarkJ] = ind2sub(size(Rfdr),find(Rfdr==1)); % For drawing the...
   % ...stars (*) in each cell  

%% To display and save results ============================================

    % Figure for average group matrix
    figure('Colormap',redbluecmap)
    Zavg(Zavg == Inf) = 0;
    imagesc(tril(Zavg)); colorbar;...
    % title([ 'INC for ' Frequency...
    % ' (pval < 0.05 and q < 0.05)']);
    caxis ([-1, 1]);
    c = colorbar;
    c.Label.String = 'z values';
    text(MarkJ,MarkI,{'*'},'fontsize',12); % For drawing the stars (*)...
    % ...in each cell
    set(gca,'Xtick',(1:numel(intNets)))
    set(gca,'Ytick',(1:numel(intNets)))
    set(gca, 'xticklabel', NetNames, 'FontSize', 14);
    set(gca, 'yticklabel', NetNames, 'FontSize', 14);

    if saveVars == 1
    saveas(gcf,[Frequency '_INC_Zvals_' num2str(length (intNets))...
        '_Nets_for_' date '.png']);
    end
   
    %% To have have headings of networks for the result tables ============
    
    % Create a vertical vector of the networks
    vertintNets = intNets';
    Networksvert = [xm; vertintNets];
    
    % Put the networks above and on the left of the result matrices
    Names_Zavg = [intNets; Zavg];
    Names_Zavg = [Networksvert, Names_Zavg];
    Names_pval = [intNets; pval];
    Names_pval = [Networksvert, Names_pval];
    Names_statsAll.df = [intNets; statsAll.df];
    Names_statsAll.df = [Networksvert, Names_statsAll.df];
    Names_statsAll.sd = [intNets; statsAll.sd];
    Names_statsAll.sd = [Networksvert, Names_statsAll.sd];
     
    %% To extract Result values ===========================================

    results{1}     = Names_Zavg;
    results{2}     = Names_pval;
    pvalTemp       = zeros(size(pval));
    pvalTemp(ValidVec(pval2 < 0.05 & q < 0.05)) = ...
        -log10(pval(ValidVec(pval2 < 0.05 & q < 0.05)));
    results{3}     = pvalTemp;
    results{4}     = Names_statsAll;
    results{5}     = Names_ExtractZval;
   
    %% To save results values =============================================
    
    if saveVars == 1
    save (['resultsINC_' date]);
    xlswrite(([Frequency '_INC_for_' num2str(length (intNets))...
        '_Nets_' 'MeanZval_' date '.xlsx']),results{1})
    xlswrite(([Frequency '_INC_for_' num2str(length (intNets))...
        '_Nets_' 'pValUncorr_' date '.xlsx']),results{2})
    xlswrite(([Frequency '_INC_for_' num2str(length (intNets))...
        '_Nets_' 'Log10_pVal_Correct_' date '.xlsx']),results{3})
    xlswrite(([Frequency '_INC_for_' num2str(length (intNets))...
        '_Nets_' 'ttest_df_' date '.xlsx']),results{4}.df)
    xlswrite(([Frequency '_INC_for_' num2str(length (intNets))...
        '_Nets_' 'ttest_SD_' date '.xlsx']),results{4}.sd)
    % Convert cell to a table and use first row as variable names
    T = cell2table(results{5}(2:end,2:end),'VariableNames',...
        results{5}(1,2:end));
    % Write the table to a CSV file
    writetable(T, [Frequency '_INC_for_' num2str(length (intNets))...
        '_Nets_' 'Z_vals_' date '.csv'])
    end
   
saveVars;
disp '********** Finished **********'
