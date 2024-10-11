%% Script to compute HC and FEP MIND edges, degrees, and their effect.

% Copyright (C) 2024 University of Seville

% Written by Natalia García San Martín (ngarcia1@us.es)

% This file is part of Neurobiology MIND Psychosis toolkit.
%
% Neurobiology MIND Mapping toolkit is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
%
% Neurobiology MIND Mapping toolkit is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Neurobiology MIND Mapping toolkit. If not, see 
% <https://www.gnu.org/licenses/>.

clear
close all
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\';
read = 'Code\2. MIND\Data\MIND_networks_PAFIP\aparc\';
write_degree = 'Code\2. MIND\Data\degree\';
write_edges = 'Code\2. MIND\Data\edges\';
cd([location,write_degree])


if exist([location,write_degree,'degree_68_CN.csv'])
    degree_68_CN = readtable([location,write_degree,'degree_68_CN.csv'],"ReadRowNames",true);
    degree_68_FEP = readtable([location,write_degree,'degree_68_FEP.csv'],"ReadRowNames",true);
    load([location,write_edges,'edge_68_CN.mat'])
    load([location,write_edges,'edge_68_FEP.mat'])

else
    files = dir([location,read]);
    
    files = files(3:end); % all sessions
    files = files(contains({files.name},'_001')); % session 1
    
    i_FEP = 1;
    i_CN = 1;
    FEP_subjects = {};
    CN_subjects = {};
    for i = 1:length(files)
        MIND = readtable([[location,read],files(i).name],'ReadRowNames',true,'ReadVariableNames',true);
       
        % Sort alphabetically
        MIND = sortrows(MIND,'RowNames');
        col_order = sort(MIND.Properties.VariableNames);
        MIND = MIND(:,col_order);
        
        % Averaged hemispheres
        lh = MIND(contains(MIND.Properties.RowNames,'lh_'),contains(MIND.Properties.VariableNames,'lh_'));
        rh = MIND(contains(MIND.Properties.RowNames,'rh_'),contains(MIND.Properties.VariableNames,'rh_'));
        
        try
            mean_MIND = (lh{:,:}+rh{:,:})/2;
        catch
            lacking_region = setdiff(replace(rh.Properties.VariableNames,'rh_',''),replace(lh.Properties.VariableNames,'lh_',''));
            index_region = find(strcmp(rh.Properties.VariableNames,['rh_',lacking_region{:}]));
            lacking_region = ['lh_',lacking_region{:}];
            lh = [lh(:,1:index_region-1), array2table(nan(length(mean_MIND)-1,1),"VariableNames",{lacking_region}), lh(:,index_region:end)];
            lh = [lh(1:index_region-1,:); array2table(nan(1,length(mean_MIND)),"RowNames",{lacking_region},"VariableNames",lh.Properties.VariableNames); lh(index_region:end,:)];
            lh(:,index_region) = rh(:,index_region);
            lh(index_region,:) = rh(index_region,:);       
            mean_MIND = (lh{:,:}+rh{:,:})/2;

            if contains(lacking_region,'lh_')
                MIND = [MIND(:,1:index_region-1), array2table(nan(width(MIND),1),"VariableNames",{lacking_region}), MIND(:,index_region:end)];
                MIND = [MIND(1:index_region-1,:); array2table(nan(1,width(MIND)),"RowNames",{lacking_region},"VariableNames",MIND.Properties.VariableNames); MIND(index_region:end,:)];
                MIND(:,index_region) = MIND(:,width(mean_MIND)+index_region);
                MIND(index_region,:) = MIND(width(mean_MIND)+index_region,:);       
            else
                MIND = [MIND(:,1:width(mean_MIND)+index_region-1), array2table(nan(width(MIND),1),"VariableNames",{lacking_region}), MIND(:,width(mean_MIND)+index_region:end)];
                MIND = [MIND(1:width(mean_MIND)+index_region-1,:); array2table(nan(1,width(MIND)),"RowNames",{lacking_region},"VariableNames",MIND.Properties.VariableNames); MIND(width(mean_MIND)+index_region:end,:)];
                MIND(:,width(mean_MIND)+index_region) = MIND(:,index_region);
                MIND(width(mean_MIND)+index_region,:) = MIND(index_region,:); 
            end

        end
    
        % FEP 
        if str2double(regexp(files(i).name, '(\d+)_', 'tokens', 'once')) < 1000
    
            FEP_subjects = [FEP_subjects,replace(files(i).name,'.csv','')];

            % Edge
            edge_68_FEP{i_FEP,:} = MIND;
    
            % Degree
            degree_68_FEP(i_FEP,:) = mean(MIND{:,:});
                
    
            i_FEP = i_FEP + 1;
        
        % CN
        else
            CN_subjects = [CN_subjects,replace(files(i).name,'.csv','')];
    
            % Edge
            edge_68_CN{i_CN,:} = MIND;

            % Degree
            degree_68_CN(i_CN,:) = mean(MIND{:,:});
             
    
            i_CN = i_CN + 1;
        end
    end
    
    degree_68_FEP = array2table(degree_68_FEP,"RowNames",FEP_subjects,"VariableNames",MIND.Properties.VariableNames);
    degree_68_CN = array2table(degree_68_CN,"RowNames",CN_subjects,"VariableNames",MIND.Properties.VariableNames);
    
    edge_68_FEP = [FEP_subjects',edge_68_FEP];
    edge_68_CN = [CN_subjects',edge_68_CN];
    
    writetable(degree_68_FEP,[[location,write_degree],'\degree_68_FEP.csv'],'WriteRowNames',true);
    writetable(degree_68_CN,[[location,write_degree],'\degree_68_CN.csv'],'WriteRowNames',true);
    save([location,write_degree,'degree_68_FEP.mat'],'degree_68_FEP')
    save([location,write_degree,'degree_68_CN.mat'],'degree_68_CN')
    save([location,write_edges,'edge_68_FEP.mat'],'edge_68_FEP')
    save([location,write_edges,'edge_68_CN.mat'],'edge_68_CN')   

end


% EFFECT SIZES DEGREE
degree_68_FEP.Properties.VariableNames = replace(replace(degree_68_FEP.Properties.VariableNames,'lh','L'),'rh','R');
degree_68_FEP.Properties.VariableNames = cellfun(@(x) [x, '_grayavg'], degree_68_FEP.Properties.VariableNames, 'UniformOutput', false);
degree_68_CN.Properties.VariableNames = replace(replace(degree_68_CN.Properties.VariableNames,'lh','L'),'rh','R');
degree_68_CN.Properties.VariableNames = cellfun(@(x) [x, '_grayavg'], degree_68_CN.Properties.VariableNames, 'UniformOutput', false);

centiles = readtable([location,'\datasets\PAFIP\centiles\centiles_CorticalMeasuresENIGMA_GrayAvg.csv'],"ReadRowNames",true);
centiles = centiles(centiles.session == 1,:);
centiles_FEP = centiles(strcmp(centiles.dx,'PS'),:);

centiles_cognition = readtable([location,'\datasets\PAFIP\centiles\ClinicalAndCentilesInformationDataBase_11_12_2023.csv'],"ReadRowNames",true);
centiles_cognition = centiles_cognition(centiles_cognition.Assessment == 1 & ismember(centiles_cognition.Subject,centiles_FEP.participant),:);
centiles_cognition = [centiles_cognition(:,{'Subject','Global_Cognitive_Functioning_average'}), array2table(str2double(centiles_cognition{:,{'BPRS_Total','SANS_Total','SAPS_Total','CPZ_equivalent'}}),'VariableNames',{'BPRS_Total','SANS_Total','SAPS_Total','CPZ_equivalent'})];
centiles_cognition = centiles_cognition(degree_68_FEP.Properties.RowNames,:);

threshold = 'median'; % change to compute the 25% or the median effsizes
threshold = '25'; % change to compute the 25% or the median effsizes
if ~strcmp(threshold,'median')
    [~, idx_cognition] = sort(centiles_cognition.Global_Cognitive_Functioning_average);
    [~, idx_BPRS] = sort(centiles_cognition.BPRS_Total);
    [~, idx_SAPS] = sort(centiles_cognition.SAPS_Total);
    [~, idx_SANS] = sort(centiles_cognition.SANS_Total);

    threshold_effsizes_low_high = round(height(degree_68_FEP)*0.25); % 25% top- and bottom- scored
    
    num_participants = num2str(threshold_effsizes_low_high);
else
    num_participants = '';
end

for i = 1:width(degree_68_CN)
    effsizes_degree_68(:,i) = computeCohen_d(degree_68_FEP{:,i},degree_68_CN{:,i});

    effsizes_degree_low_cognition(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.Global_Cognitive_Functioning_average < median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),i},degree_68_CN{:,i});
    effsizes_degree_high_cognition(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.Global_Cognitive_Functioning_average >= median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),i},degree_68_CN{:,i});
    effsizes_degree_low_BPRS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.BPRS_Total < median(centiles_cognition.BPRS_Total,'omitnan'),i},degree_68_CN{:,i});
    effsizes_degree_high_BPRS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.BPRS_Total >= median(centiles_cognition.BPRS_Total,'omitnan'),i},degree_68_CN{:,i});
    effsizes_degree_low_SAPS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.SAPS_Total < median(centiles_cognition.SAPS_Total,'omitnan'),i},degree_68_CN{:,i});
    effsizes_degree_high_SAPS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.SAPS_Total >= median(centiles_cognition.SAPS_Total,'omitnan'),i},degree_68_CN{:,i});
    effsizes_degree_low_SANS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.SANS_Total < median(centiles_cognition.SANS_Total,'omitnan'),i},degree_68_CN{:,i});
    effsizes_degree_high_SANS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.SANS_Total >= median(centiles_cognition.SANS_Total,'omitnan'),i},degree_68_CN{:,i});


    if ~strcmp(threshold,'median')
        % 25% top- and bottom- scored
        effsizes_degree_low_high_cognition(:,i) = computeCohen_d(degree_68_FEP{idx_cognition(1:threshold_effsizes_low_high),i},degree_68_FEP{idx_cognition(length(idx_cognition):-1:(length(idx_cognition)-threshold_effsizes_low_high+1)),i});
        effsizes_degree_low_high_BPRS(:,i) = computeCohen_d(degree_68_FEP{idx_BPRS(length(idx_BPRS):-1:(length(idx_BPRS)-threshold_effsizes_low_high+1)),i},degree_68_FEP{idx_BPRS(1:threshold_effsizes_low_high),i});
        effsizes_degree_low_high_SAPS(:,i) = computeCohen_d(degree_68_FEP{idx_SAPS(length(idx_SAPS):-1:(length(idx_SAPS)-threshold_effsizes_low_high+1)),i},degree_68_FEP{idx_SAPS(1:threshold_effsizes_low_high),i});
        effsizes_degree_low_high_SANS(:,i) = computeCohen_d(degree_68_FEP{idx_SANS(length(idx_SANS):-1:(length(idx_SANS)-threshold_effsizes_low_high+1)),i},degree_68_FEP{idx_SANS(1:threshold_effsizes_low_high),i});
    
        [h,p_low_high_cognition] = ttest2(degree_68_FEP{idx_cognition(1:threshold_effsizes_low_high),i},degree_68_FEP{idx_cognition(length(idx_cognition):-1:(length(idx_cognition)-threshold_effsizes_low_high+1)),i});
        pval_effsizes_edges_68_low_high_cognition(:,i) = p_low_high_cognition;
        [h,p_low_high_BPRS] = ttest2(degree_68_FEP{idx_BPRS(length(idx_BPRS):-1:(length(idx_BPRS)-threshold_effsizes_low_high+1)),i},degree_68_FEP{idx_BPRS(1:threshold_effsizes_low_high),i});
        pval_effsizes_edges_68_low_high_BPRS(:,i) = p_low_high_BPRS;
        [h,p_low_high_SAPS] = ttest2(degree_68_FEP{idx_SAPS(length(idx_SAPS):-1:(length(idx_SAPS)-threshold_effsizes_low_high+1)),i},degree_68_FEP{idx_SAPS(1:threshold_effsizes_low_high),i});
        pval_effsizes_edges_68_low_high_SAPS(:,i) = p_low_high_SAPS;
        [h,p_low_high_SANS] = ttest2(degree_68_FEP{idx_SANS(length(idx_SANS):-1:(length(idx_SANS)-threshold_effsizes_low_high+1)),i},degree_68_FEP{idx_SANS(1:threshold_effsizes_low_high),i});
        pval_effsizes_edges_68_low_high_SANS(:,i) = p_low_high_SANS;

    else
        % median scored
        effsizes_degree_low_high_cognition(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.Global_Cognitive_Functioning_average < median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),i},degree_68_FEP{centiles_cognition.Global_Cognitive_Functioning_average >= median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),i});
        effsizes_degree_low_high_BPRS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.BPRS_Total >= median(centiles_cognition.BPRS_Total,'omitnan'),i},degree_68_FEP{centiles_cognition.BPRS_Total < median(centiles_cognition.BPRS_Total,'omitnan'),i});
        effsizes_degree_low_high_SAPS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.SAPS_Total >= median(centiles_cognition.SAPS_Total,'omitnan'),i},degree_68_FEP{centiles_cognition.SAPS_Total < median(centiles_cognition.SAPS_Total,'omitnan'),i});
        effsizes_degree_low_high_SANS(:,i) = computeCohen_d(degree_68_FEP{centiles_cognition.SANS_Total >= median(centiles_cognition.SANS_Total,'omitnan'),i},degree_68_FEP{centiles_cognition.SANS_Total < median(centiles_cognition.SANS_Total,'omitnan'),i});
    
        [h,p_low_high_cognition] = ttest2(degree_68_FEP{centiles_cognition.Global_Cognitive_Functioning_average < median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),i},degree_68_FEP{centiles_cognition.Global_Cognitive_Functioning_average >= median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),i});
        pval_effsizes_edges_68_low_high_cognition(:,i) = p_low_high_cognition;
        [h,p_low_high_BPRS] = ttest2(degree_68_FEP{centiles_cognition.BPRS_Total < median(centiles_cognition.BPRS_Total,'omitnan'),i},degree_68_FEP{centiles_cognition.BPRS_Total >= median(centiles_cognition.BPRS_Total,'omitnan'),i});
        pval_effsizes_edges_68_low_high_BPRS(:,i) = p_low_high_BPRS;
        [h,p_low_high_SAPS] = ttest2(degree_68_FEP{centiles_cognition.SAPS_Total < median(centiles_cognition.SAPS_Total,'omitnan'),i},degree_68_FEP{centiles_cognition.SAPS_Total >= median(centiles_cognition.SAPS_Total,'omitnan'),i});
        pval_effsizes_edges_68_low_high_SAPS(:,i) = p_low_high_SAPS;
        [h,p_low_high_SANS] = ttest2(degree_68_FEP{centiles_cognition.SANS_Total < median(centiles_cognition.SANS_Total,'omitnan'),i},degree_68_FEP{centiles_cognition.SANS_Total >= median(centiles_cognition.SANS_Total,'omitnan'),i});
        pval_effsizes_edges_68_low_high_SANS(:,i) = p_low_high_SANS;

    end
end
writetable(array2table(effsizes_degree_68,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68.csv'],'WriteRowNames',true);
save('effsizes_degree_68',"effsizes_degree_68")

writetable(array2table(effsizes_degree_low_cognition,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_low_cognition.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_high_cognition,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_high_cognition.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_low_BPRS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_low_BPRS.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_high_BPRS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_high_BPRS.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_low_SAPS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_low_SAPS.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_high_SAPS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_high_SAPS.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_low_SANS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_low_SANS.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_high_SANS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_high_SANS.csv'],'WriteRowNames',true);

pval_eff_degree_68_low_high_cognition = mafdr(pval_effsizes_edges_68_low_high_cognition,'BHFDR',true);
effsizes_degree_68_low_high_significant = effsizes_degree_low_high_cognition;
effsizes_degree_68_low_high_significant(pval_eff_degree_68_low_high_cognition >= 0.05) = NaN;
effsizes_degree_68_low_high_significant = array2table(effsizes_degree_68_low_high_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
writetable(effsizes_degree_68_low_high_significant,[[location,write_degree],'\effsizes_degree_68_low_high_',num_participants,'_cognition_ttest_significant.csv'],'WriteRowNames',true);

pval_eff_degree_68_low_high_BPRS = mafdr(pval_effsizes_edges_68_low_high_BPRS,'BHFDR',true);
effsizes_degree_68_low_high_significant = effsizes_degree_low_high_BPRS;
effsizes_degree_68_low_high_significant(pval_effsizes_edges_68_low_high_BPRS >= 0.05) = NaN;
effsizes_degree_68_low_high_significant = array2table(effsizes_degree_68_low_high_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
writetable(effsizes_degree_68_low_high_significant,[[location,write_degree],'\effsizes_degree_68_low_high_',num_participants,'_BPRS_ttest_significant.csv'],'WriteRowNames',true);

pval_eff_degree_68_low_high_SAPS = mafdr(pval_effsizes_edges_68_low_high_SAPS,'BHFDR',true);
effsizes_degree_68_low_high_significant = effsizes_degree_low_high_SAPS;
effsizes_degree_68_low_high_significant(pval_effsizes_edges_68_low_high_SAPS >= 0.05) = NaN;
effsizes_degree_68_low_high_significant = array2table(effsizes_degree_68_low_high_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
writetable(effsizes_degree_68_low_high_significant,[[location,write_degree],'\effsizes_degree_68_low_high_',num_participants,'_SAPS_ttest_significant.csv'],'WriteRowNames',true);

pval_eff_degree_68_low_high_SANS = mafdr(pval_effsizes_edges_68_low_high_SANS,'BHFDR',true);
effsizes_degree_68_low_high_significant = effsizes_degree_low_high_SANS;
effsizes_degree_68_low_high_significant(pval_effsizes_edges_68_low_high_SANS >= 0.05) = NaN;
effsizes_degree_68_low_high_significant = array2table(effsizes_degree_68_low_high_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
writetable(effsizes_degree_68_low_high_significant,[[location,write_degree],'\effsizes_degree_68_low_high_',num_participants,'_SANS_ttest_significant.csv'],'WriteRowNames',true);

writetable(array2table(effsizes_degree_low_high_cognition,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_low_high_',num_participants,'_cognition.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_low_high_BPRS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_low_high_',num_participants,'_BPRS.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_low_high_SAPS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_low_high_',num_participants,'_SAPS.csv'],'WriteRowNames',true);
writetable(array2table(effsizes_degree_low_high_SANS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames","effsizes"),[[location,write_degree],'\effsizes_degree_68_low_high_',num_participants,'_SANS.csv'],'WriteRowNames',true);

nperm = 1000;
if ~exist([location,write_degree,'effsizes_degree_68_significant.csv'])
    for iperm = 1:nperm
        [selection_dx_1,selection_dx_2] = mix_dx(degree_68_FEP,degree_68_CN);

        [selection_dx_1_low_cognition,selection_dx_2_low_cognition] = mix_dx(degree_68_FEP(centiles_cognition.Global_Cognitive_Functioning_average < median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),:),degree_68_CN);
        [selection_dx_1_high_cognition,selection_dx_2_high_cognition] = mix_dx(degree_68_FEP(centiles_cognition.Global_Cognitive_Functioning_average >= median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),:),degree_68_CN);    
        [selection_dx_1_low_BPRS,selection_dx_2_low_BPRS] = mix_dx(degree_68_FEP(centiles_cognition.BPRS_Total < median(centiles_cognition.BPRS_Total,'omitnan'),:),degree_68_CN);
        [selection_dx_1_high_BPRS,selection_dx_2_high_BPRS] = mix_dx(degree_68_FEP(centiles_cognition.BPRS_Total >= median(centiles_cognition.BPRS_Total,'omitnan'),:),degree_68_CN);    
        [selection_dx_1_low_SAPS,selection_dx_2_low_SAPS] = mix_dx(degree_68_FEP(centiles_cognition.SAPS_Total < median(centiles_cognition.SAPS_Total,'omitnan'),:),degree_68_CN);
        [selection_dx_1_high_SAPS,selection_dx_2_high_SAPS] = mix_dx(degree_68_FEP(centiles_cognition.SAPS_Total >= median(centiles_cognition.SAPS_Total,'omitnan'),:),degree_68_CN);    
        [selection_dx_1_low_SANS,selection_dx_2_low_SANS] = mix_dx(degree_68_FEP(centiles_cognition.SANS_Total < median(centiles_cognition.SANS_Total,'omitnan'),:),degree_68_CN);
        [selection_dx_1_high_SANS,selection_dx_2_high_SANS] = mix_dx(degree_68_FEP(centiles_cognition.SANS_Total >= median(centiles_cognition.SANS_Total,'omitnan'),:),degree_68_CN);    

        for ir = 1:width(degree_68_CN)
            effsize_degree_perm_68(ir) = computeCohen_d(selection_dx_1(:,ir),selection_dx_2(:,ir));

            effsize_degree_perm_68_low_cognition(ir) = computeCohen_d(selection_dx_1_low_cognition(:,ir),selection_dx_2_low_cognition(:,ir));
            effsize_degree_perm_68_high_cognition(ir) = computeCohen_d(selection_dx_1_high_cognition(:,ir),selection_dx_2_high_cognition(:,ir));
            effsize_degree_perm_68_low_BPRS(ir) = computeCohen_d(selection_dx_1_low_BPRS(:,ir),selection_dx_2_low_BPRS(:,ir));
            effsize_degree_perm_68_high_BPRS(ir) = computeCohen_d(selection_dx_1_high_BPRS(:,ir),selection_dx_2_high_BPRS(:,ir));
            effsize_degree_perm_68_low_SAPS(ir) = computeCohen_d(selection_dx_1_low_SAPS(:,ir),selection_dx_2_low_SAPS(:,ir));
            effsize_degree_perm_68_high_SAPS(ir) = computeCohen_d(selection_dx_1_high_SAPS(:,ir),selection_dx_2_high_SAPS(:,ir));
            effsize_degree_perm_68_low_SANS(ir) = computeCohen_d(selection_dx_1_low_SANS(:,ir),selection_dx_2_low_SANS(:,ir));
            effsize_degree_perm_68_high_SANS(ir) = computeCohen_d(selection_dx_1_high_SANS(:,ir),selection_dx_2_high_SANS(:,ir));

        end
        effsizes_degree_perm_68(iperm,:) = effsize_degree_perm_68;

        effsizes_degree_perm_68_low_cognition(iperm,:) = effsize_degree_perm_68_low_cognition;
        effsizes_degree_perm_68_high_cognition(iperm,:) = effsize_degree_perm_68_high_cognition;
        effsizes_degree_perm_68_low_BPRS(iperm,:) = effsize_degree_perm_68_low_BPRS;
        effsizes_degree_perm_68_high_BPRS(iperm,:) = effsize_degree_perm_68_high_BPRS;
        effsizes_degree_perm_68_low_SAPS(iperm,:) = effsize_degree_perm_68_low_SAPS;
        effsizes_degree_perm_68_high_SAPS(iperm,:) = effsize_degree_perm_68_high_SAPS;
        effsizes_degree_perm_68_low_SANS(iperm,:) = effsize_degree_perm_68_low_SANS;
        effsizes_degree_perm_68_high_SANS(iperm,:) = effsize_degree_perm_68_high_SANS;


    end
    
    % Test if it explains more variance that expected by chance
    pval_eff_degree_68 = sum(abs(effsizes_degree_perm_68) > abs(effsizes_degree_68)) / nperm;
    pval_eff_degree_68 = mafdr(pval_eff_degree_68,'BHFDR',true);
    effsizes_degree_68_significant = effsizes_degree_68;
    effsizes_degree_68_significant(pval_eff_degree_68 >= 0.05) = NaN;
    effsizes_degree_68_significant = array2table(effsizes_degree_68_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_significant,[[location,write_degree],'\effsizes_degree_68_significant.csv'],'WriteRowNames',true);

    pval_eff_degree_68_low = sum(abs(effsizes_degree_perm_68_low_cognition) > abs(effsizes_degree_low_cognition)) / nperm;
    pval_eff_degree_68_low = mafdr(pval_eff_degree_68_low,'BHFDR',true);
    effsizes_degree_68_low_significant = effsizes_degree_low_cognition;
    effsizes_degree_68_low_significant(pval_eff_degree_68_low >= 0.05) = NaN;
    effsizes_degree_68_low_significant = array2table(effsizes_degree_68_low_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_low_significant,[[location,write_degree],'\effsizes_degree_68_low_cognition_significant.csv'],'WriteRowNames',true);
    
    pval_eff_degree_68_high = sum(abs(effsizes_degree_perm_68_high_cognition) > abs(effsizes_degree_high_cognition)) / nperm;
    pval_eff_degree_68_high = mafdr(pval_eff_degree_68_high,'BHFDR',true);
    effsizes_degree_68_high_significant = effsizes_degree_high_cognition;
    effsizes_degree_68_high_significant(pval_eff_degree_68_high >= 0.05) = NaN;
    effsizes_degree_68_high_significant = array2table(effsizes_degree_68_high_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_high_significant,[[location,write_degree],'\effsizes_degree_68_high_cognition_significant.csv'],'WriteRowNames',true);
        
    pval_eff_degree_68_low = sum(abs(effsizes_degree_perm_68_low_BPRS) > abs(effsizes_degree_low_BPRS)) / nperm;
    pval_eff_degree_68_low = mafdr(pval_eff_degree_68_low,'BHFDR',true);
    effsizes_degree_68_low_significant = effsizes_degree_low_BPRS;
    effsizes_degree_68_low_significant(pval_eff_degree_68_low >= 0.05) = NaN;
    effsizes_degree_68_low_significant = array2table(effsizes_degree_68_low_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_low_significant,[[location,write_degree],'\effsizes_degree_68_low_BPRS_significant.csv'],'WriteRowNames',true);
    
    pval_eff_degree_68_high = sum(abs(effsizes_degree_perm_68_high_BPRS) > abs(effsizes_degree_high_BPRS)) / nperm;
    pval_eff_degree_68_high = mafdr(pval_eff_degree_68_high,'BHFDR',true);
    effsizes_degree_68_high_significant = effsizes_degree_high_BPRS;
    effsizes_degree_68_high_significant(pval_eff_degree_68_high >= 0.05) = NaN;
    effsizes_degree_68_high_significant = array2table(effsizes_degree_68_high_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_high_significant,[[location,write_degree],'\effsizes_degree_68_high_BPRS_significant.csv'],'WriteRowNames',true);
    
    pval_eff_degree_68_low = sum(abs(effsizes_degree_perm_68_low_SAPS) > abs(effsizes_degree_low_SAPS)) / nperm;
    pval_eff_degree_68_low = mafdr(pval_eff_degree_68_low,'BHFDR',true);
    effsizes_degree_68_low_significant = effsizes_degree_low_SAPS;
    effsizes_degree_68_low_significant(pval_eff_degree_68_low >= 0.05) = NaN;
    effsizes_degree_68_low_significant = array2table(effsizes_degree_68_low_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_low_significant,[[location,write_degree],'\effsizes_degree_68_low_SAPS_significant.csv'],'WriteRowNames',true);
    
    pval_eff_degree_68_high = sum(abs(effsizes_degree_perm_68_high_SAPS) > abs(effsizes_degree_high_SAPS)) / nperm;
    pval_eff_degree_68_high = mafdr(pval_eff_degree_68_high,'BHFDR',true);
    effsizes_degree_68_high_significant = effsizes_degree_high_SAPS;
    effsizes_degree_68_high_significant(pval_eff_degree_68_high >= 0.05) = NaN;
    effsizes_degree_68_high_significant = array2table(effsizes_degree_68_high_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_high_significant,[[location,write_degree],'\effsizes_degree_68_high_SAPS_significant.csv'],'WriteRowNames',true);
    
    pval_eff_degree_68_low = sum(abs(effsizes_degree_perm_68_low_SANS) > abs(effsizes_degree_low_SANS)) / nperm;
    pval_eff_degree_68_low = mafdr(pval_eff_degree_68_low,'BHFDR',true);
    effsizes_degree_68_low_significant = effsizes_degree_low_SANS;
    effsizes_degree_68_low_significant(pval_eff_degree_68_low >= 0.05) = NaN;
    effsizes_degree_68_low_significant = array2table(effsizes_degree_68_low_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_low_significant,[[location,write_degree],'\effsizes_degree_68_low_SANS_significant.csv'],'WriteRowNames',true);
    
    pval_eff_degree_68_high = sum(abs(effsizes_degree_perm_68_high_SANS) > abs(effsizes_degree_high_SANS)) / nperm;
    pval_eff_degree_68_high = mafdr(pval_eff_degree_68_high,'BHFDR',true);
    effsizes_degree_68_high_significant = effsizes_degree_high_SANS;
    effsizes_degree_68_high_significant(pval_eff_degree_68_high >= 0.05) = NaN;
    effsizes_degree_68_high_significant = array2table(effsizes_degree_68_high_significant,'RowNames',"effsizes",'VariableNames',degree_68_CN.Properties.VariableNames);
    writetable(effsizes_degree_68_high_significant,[[location,write_degree],'\effsizes_degree_68_high_SANS_significant.csv'],'WriteRowNames',true);


else
    effsizes_degree_68_significant = readtable("effsizes_degree_68_significant.csv","ReadRowNames",true);
end



% EFFSIZES EDGES
if  ~exist([location,write_edges,'effsizes_edges_68.csv']) || ~exist([location,write_edges,'effsizes_edges_68_low_cognition.csv'])
    edge_68_FEP_low_cognition = edge_68_FEP(centiles_cognition.Global_Cognitive_Functioning_average < median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),:);
    edge_68_FEP_high_cognition = edge_68_FEP(centiles_cognition.Global_Cognitive_Functioning_average >= median(centiles_cognition.Global_Cognitive_Functioning_average,'omitnan'),:);
    edge_68_FEP_low_BPRS = edge_68_FEP(centiles_cognition.BPRS_Total < median(centiles_cognition.BPRS_Total,'omitnan'),:);
    edge_68_FEP_high_BPRS = edge_68_FEP(centiles_cognition.BPRS_Total >= median(centiles_cognition.BPRS_Total,'omitnan'),:);
    edge_68_FEP_low_SAPS = edge_68_FEP(centiles_cognition.SAPS_Total < median(centiles_cognition.SAPS_Total,'omitnan'),:);
    edge_68_FEP_high_SAPS = edge_68_FEP(centiles_cognition.SAPS_Total >= median(centiles_cognition.SAPS_Total,'omitnan'),:);
    edge_68_FEP_low_SANS = edge_68_FEP(centiles_cognition.SANS_Total < median(centiles_cognition.SANS_Total,'omitnan'),:);
    edge_68_FEP_high_SANS = edge_68_FEP(centiles_cognition.SANS_Total>= median(centiles_cognition.SANS_Total,'omitnan'),:);

    for i = 1:width(edge_68_CN{1,2})
        for  j = 1:width(edge_68_CN{1,2})
            for k = 1:length(edge_68_CN)
                edge_subject_CN = table2array(edge_68_CN{k,2});
                edge_region_subject_CN(k) = edge_subject_CN(i,j);
            end
            for l = 1:length(edge_68_FEP)
                edge_subject_FEP = table2array(edge_68_FEP{l,2});
                edge_region_subject_FEP(l) = edge_subject_FEP(i,j);
            end
            for l = 1:length(edge_68_FEP_low_cognition)
                edge_subject_FEP_low_cognition = table2array(edge_68_FEP_low_cognition{l,2});
                edge_region_subject_FEP_low_cognition(l) = edge_subject_FEP_low_cognition(i,j);
            end
            for l = 1:length(edge_68_FEP_high_cognition)
                edge_subject_FEP_high_cognition = table2array(edge_68_FEP_high_cognition{l,2});
                edge_region_subject_FEP_high_cognition(l) = edge_subject_FEP_high_cognition(i,j);
            end
            for l = 1:length(edge_68_FEP_low_BPRS)
                edge_subject_FEP_low_BPRS = table2array(edge_68_FEP_low_BPRS{l,2});
                edge_region_subject_FEP_low_BPRS(l) = edge_subject_FEP_low_BPRS(i,j);
            end
            for l = 1:length(edge_68_FEP_high_BPRS)
                edge_subject_FEP_high_BPRS = table2array(edge_68_FEP_high_BPRS{l,2});
                edge_region_subject_FEP_high_BPRS(l) = edge_subject_FEP_high_BPRS(i,j);
            end
            for l = 1:length(edge_68_FEP_low_SAPS)
                edge_subject_FEP_low_SAPS = table2array(edge_68_FEP_low_SAPS{l,2});
                edge_region_subject_FEP_low_SAPS(l) = edge_subject_FEP_low_SAPS(i,j);
            end
            for l = 1:length(edge_68_FEP_high_SAPS)
                edge_subject_FEP_high_SAPS = table2array(edge_68_FEP_high_SAPS{l,2});
                edge_region_subject_FEP_high_SAPS(l) = edge_subject_FEP_high_SAPS(i,j);
            end
            for l = 1:length(edge_68_FEP_low_SANS)
                edge_subject_FEP_low_SANS = table2array(edge_68_FEP_low_SANS{l,2});
                edge_region_subject_FEP_low_SANS(l) = edge_subject_FEP_low_SANS(i,j);
            end
            for l = 1:length(edge_68_FEP_high_SANS)
                edge_subject_FEP_high_SANS = table2array(edge_68_FEP_high_SANS{l,2});
                edge_region_subject_FEP_high_SANS(l) = edge_subject_FEP_high_SANS(i,j);
            end
            effsizes_edges_68(i,j) = computeCohen_d(edge_region_subject_FEP,edge_region_subject_CN);
            [h,p] = ttest2(edge_region_subject_FEP,edge_region_subject_CN);
            pval_effsizes_edges_68(i,j) = p;

            effsizes_edges_68_low_cognition(i,j) = computeCohen_d(edge_region_subject_FEP_low_cognition,edge_region_subject_CN);
            effsizes_edges_68_high_cognition(i,j) = computeCohen_d(edge_region_subject_FEP_high_cognition,edge_region_subject_CN);
            effsizes_edges_68_low_BPRS(i,j) = computeCohen_d(edge_region_subject_FEP_low_BPRS,edge_region_subject_CN);
            effsizes_edges_68_high_BPRS(i,j) = computeCohen_d(edge_region_subject_FEP_high_BPRS,edge_region_subject_CN);
            effsizes_edges_68_low_SAPS(i,j) = computeCohen_d(edge_region_subject_FEP_low_SAPS,edge_region_subject_CN);
            effsizes_edges_68_high_SAPS(i,j) = computeCohen_d(edge_region_subject_FEP_high_SAPS,edge_region_subject_CN);
            effsizes_edges_68_low_SANS(i,j) = computeCohen_d(edge_region_subject_FEP_low_SANS,edge_region_subject_CN);
            effsizes_edges_68_high_SANS(i,j) = computeCohen_d(edge_region_subject_FEP_high_SANS,edge_region_subject_CN);
            
            [h,p_low] = ttest2(edge_region_subject_FEP_low_cognition,edge_region_subject_CN);
            [h,p_high] = ttest2(edge_region_subject_FEP_high_cognition,edge_region_subject_CN);
            pval_effsizes_edges_68_low_cognition(i,j) = p_low;
            pval_effsizes_edges_68_high_cognition(i,j) = p_high;

            [h,p_low] = ttest2(edge_region_subject_FEP_low_BPRS,edge_region_subject_CN);
            [h,p_high] = ttest2(edge_region_subject_FEP_high_BPRS,edge_region_subject_CN);
            pval_effsizes_edges_68_low_BPRS(i,j) = p_low;
            pval_effsizes_edges_68_high_BPRS(i,j) = p_high;

            [h,p_low] = ttest2(edge_region_subject_FEP_low_SAPS,edge_region_subject_CN);
            [h,p_high] = ttest2(edge_region_subject_FEP_high_SAPS,edge_region_subject_CN);
            pval_effsizes_edges_68_low_SAPS(i,j) = p_low;
            pval_effsizes_edges_68_high_SAPS(i,j) = p_high;

            [h,p_low] = ttest2(edge_region_subject_FEP_low_SANS,edge_region_subject_CN);
            [h,p_high] = ttest2(edge_region_subject_FEP_high_SANS,edge_region_subject_CN);
            pval_effsizes_edges_68_low_SANS(i,j) = p_low;
            pval_effsizes_edges_68_high_SANS(i,j) = p_high;
        end
    end
    writetable(array2table(effsizes_edges_68,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68.csv'],'WriteRowNames',true);
    save([location,write_edges,'effsizes_edges_68'],"effsizes_edges_68")
    save([location,write_edges,'pval_effsizes_edges_68'],"pval_effsizes_edges_68")

    writetable(array2table(effsizes_edges_68_low_cognition,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68_low_cognition.csv'],'WriteRowNames',true);
    writetable(array2table(effsizes_edges_68_high_cognition,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68_high_cognition.csv'],'WriteRowNames',true);
    save([location,write_edges,'effsizes_edges_68_low_cognition'],"effsizes_edges_68_low_cognition")
    save([location,write_edges,'pval_effsizes_edges_68_low_cognition'],"pval_effsizes_edges_68_low_cognition")
    save([location,write_edges,'effsizes_edges_68_high_cognition'],"effsizes_edges_68_high_cognition")
    save([location,write_edges,'pval_effsizes_edges_68_high_cognition'],"pval_effsizes_edges_68_high_cognition")

    writetable(array2table(effsizes_edges_68_low_BPRS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68_low_BPRS.csv'],'WriteRowNames',true);
    writetable(array2table(effsizes_edges_68_high_BPRS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68_high_BPRS.csv'],'WriteRowNames',true);
    save([location,write_edges,'effsizes_edges_68_low_BPRS'],"effsizes_edges_68_low_BPRS")
    save([location,write_edges,'pval_effsizes_edges_68_low_BPRS'],"pval_effsizes_edges_68_low_BPRS")
    save([location,write_edges,'effsizes_edges_68_high_BPRS'],"effsizes_edges_68_high_BPRS")
    save([location,write_edges,'pval_effsizes_edges_68_high_BPRS'],"pval_effsizes_edges_68_high_BPRS")

    writetable(array2table(effsizes_edges_68_low_SAPS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68_low_SAPS.csv'],'WriteRowNames',true);
    writetable(array2table(effsizes_edges_68_high_SAPS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68_high_SAPS.csv'],'WriteRowNames',true);
    save([location,write_edges,'effsizes_edges_68_low_SAPS'],"effsizes_edges_68_low_SAPS")
    save([location,write_edges,'pval_effsizes_edges_68_low_SAPS'],"pval_effsizes_edges_68_low_SAPS")
    save([location,write_edges,'effsizes_edges_68_high_SAPS'],"effsizes_edges_68_high_SAPS")
    save([location,write_edges,'pval_effsizes_edges_68_high_SAPS'],"pval_effsizes_edges_68_high_SAPS")

    writetable(array2table(effsizes_edges_68_low_SANS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68_low_SANS.csv'],'WriteRowNames',true);
    writetable(array2table(effsizes_edges_68_high_SANS,"VariableNames",degree_68_CN.Properties.VariableNames,"RowNames",degree_68_CN.Properties.VariableNames),[location,write_edges,'effsizes_edges_68_high_SANS.csv'],'WriteRowNames',true);
    save([location,write_edges,'effsizes_edges_68_low_SANS'],"effsizes_edges_68_low_SANS")
    save([location,write_edges,'pval_effsizes_edges_68_low_SANS'],"pval_effsizes_edges_68_low_SANS")
    save([location,write_edges,'effsizes_edges_68_high_SANS'],"effsizes_edges_68_high_SANS")
    save([location,write_edges,'pval_effsizes_edges_68_high_SANS'],"pval_effsizes_edges_68_high_SANS")
else
    load("effsizes_edges_68.mat")
    load("pval_effsizes_edges_68.mat")

    load("effsizes_edges_68_low_cognition.mat")
    load("pval_effsizes_edges_68_low_cognition.mat")
    load("effsizes_edges_68_high_cognition.mat")
    load("pval_effsizes_edges_68_high_cognition.mat")

    load("effsizes_edges_68_low_BPRS.mat")
    load("pval_effsizes_edges_68_low_BPRS.mat")
    load("effsizes_edges_68_high_BPRS.mat")
    load("pval_effsizes_edges_68_high_BPRS.mat")

    load("effsizes_edges_68_low_SAPS.mat")
    load("pval_effsizes_edges_68_low_SAPS.mat")
    load("effsizes_edges_68_high_SAPS.mat")
    load("pval_effsizes_edges_68_high_SAPS.mat")

    load("effsizes_edges_68_low_SANS.mat")
    load("pval_effsizes_edges_68_low_SANS.mat")
    load("effsizes_edges_68_high_SANS.mat")
    load("pval_effsizes_edges_68_high_SANS.mat")

end

pval_effsizes_col = pval_effsizes_edges_68(find(triu(ones(size(pval_effsizes_edges_68)),1)));
pval_effsizes_corrected_col = mafdr(pval_effsizes_col,'BHFDR',true);
pval_effsizes_edges_corrected = zeros(size(pval_effsizes_edges_68));
upper_triangular_index = find(triu(ones(size(pval_effsizes_edges_68)), 1));
pval_effsizes_edges_corrected(upper_triangular_index) = pval_effsizes_corrected_col;
pval_effsizes_edges_corrected = pval_effsizes_edges_corrected + pval_effsizes_edges_corrected';
effsizes_edges_68_significant = effsizes_edges_68;
effsizes_edges_68_significant(pval_effsizes_edges_corrected >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected'],"pval_effsizes_edges_corrected")

pval_effsizes_low_col = pval_effsizes_edges_68_low_cognition(find(triu(ones(size(pval_effsizes_edges_68_low_cognition)),1)));
pval_effsizes_corrected_low_col = mafdr(pval_effsizes_low_col,'BHFDR',true);
pval_effsizes_edges_corrected_low = zeros(size(pval_effsizes_edges_68_low_cognition));
upper_triangular_index_low = find(triu(ones(size(pval_effsizes_edges_68_low_cognition)), 1));
pval_effsizes_edges_corrected_low(upper_triangular_index_low) = pval_effsizes_corrected_low_col;
pval_effsizes_edges_corrected_low = pval_effsizes_edges_corrected_low + pval_effsizes_edges_corrected_low';
effsizes_edges_68_low_cognition_significant = effsizes_edges_68_low_cognition;
effsizes_edges_68_low_cognition_significant(pval_effsizes_edges_corrected_low >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected_low_cognition'],"pval_effsizes_edges_corrected_low")

pval_effsizes_high_col = pval_effsizes_edges_68_high_cognition(find(triu(ones(size(pval_effsizes_edges_68_high_cognition)),1)));
pval_effsizes_corrected_high_col = mafdr(pval_effsizes_high_col,'BHFDR',true);
pval_effsizes_edges_corrected_high = zeros(size(pval_effsizes_edges_68_high_cognition));
upper_triangular_index_high = find(triu(ones(size(pval_effsizes_edges_68_high_cognition)), 1));
pval_effsizes_edges_corrected_high(upper_triangular_index_high) = pval_effsizes_corrected_high_col;
pval_effsizes_edges_corrected_high = pval_effsizes_edges_corrected_high + pval_effsizes_edges_corrected_high';
effsizes_edges_68_high_cognition_significant = effsizes_edges_68_high_cognition;
effsizes_edges_68_high_cognition_significant(pval_effsizes_edges_corrected_high >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected_high_cognition'],"pval_effsizes_edges_corrected_high")

pval_effsizes_low_col = pval_effsizes_edges_68_low_BPRS(find(triu(ones(size(pval_effsizes_edges_68_low_BPRS)),1)));
pval_effsizes_corrected_low_col = mafdr(pval_effsizes_low_col,'BHFDR',true);
pval_effsizes_edges_corrected_low = zeros(size(pval_effsizes_edges_68_low_BPRS));
upper_triangular_index_low = find(triu(ones(size(pval_effsizes_edges_68_low_BPRS)), 1));
pval_effsizes_edges_corrected_low(upper_triangular_index_low) = pval_effsizes_corrected_low_col;
pval_effsizes_edges_corrected_low = pval_effsizes_edges_corrected_low + pval_effsizes_edges_corrected_low';
effsizes_edges_68_low_BPRS_significant = effsizes_edges_68_low_BPRS;
effsizes_edges_68_low_BPRS_significant(pval_effsizes_edges_corrected_low >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected_low_BPRS'],"pval_effsizes_edges_corrected_low")

pval_effsizes_high_col = pval_effsizes_edges_68_high_BPRS(find(triu(ones(size(pval_effsizes_edges_68_high_BPRS)),1)));
pval_effsizes_corrected_high_col = mafdr(pval_effsizes_high_col,'BHFDR',true);
pval_effsizes_edges_corrected_high = zeros(size(pval_effsizes_edges_68_high_BPRS));
upper_triangular_index_high = find(triu(ones(size(pval_effsizes_edges_68_high_BPRS)), 1));
pval_effsizes_edges_corrected_high(upper_triangular_index_high) = pval_effsizes_corrected_high_col;
pval_effsizes_edges_corrected_high = pval_effsizes_edges_corrected_high + pval_effsizes_edges_corrected_high';
effsizes_edges_68_high_BPRS_significant = effsizes_edges_68_high_BPRS;
effsizes_edges_68_high_BPRS_significant(pval_effsizes_edges_corrected_high >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected_high_BPRS'],"pval_effsizes_edges_corrected_high")

pval_effsizes_low_col = pval_effsizes_edges_68_low_SAPS(find(triu(ones(size(pval_effsizes_edges_68_low_SAPS)),1)));
pval_effsizes_corrected_low_col = mafdr(pval_effsizes_low_col,'BHFDR',true);
pval_effsizes_edges_corrected_low = zeros(size(pval_effsizes_edges_68_low_SAPS));
upper_triangular_index_low = find(triu(ones(size(pval_effsizes_edges_68_low_SAPS)), 1));
pval_effsizes_edges_corrected_low(upper_triangular_index_low) = pval_effsizes_corrected_low_col;
pval_effsizes_edges_corrected_low = pval_effsizes_edges_corrected_low + pval_effsizes_edges_corrected_low';
effsizes_edges_68_low_SAPS_significant = effsizes_edges_68_low_SAPS;
effsizes_edges_68_low_SAPS_significant(pval_effsizes_edges_corrected_low >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected_low_SAPS'],"pval_effsizes_edges_corrected_low")

pval_effsizes_high_col = pval_effsizes_edges_68_high_SAPS(find(triu(ones(size(pval_effsizes_edges_68_high_SAPS)),1)));
pval_effsizes_corrected_high_col = mafdr(pval_effsizes_high_col,'BHFDR',true);
pval_effsizes_edges_corrected_high = zeros(size(pval_effsizes_edges_68_high_SAPS));
upper_triangular_index_high = find(triu(ones(size(pval_effsizes_edges_68_high_SAPS)), 1));
pval_effsizes_edges_corrected_high(upper_triangular_index_high) = pval_effsizes_corrected_high_col;
pval_effsizes_edges_corrected_high = pval_effsizes_edges_corrected_high + pval_effsizes_edges_corrected_high';
effsizes_edges_68_high_SAPS_significant = effsizes_edges_68_high_SAPS;
effsizes_edges_68_high_SAPS_significant(pval_effsizes_edges_corrected_high >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected_high_SAPS'],"pval_effsizes_edges_corrected_high")

pval_effsizes_low_col = pval_effsizes_edges_68_low_SANS(find(triu(ones(size(pval_effsizes_edges_68_low_SANS)),1)));
pval_effsizes_corrected_low_col = mafdr(pval_effsizes_low_col,'BHFDR',true);
pval_effsizes_edges_corrected_low = zeros(size(pval_effsizes_edges_68_low_SANS));
upper_triangular_index_low = find(triu(ones(size(pval_effsizes_edges_68_low_SANS)), 1));
pval_effsizes_edges_corrected_low(upper_triangular_index_low) = pval_effsizes_corrected_low_col;
pval_effsizes_edges_corrected_low = pval_effsizes_edges_corrected_low + pval_effsizes_edges_corrected_low';
effsizes_edges_68_low_SANS_significant = effsizes_edges_68_low_SANS;
effsizes_edges_68_low_SANS_significant(pval_effsizes_edges_corrected_low >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected_low_SANS'],"pval_effsizes_edges_corrected_low")

pval_effsizes_high_col = pval_effsizes_edges_68_high_SANS(find(triu(ones(size(pval_effsizes_edges_68_high_SANS)),1)));
pval_effsizes_corrected_high_col = mafdr(pval_effsizes_high_col,'BHFDR',true);
pval_effsizes_edges_corrected_high = zeros(size(pval_effsizes_edges_68_high_SANS));
upper_triangular_index_high = find(triu(ones(size(pval_effsizes_edges_68_high_SANS)), 1));
pval_effsizes_edges_corrected_high(upper_triangular_index_high) = pval_effsizes_corrected_high_col;
pval_effsizes_edges_corrected_high = pval_effsizes_edges_corrected_high + pval_effsizes_edges_corrected_high';
effsizes_edges_68_high_SANS_significant = effsizes_edges_68_high_SANS;
effsizes_edges_68_high_SANS_significant(pval_effsizes_edges_corrected_high >= 0.05) = NaN;
save([location,write_edges,'pval_effsizes_edges_corrected_high_SANS'],"pval_effsizes_edges_corrected_high")

effsizes_edges_68_significant = array2table(effsizes_edges_68_significant,"RowNames",effsizes_degree_68_significant.Properties.VariableNames,"VariableNames",effsizes_degree_68_significant.Properties.VariableNames);
effsizes_edges_68_significant.Properties.VariableNames(contains(effsizes_edges_68_significant.Properties.VariableNames,'entoRinal')) = {'L_entorhinal_grayavg','R_entorhinal_grayavg'};
effsizes_edges_68_significant.Properties.RowNames(contains(effsizes_edges_68_significant.Properties.RowNames,'entoRinal')) = {'L_entorhinal_grayavg','R_entorhinal_grayavg'};


% plot by lobes
lobes.frontal = {'superiorfrontal','rostralmiddlefrontal','caudalmiddlefrontal','parsopercularis','parstriangularis', 'parsorbitalis'...
    'lateralorbitofrontal','medialorbitofrontal','precentral','paracentral','frontalpole'};
lobes.parietal = {'superiorparietal','inferiorparietal','supramarginal','postcentral','precuneus'};
lobes.temporal = {'superiortemporal', 'middletemporal','inferiortemporal','bankssts','fusiform','transversetemporal','entorhinal','temporalpole','parahippocampal'};
lobes.occipital = {'lateraloccipital','lingual','cuneus','pericalcarine'};
lobes.cingulate = {'rostralanteriorcingulate','caudalanteriorcingulate','posteriorcingulate','isthmuscingulate'};
lobes.insula = {'insula'};
mesulam = readtable([location,'Molecular\parcellations\Transform_aparc_TO_mesulam.csv']);
dk = array2table(mesulam.label1,"VariableNames",{'label1'});
dk(ismember(dk.label1,lobes.frontal),"label2") = {'frontal'};
dk(ismember(dk.label1,lobes.parietal),"label2") = {'parietal'};
dk(ismember(dk.label1,lobes.temporal),"label2") = {'temporal'};
dk(ismember(dk.label1,lobes.occipital),"label2") = {'occipital'};
dk(ismember(dk.label1,lobes.cingulate),"label2") = {'cingulate'};
dk(ismember(dk.label1,lobes.insula),"label2") = {'insula'};

dk = dk(~(strcmp(dk.label1,'no label')|strcmp(dk.label1,'unknown')),:);
atlas = dk;
regions = unique(atlas.label2,'stable');
[~,idx] = ismember(atlas.label2,regions);
atlas.idx = idx;
atlas_ordered = sortrows(atlas,'idx');
regions_ordered = unique(atlas_ordered.label1,'stable');
for i = 1:length(regions_ordered)
    regions_ordered_68(2*i-1:2*i) = {['L_',regions_ordered{i},'_grayavg'],['R_',regions_ordered{i},'_grayavg']}; 
end

effsizes_edges_68_significant = effsizes_edges_68_significant(regions_ordered_68,regions_ordered_68);

figure;
imagesc(effsizes_edges_68_significant{:,:})
title('Effect sizes edges')
xticks(1:height(atlas))
xticklabels(replace(replace(effsizes_edges_68_significant.Properties.VariableNames,'_',' '),'grayavg',''))
yticks(1:height(atlas))
yticklabels(replace(replace(effsizes_edges_68_significant.Properties.VariableNames,'_',' '),'grayavg',''))
colormap([[0 0 0]; parula])
clim([min(min(effsizes_edges_68)),max(max(effsizes_edges_68))])
colormaps = parula(length(regions));
for i = 1:width(effsizes_edges_68)
    text(i,min(ylim)+105,replace(effsizes_edges_68_significant.Properties.VariableNames(i),'_',' '),'Color',colormaps(atlas_ordered{i,'idx'},:),'Rotation',90,'FontSize', 14)
end
set(gca,'XTickLabel',[]);
ax = gca;
ax.Position(2) = ax.Position(2) + 0.25;
ax.Position(4) = ax.Position(4) - 0.2;

hold on

% Coordenadas para los cuadros de color y las etiquetas
xColorBox = 0.2; % Posición X de los cuadros de color
xLabel = 0.23; % Posición X de las etiquetas
yStart = 0.17; % Posición Y inicial
yStep = 0.027; % Espaciado entre elementos

% Dibujar los cuadros de color y etiquetas
for i = 1:length(regions)
    % Dibujar cuadro de color
    annotation('rectangle', [xColorBox yStart-(i-1)*yStep 0.02 0.02], 'FaceColor', colormaps(i,:), 'EdgeColor', 'none');
    
    % Agregar etiqueta
    annotation('textbox', [xLabel yStart-(i-1)*yStep 0.9 0.03], 'String', regions{i}, 'EdgeColor', 'none', 'FontSize', 11);
end

% Ajustar los ejes y remover
axis off;


figure;
imagesc(effsizes_edges_68_significant{:,:})
title('Effect sizes edges')
colormap([[0 0 0]; parula])
clim([min(min(effsizes_edges_68)),max(max(effsizes_edges_68))])
x_pos_label = [0 histcounts(atlas.idx)];
for i = 1:length(regions)
    text(sum(x_pos_label(1:i)) + x_pos_label(i+1)/2,max(ylim)+3,regions(i),'Color',colormaps(i,:),'Rotation',0,'FontSize', 21,'HorizontalAlignment','center')
    if i > 1
        line([sum(x_pos_label(1:i))+0.5 sum(x_pos_label(1:i))+0.5], ylim, 'Color', 'w', 'LineWidth', 2);
        line([1-0.5,height(effsizes_edges_68_significant)+0.5],[sum(x_pos_label(1:i))+0.5 sum(x_pos_label(1:i))+0.5], ylim, 'Color', 'w', 'LineWidth', 2);
    end
end
set(gca,'XTickLabel',[],'YTickLabel',[]);
