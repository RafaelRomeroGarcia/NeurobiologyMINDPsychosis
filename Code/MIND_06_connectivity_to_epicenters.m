%% Script to compute regions connected to disease epicenters

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

location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\2. MIND\Data\';
mkdir([location,'connectivity_to_epicenters'])
directory = [location,'connectivity_to_epicenters']; % change as desired

cd(directory)


% Load dataset
load("degree_68_CN.mat")
load("degree_68_FEP.mat")
load("edge_CN_68.mat")
load("edge_FEP_68.mat")

centiles = readtable([location,'centiles\centiles_CorticalMeasuresENIGMA_GrayAvg.csv'],"ReadRowNames",true);
centiles = centiles(centiles.session == 1,:);
centiles_CN = centiles(strcmp(centiles.dx,'CN'),:);
centiles_FEP = centiles(strcmp(centiles.dx,'PS'),:);

index_1_centiles = find(strcmp(centiles.Properties.VariableNames,'lh_GM_bankssts'));
index_2_centiles = find(strcmp(centiles.Properties.VariableNames,'rh_GM_transversetemporal'));

degree_CN_68 = degree_68_CN(centiles_CN.Properties.RowNames,:); 
degree_FEP_68 = degree_68_FEP(centiles_FEP.Properties.RowNames,:); 

edge_CN_68 = edge_CN_68(ismember(edge_CN_68(:,1),centiles_CN.Properties.RowNames),:);

centiles_cognition = readtable([location,'centiles\ClinicalAndCentilesInformationDataBase_11_12_2023.csv'],"ReadRowNames",true);
centiles_cognition_FEP_and_CN = centiles_cognition(centiles_cognition.Assessment == 1 ,:);
centiles_cognition_CN = centiles_cognition_FEP_and_CN(centiles_CN.Properties.RowNames(ismember(centiles_CN.Properties.RowNames,centiles_cognition_FEP_and_CN.Properties.RowNames)),:);
centiles_cognition = centiles_cognition(centiles_cognition.Assessment == 1 & ismember(centiles_cognition.Subject,centiles_FEP.participant),:);
centiles_cognition = [centiles_cognition(:,{'Subject','Global_Cognitive_Functioning_average'}), array2table(str2double(centiles_cognition{:,{'BPRS_Total','SANS_Total','SAPS_Total','CPZ_equivalent'}}),'VariableNames',{'BPRS_Total','SANS_Total','SAPS_Total','CPZ_equivalent'})];
centiles_cognition_FEP_and_CN = [centiles_cognition_FEP_and_CN(:,{'Subject','Global_Cognitive_Functioning_average'}), array2table(str2double(centiles_cognition_FEP_and_CN{:,{'BPRS_Total','SANS_Total','SAPS_Total','CPZ_equivalent'}}),'VariableNames',{'BPRS_Total','SANS_Total','SAPS_Total','CPZ_equivalent'})];


% Centiles vs MIND CN
edges_CN_averaged = zeros(68);
nperm=10000;
for i = 1:length(edge_CN_68)
    for j = 1:68
        subject_region_MIND = edge_CN_68{i,2}{j,:}; 
        [epicenters_centiles_MIND_subjects_CN(i,j),p_epicenters_centiles_MIND_subjects_CN(i,j)] = corr(subject_region_MIND',mean(centiles_CN{:,index_1_centiles:index_2_centiles})');
    end
    edges_CN_averaged = edges_CN_averaged + edge_CN_68{i,2}{:,:};
end

edges_CN_averaged = edges_CN_averaged/length(edge_CN_68);

epicenters_centiles_MIND_CN = mean(epicenters_centiles_MIND_subjects_CN);
writetable(array2table(epicenters_centiles_MIND_CN,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_CN.csv')

p_epicenters_centiles_MIND_CN = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_CN .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_CN.^2))), 68-1));
p_epicenters_centiles_MIND_CN = mafdr(p_epicenters_centiles_MIND_CN,'BHFDR','true');
epicenters_centiles_MIND_CN_significant = epicenters_centiles_MIND_CN;
epicenters_centiles_MIND_CN_significant(p_epicenters_centiles_MIND_CN>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_CN_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_CN.csv')


figure;
imagesc(epicenters_centiles_MIND_subjects_CN)
title('Epicenters centiles MIND subjects CN')


% Centiles vs MIND FEP
edges_FEP_averaged = zeros(68);
for i = 1:length(edge_FEP_68)
    for j = 1:68
        subject_region_MIND_FEP = edge_FEP_68{i,2}{j,:}; 
        [epicenters_centiles_MIND_subjects_FEP_FEP(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{:,index_1_centiles:index_2_centiles})');
   
        [epicenters_centiles_MIND_subjects_FEP_FEP_low_cognition(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP_low_cognition(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{centiles_cognition{:,"Global_Cognitive_Functioning_average"} < median(centiles_cognition{:,"Global_Cognitive_Functioning_average"},'omitnan'),index_1_centiles:index_2_centiles})');       
        [epicenters_centiles_MIND_subjects_FEP_FEP_high_cognition(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP_high_cognition(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{centiles_cognition{:,"Global_Cognitive_Functioning_average"} >= median(centiles_cognition{:,"Global_Cognitive_Functioning_average"},'omitnan'),index_1_centiles:index_2_centiles})');          
        [epicenters_centiles_MIND_subjects_FEP_FEP_low_BPRS(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP_low_BPRS(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{centiles_cognition{:,"BPRS_Total"} < median(centiles_cognition{:,"BPRS_Total"},'omitnan'),index_1_centiles:index_2_centiles})');            
        [epicenters_centiles_MIND_subjects_FEP_FEP_high_BPRS(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP_high_BPRS(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{centiles_cognition{:,"BPRS_Total"} >= median(centiles_cognition{:,"BPRS_Total"},'omitnan'),index_1_centiles:index_2_centiles})');            

        [epicenters_centiles_MIND_subjects_FEP_FEP_low_SAPS(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP_low_SAPS(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{centiles_cognition{:,"SAPS_Total"} < median(centiles_cognition{:,"SAPS_Total"},'omitnan'),index_1_centiles:index_2_centiles})');       
        [epicenters_centiles_MIND_subjects_FEP_FEP_high_SAPS(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP_high_SAPS(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{centiles_cognition{:,"SAPS_Total"} >= median(centiles_cognition{:,"SAPS_Total"},'omitnan'),index_1_centiles:index_2_centiles})');          
        [epicenters_centiles_MIND_subjects_FEP_FEP_low_SANS(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP_low_SANS(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{centiles_cognition{:,"SANS_Total"} < median(centiles_cognition{:,"SANS_Total"},'omitnan'),index_1_centiles:index_2_centiles})');            
        [epicenters_centiles_MIND_subjects_FEP_FEP_high_SANS(i,j),p_epicenters_centiles_MIND_subjects_FEP_FEP_high_SANS(i,j)] = corr(subject_region_MIND_FEP',mean(centiles_FEP{centiles_cognition{:,"SANS_Total"} >= median(centiles_cognition{:,"SANS_Total"},'omitnan'),index_1_centiles:index_2_centiles})');            
    end

    edges_FEP_averaged = edges_FEP_averaged + edge_FEP_68{i,2}{:,:};
end

edges_FEP_averaged = edges_FEP_averaged/length(edge_FEP_68);


epicenters_centiles_MIND_FEP = mean(epicenters_centiles_MIND_subjects_FEP_FEP);
writetable(array2table(epicenters_centiles_MIND_FEP,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP.csv')

epicenters_centiles_MIND_FEP_FEP_low_cognition = mean(epicenters_centiles_MIND_subjects_FEP_FEP_low_cognition);
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_low_cognition,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP_low_cognition.csv')
epicenters_centiles_MIND_FEP_FEP_high_cognition = mean(epicenters_centiles_MIND_subjects_FEP_FEP_high_cognition);
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_high_cognition,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP_high_cognition.csv')
epicenters_centiles_MIND_FEP_FEP_low_BPRS = mean(epicenters_centiles_MIND_subjects_FEP_FEP_low_BPRS);
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_low_BPRS,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP_low_BPRS.csv')
epicenters_centiles_MIND_FEP_FEP_high_BPRS = mean(epicenters_centiles_MIND_subjects_FEP_FEP_high_BPRS);
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_high_BPRS,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP_high_BPRS.csv')

epicenters_centiles_MIND_FEP_FEP_low_SAPS = mean(epicenters_centiles_MIND_subjects_FEP_FEP_low_SAPS);
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_low_SAPS,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP_low_SAPS.csv')
epicenters_centiles_MIND_FEP_FEP_high_SAPS = mean(epicenters_centiles_MIND_subjects_FEP_FEP_high_SAPS);
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_high_SAPS,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP_high_SAPS.csv')
epicenters_centiles_MIND_FEP_FEP_low_SANS = mean(epicenters_centiles_MIND_subjects_FEP_FEP_low_SANS);
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_low_SANS,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP_low_SANS.csv')
epicenters_centiles_MIND_FEP_FEP_high_SANS = mean(epicenters_centiles_MIND_subjects_FEP_FEP_high_SANS);
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_high_SANS,"VariableNames",degree_CN_68.Properties.VariableNames),'epicenters_centiles_MIND_FEP_FEP_high_SANS.csv')

p_epicenters_centiles_MIND_FEP_FEP = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP = mafdr(p_epicenters_centiles_MIND_FEP_FEP,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_significant = epicenters_centiles_MIND_FEP;
epicenters_centiles_MIND_FEP_FEP_significant(p_epicenters_centiles_MIND_FEP_FEP>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP.csv')

p_epicenters_centiles_MIND_FEP_FEP_low_cognition = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP_FEP_low_cognition .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP_FEP_low_cognition.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP_low_cognition = mafdr(p_epicenters_centiles_MIND_FEP_FEP_low_cognition,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_low_cognition_significant = epicenters_centiles_MIND_FEP_FEP_low_cognition;
epicenters_centiles_MIND_FEP_FEP_low_cognition_significant(p_epicenters_centiles_MIND_FEP_FEP_low_cognition>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_low_cognition_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP_low_cognition.csv')

p_epicenters_centiles_MIND_FEP_FEP_high_cognition = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP_FEP_high_cognition .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP_FEP_high_cognition.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP_high_cognition = mafdr(p_epicenters_centiles_MIND_FEP_FEP_high_cognition,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_high_cognition_significant = epicenters_centiles_MIND_FEP_FEP_high_cognition;
epicenters_centiles_MIND_FEP_FEP_high_cognition_significant(p_epicenters_centiles_MIND_FEP_FEP_high_cognition>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_high_cognition_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP_high_cognition.csv')

p_epicenters_centiles_MIND_FEP_FEP_low_BPRS = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP_FEP_low_BPRS .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP_FEP_low_BPRS.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP_low_BPRS = mafdr(p_epicenters_centiles_MIND_FEP_FEP_low_BPRS,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_low_BPRS_significant = epicenters_centiles_MIND_FEP_FEP_low_BPRS;
epicenters_centiles_MIND_FEP_FEP_low_BPRS_significant(p_epicenters_centiles_MIND_FEP_FEP_low_BPRS>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_low_BPRS_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP_low_BPRS.csv')

p_epicenters_centiles_MIND_FEP_FEP_high_BPRS = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP_FEP_high_BPRS .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP_FEP_high_BPRS.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP_high_BPRS = mafdr(p_epicenters_centiles_MIND_FEP_FEP_high_BPRS,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_high_BPRS_significant = epicenters_centiles_MIND_FEP_FEP_high_BPRS;
epicenters_centiles_MIND_FEP_FEP_high_BPRS_significant(p_epicenters_centiles_MIND_FEP_FEP_high_BPRS>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_high_BPRS_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP_high_BPRS.csv')

p_epicenters_centiles_MIND_FEP_FEP_low_SAPS = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP_FEP_low_SAPS .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP_FEP_low_SAPS.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP_low_SAPS = mafdr(p_epicenters_centiles_MIND_FEP_FEP_low_SAPS,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_low_SAPS_significant = epicenters_centiles_MIND_FEP_FEP_low_SAPS;
epicenters_centiles_MIND_FEP_FEP_low_SAPS_significant(p_epicenters_centiles_MIND_FEP_FEP_low_SAPS>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_low_SAPS_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP_low_SAPS.csv')

p_epicenters_centiles_MIND_FEP_FEP_high_SAPS = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP_FEP_high_SAPS .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP_FEP_high_SAPS.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP_high_SAPS = mafdr(p_epicenters_centiles_MIND_FEP_FEP_high_SAPS,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_high_SAPS_significant = epicenters_centiles_MIND_FEP_FEP_high_SAPS;
epicenters_centiles_MIND_FEP_FEP_high_SAPS_significant(p_epicenters_centiles_MIND_FEP_FEP_high_SAPS>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_high_SAPS_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP_high_SAPS.csv')

p_epicenters_centiles_MIND_FEP_FEP_low_SANS = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP_FEP_low_SANS .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP_FEP_low_SANS.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP_low_SANS = mafdr(p_epicenters_centiles_MIND_FEP_FEP_low_SANS,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_low_SANS_significant = epicenters_centiles_MIND_FEP_FEP_low_SANS;
epicenters_centiles_MIND_FEP_FEP_low_SANS_significant(p_epicenters_centiles_MIND_FEP_FEP_low_SANS>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_low_SANS_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP_low_SANS.csv')

p_epicenters_centiles_MIND_FEP_FEP_high_SANS = 2 * (1 - tcdf(abs(epicenters_centiles_MIND_FEP_FEP_high_SANS .* sqrt((68-1 - 2) ./ (1 - epicenters_centiles_MIND_FEP_FEP_high_SANS.^2))), 68-1));
p_epicenters_centiles_MIND_FEP_FEP_high_SANS = mafdr(p_epicenters_centiles_MIND_FEP_FEP_high_SANS,'BHFDR','true');
epicenters_centiles_MIND_FEP_FEP_high_SANS_significant = epicenters_centiles_MIND_FEP_FEP_high_SANS;
epicenters_centiles_MIND_FEP_FEP_high_SANS_significant(p_epicenters_centiles_MIND_FEP_FEP_high_SANS>0.05) = nan;
writetable(array2table(epicenters_centiles_MIND_FEP_FEP_high_SANS_significant,"VariableNames",degree_FEP_68.Properties.VariableNames),'significant_epicenters_centiles_MIND_FEP_FEP_high_SANS.csv')