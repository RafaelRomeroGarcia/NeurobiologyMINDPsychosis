%% Script to find associations between MIND and (1) maturational features 
%% including centiles, volume peak, and velocity peak; 
%% (2) psychosis co-vulnerability; and (3) SA axis

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

% MIND
degree_68_FEP = readtable([location,'degree\degree_68_FEP.csv'],ReadRowNames=true);
degree_68_CN = readtable([location,'degree\degree_68_CN.csv'],ReadRowNames=true);

degree_68_FEP.Properties.VariableNames = replace(replace(degree_68_FEP.Properties.VariableNames,'lh','L'),'rh','R');
degree_68_FEP.Properties.VariableNames = cellfun(@(x) [x, '_grayavg'], degree_68_FEP.Properties.VariableNames, 'UniformOutput', false);
degree_68_CN.Properties.VariableNames = replace(replace(degree_68_CN.Properties.VariableNames,'lh','L'),'rh','R');
degree_68_CN.Properties.VariableNames = cellfun(@(x) [x, '_grayavg'], degree_68_CN.Properties.VariableNames, 'UniformOutput', false);

load("edge_CN.mat")
load("edge_FEP.mat")
load("edge_CN_68.mat")
load("edge_FEP_68.mat")

% EFFECT SIZES
effsizes_degree_68 = readtable([location,'degree\effsizes_degree_68.csv'],ReadRowNames=true);
effsizes_degree_68 = effsizes_degree_68{:,:};
effsizes_degree_68_significant = readtable([location,'degree\effsizes_degree_68_significant.csv'],ReadRowNames=true);
effsizes_degree_68_significant = effsizes_degree_68_significant{:,:};

% SPIN TEST
X_rotated = csvread([location,'perm_sphere_10000_DK.csv'])';
X_rot_ordered = X_rotated(:,:);
X_rot_ordered_34 = X_rotated(1:34,:);

% CENTILES
centiles_old = readtable([location,'centiles\centilesExtendidos_Nat_final.csv'],ReadRowNames=true);
centiles_old.Properties.RowNames = arrayfun(@(p, s) sprintf('%d_%s%02d', p, '0', s), centiles_old.participant, centiles_old.session, 'UniformOutput', false);

total_centiles_68 = readtable([location,'centiles\centiles_CorticalMeasuresENIGMA_GrayAvg.csv'],ReadRowNames=true);
total_centiles_68.Properties.VariableNames = replace(total_centiles_68.Properties.VariableNames,'rh_GM','R');
total_centiles_68.Properties.VariableNames = replace(total_centiles_68.Properties.VariableNames,'lh_GM','L');
index_1_centiles = find(strcmp(total_centiles_68.Properties.VariableNames,'L_bankssts'));
index_2_centiles = find(strcmp(total_centiles_68.Properties.VariableNames,'R_transversetemporal'));
total_centiles_68.Properties.VariableNames(index_1_centiles:end) = cellfun(@(x) [x, '_grayavg'], total_centiles_68.Properties.VariableNames(index_1_centiles:end), 'UniformOutput', false);


% VOLUMES (to check dependency)
volume_68 = readtable([location,'volumes\CorticalMeasuresENIGMA_GrayAvg.csv'],"ReadRowNames",true);

% sort columns alphabetically
volume_68 = volume_68(:,[total_centiles_68.Properties.VariableNames(index_1_centiles:index_2_centiles)]); 

% sort subjects
volume_68 = volume_68(ismember(volume_68.Properties.RowNames,total_centiles_68.Properties.RowNames),:);
volume_68 = volume_68(total_centiles_68.Properties.RowNames,:);

% Session 1
volume_68 = volume_68(total_centiles_68.session==1,:);
centiles_68 = total_centiles_68(total_centiles_68.session==1,:);
volume_68_FEP = volume_68(strcmp(centiles_68.dx,'PS'),:);
volume_68_CN = volume_68(strcmp(centiles_68.dx,'CN'),:);

centiles_68_FEP = centiles_68(strcmp(centiles_68.dx,'PS') & centiles_68.session == 1,:);
centiles_68_CN = centiles_68(strcmp(centiles_68.dx,'CN') & centiles_68.session == 1,:);
centiles_FEP_1_5T_68 = centiles_68(strcmp(centiles_68.study,'PAFIP_1.5T') & strcmp(centiles_68.dx,'PS') & centiles_68.session == 1,:);
centiles_FEP_3T_68 = centiles_68(strcmp(centiles_68.study,'PAFIP_3T') & strcmp(centiles_68.dx,'PS') & centiles_68.session == 1,:);
centiles_CN_1_5T_68 = centiles_68(strcmp(centiles_68.study,'PAFIP_1.5T') & strcmp(centiles_68.dx,'CN') & centiles_68.session == 1,:);
centiles_CN_3T_68 = centiles_68(strcmp(centiles_68.study,'PAFIP_3T') & strcmp(centiles_68.dx,'CN') & centiles_68.session == 1,:);

degree_68_CN = degree_68_CN(ismember(degree_68_CN.Properties.RowNames,centiles_68_CN.Properties.RowNames),:);
degree_68_CN = degree_68_CN(centiles_68_CN.Properties.RowNames,:);
degree_CN_1_5T_68 = degree_68_CN(strcmp(centiles_68_CN.study,'PAFIP_1.5T'),:);
degree_CN_3T_68 = degree_68_CN(strcmp(centiles_68_CN.study,'PAFIP_3T'),:);

% Select common participants between degree and volume and centiles
degree_68_FEP = degree_68_FEP(ismember(degree_68_FEP.Properties.RowNames,centiles_68_FEP.Properties.RowNames),:);
degree_68_FEP = degree_68_FEP(centiles_68_FEP.Properties.RowNames,:);

figure;
scatter(mean(degree_68_FEP{:,:}),mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles}))
lsline
ylabel('Centiles FEP')
xlabel('Degree FEP')
r = corr(mean(degree_68_FEP{:,:})',mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})');
r_par = partialcorr(mean(degree_68_FEP{:,:})',mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',mean(volume_68_FEP{:,:})');
textString = ['r = ', sprintf('%.2f',r),'; ','r_{parcorr} = ', sprintf('%.2f',r_par) ];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.11 0.29 0.08]);

figure;
scatter(mean(degree_68_CN{:,:}),mean(centiles_68_CN{:,index_1_centiles:index_2_centiles}))
lsline
ylabel('Centiles CN')
xlabel('Degree CN')
r = corr(mean(degree_68_CN{:,:})',mean(centiles_68_CN{:,index_1_centiles:index_2_centiles})');
r_par = partialcorr(mean(degree_68_CN{:,:})',mean(centiles_68_CN{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{:,:})');
textString = ['r = ', sprintf('%.2f',r),'; ','r_{parcorr} = ', sprintf('%.2f',r_par) ];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.11 0.29 0.08]);


% Site effect
figure;
scatter(mean(degree_CN_1_5T_68{:,:}),mean(centiles_CN_1_5T_68{:,index_1_centiles:index_2_centiles}))
lsline
ylabel('Centiles CN 1.5T')
xlabel('Degree CN 1.5T')
r = corr(mean(degree_CN_1_5T_68{:,:})',mean(centiles_CN_1_5T_68{:,index_1_centiles:index_2_centiles})');
r_par = partialcorr(mean(degree_CN_1_5T_68{:,:})',mean(centiles_CN_1_5T_68{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{strcmp(centiles_68_CN.study,'PAFIP_1.5T'),:})');
textString = ['r = ', sprintf('%.2f',r),'; ','r_{parcorr} = ', sprintf('%.2f',r_par) ];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.11 0.29 0.08]);


figure;
scatter(mean(degree_CN_3T_68{:,:}),mean(centiles_CN_3T_68{:,index_1_centiles:index_2_centiles}))
lsline
ylabel('Centiles CN 3T')
xlabel('Degree CN 3T')
r = corr(mean(degree_CN_3T_68{:,:})',mean(centiles_CN_3T_68{:,index_1_centiles:index_2_centiles})');
r_par = partialcorr(mean(degree_CN_3T_68{:,:})',mean(centiles_CN_3T_68{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{strcmp(centiles_68_CN.study,'PAFIP_3T'),:})');
textString = ['r = ', sprintf('%.2f',r),'; ','r_{parcorr} = ', sprintf('%.2f',r_par) ];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.11 0.29 0.08]);

figure;
scatter(mean(degree_CN_3T_68{:,:}),mean(degree_CN_1_5T_68{:,:}),'filled')
lsline
ylabel('Degree CN 1.5T')
xlabel('Degree CN 3T')
set(gca, 'FontSize', 14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(mean(degree_CN_1_5T_68{:,:})',mean(degree_CN_3T_68{:,:})');
r_par = partialcorr(mean(degree_CN_3T_68{:,:})',mean(degree_CN_1_5T_68{:,:})',mean(volume_68_CN{:,:})');
nperm = 10000;
for iperm = 1:nperm
    degree_CN_3T_68_perm = degree_CN_3T_68(:,X_rot_ordered(:,iperm));

    r_perm(iperm) = corr(mean(degree_CN_3T_68_perm{:,:})',mean(degree_CN_1_5T_68{:,:})');

    r_par_perm(iperm) = partialcorr(mean(degree_CN_3T_68_perm{:,:})',mean(degree_CN_1_5T_68{:,:})',mean(volume_68_CN{:,:})');
    
end

pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r),'; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_corr),'\color{black}','; ','r_{parcorr} = ', sprintf('%.2f',r_par),'; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_parcorr) ];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.11 0.6 0.08]);


figure;
scatter(mean(centiles_CN_3T_68{:,index_1_centiles:index_2_centiles}),mean(centiles_CN_1_5T_68{:,index_1_centiles:index_2_centiles}))
lsline
ylabel('Centiles CN 1.5T')
xlabel('Centiles CN 3T')
r = corr(mean(centiles_CN_1_5T_68{:,index_1_centiles:index_2_centiles})',mean(centiles_CN_3T_68{:,index_1_centiles:index_2_centiles})');
r_par = partialcorr(mean(centiles_CN_3T_68{:,index_1_centiles:index_2_centiles})',mean(centiles_CN_1_5T_68{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{:,:})');
textString = ['r = ', sprintf('%.2f',r),'; ','r_{parcorr} = ', sprintf('%.2f',r_par) ];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.11 0.29 0.08]);


figure;
scatter(mean(degree_68_CN{:,:}),mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles}),'filled')
lsline
ylabel('Centiles FEP')
xlabel('Degree CN')
set(gca, 'FontSize', 14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',mean(degree_68_CN{:,:})');
r_par = partialcorr(mean(degree_68_CN{:,:})',mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{:,:})');
for iperm = 1:nperm
    centiles_FEP_perm = centiles_68_FEP(:,index_1_centiles - 1 + X_rot_ordered(:,iperm));
    volume_CN_perm = volume_68_CN(:,X_rot_ordered(:,iperm));

    degree_CN_perm = degree_68_CN(:,X_rot_ordered(:,iperm));
    r_perm(iperm) = corr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',mean(degree_CN_perm{:,:})');
%     r_perm(iperm) = corr(mean(centiles_FEP_perm{:,index_1_centiles:index_2_centiles})',mean(degree_CN{:,:})');

    r_par_perm(iperm) = partialcorr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',mean(degree_CN_perm{:,:})',mean(volume_68_CN{:,:})');
    
end

pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' P_{spin} =', sprintf('%.3f',pval_parcorr)];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.15 0.33 0.1]);


figure;
scatter(mean(degree_68_CN{:,:})',mean(volume_68_CN{:,:})')
xlabel('Degree CN')
ylabel('Volume CN')
corr(mean(degree_68_FEP{:,:})',mean(volume_68_CN{:,:})');
hold on
coefficients = polyfit(mean(degree_68_CN{:,:})', mean(volume_68_CN{:,:})', 2);
xFit = linspace(min(mean(degree_68_CN{:,:})'), max(mean(degree_68_CN{:,:})'), 100);
yFit = polyval(coefficients, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);

figure;
scatter(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{:,:})')
ylabel('Volume')
xlabel('Centiles')
corr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{:,:})');
hold on
coefficients = polyfit(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})', mean(volume_68_CN{:,:})', 1);
xFit = linspace(min(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})'), max(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})'), 100);
yFit = polyval(coefficients, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);

% Effsizes degree
figure;
scatter(effsizes_degree_68,mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles}),'filled')
lsline
ylabel('Centiles FEP')
xlabel('Effect sizes degree')
set(gca, 'FontSize', 14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',effsizes_degree_68');
r_par = partialcorr(effsizes_degree_68',mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{:,:})');
for iperm = 1:nperm
    centiles_FEP_perm = centiles_68_FEP(:,index_1_centiles - 1 + X_rot_ordered(:,iperm));
    volume_CN_perm = volume_68_CN(:,X_rot_ordered(:,iperm));
    effsizes_degree_perm = effsizes_degree_68(:,X_rot_ordered(:,iperm));
    
    r_perm(iperm) = corr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',effsizes_degree_perm');
%     r_perm(iperm) = corr(mean(centiles_FEP_perm{:,index_1_centiles:index_2_centiles})',effsizes_degree');

    r_par_perm(iperm) = partialcorr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',effsizes_degree_perm',mean(volume_68_CN{:,:})');
end
pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' p =', sprintf('%.2f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' p =', sprintf('%.3f',pval_parcorr)];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.15 0.29 0.1]);

% figure;
% scatter(effsizes_degree_68,mean(centiles_68_CN{:,index_1_centiles:index_2_centiles}))
% lsline
% ylabel('Centiles CN')
% xlabel('Effect sizes degree')
% r = corr(mean(centiles_68_CN{:,index_1_centiles:index_2_centiles})',effsizes_degree_68');
% r_par = partialcorr(effsizes_degree_68',mean(centiles_68_CN{:,index_1_centiles:index_2_centiles})',mean(volume_68_CN{:,:})');
% for iperm = 1:nperm
%     centiles_CN_perm = centiles_68_CN(:,index_1_centiles - 1 + X_rot_ordered(:,iperm));
%     volume_CN_perm = volume_68_CN(:,X_rot_ordered(:,iperm));
%     effsizes_degree_perm = effsizes_degree_68(:,X_rot_ordered(:,iperm));
%     
%     r_perm(iperm) = corr(mean(centiles_68_CN{:,index_1_centiles:index_2_centiles})',effsizes_degree_perm');
% %     r_perm(iperm) = corr(mean(centiles_CN_perm{:,index_1_centiles:index_2_centiles})',effsizes_degree');
% 
%     r_par_perm(iperm) = partialcorr(mean(centiles_68_CN{:,index_1_centiles:index_2_centiles})',effsizes_degree_perm',mean(volume_68_CN{:,:})');
% end
% pval_corr = sum(abs(r_perm) > abs(r))/nperm;
% pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
% textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' p =', sprintf('%.3f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' p =', sprintf('%.3f',pval_parcorr)];
% annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.15 0.29 0.1]);


% Site effect
figure;
scatter(mean(degree_CN_1_5T_68{:,:}),mean(centiles_FEP_1_5T_68{:,index_1_centiles:index_2_centiles}))
lsline
ylabel('Centiles FEP 1.5T')
xlabel('Degree CN 1.5T')
r = corr(mean(centiles_FEP_1_5T_68{:,index_1_centiles:index_2_centiles})',mean(degree_CN_1_5T_68{:,:})');
r_par = partialcorr(mean(degree_CN_1_5T_68{:,:})',mean(centiles_FEP_1_5T_68{:,index_1_centiles:index_2_centiles})',mean(volume_68{strcmp(centiles_68_CN.study,'PAFIP_1.5T'),:})');
textString = ['r = ', sprintf('%.2f',r),'; ','r_{parcorr} = ', sprintf('%.2f',r_par) ];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.11 0.29 0.08]);


figure;
scatter(mean(degree_CN_3T_68{:,:}),mean(centiles_FEP_3T_68{:,index_1_centiles:index_2_centiles}))
lsline
ylabel('Centiles FEP 3T')
xlabel('Degree CN 3T')
r = corr(mean(centiles_FEP_3T_68{:,index_1_centiles:index_2_centiles})',mean(degree_CN_3T_68{:,:})');
r_par = partialcorr(mean(degree_CN_3T_68{:,:})',mean(centiles_FEP_3T_68{:,index_1_centiles:index_2_centiles})',mean(volume_68{strcmp(centiles_68_CN.study,'PAFIP_3T'),:})');
textString = ['r = ', sprintf('%.2f',r),'; ','r_{parcorr} = ', sprintf('%.2f',r_par) ];
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.17 0.11 0.29 0.08]);


% VOLUMES (to check dependency)
figure;
scatter(mean(volume_68_CN{:,:}),mean(degree_68_CN{:,:}))
xlabel('Volume CN')
ylabel('Degree CN')
lsline
correlation = corr(mean(volume_68_CN{:,:})',mean(degree_68_CN{:,:})');
annotation('textbox', 'String', ['r = ',sprintf('%.2f',correlation)], 'FontSize', 20);


figure;
scatter(mean(volume_68_FEP{:,:}),mean(degree_68_CN{:,:}))
xlabel('Volume FEP')
ylabel('Degree CN')
lsline
correlation = corr(mean(volume_68_FEP{:,:})',mean(degree_68_CN{:,:})');
annotation('textbox', 'String', ['r = ',sprintf('%.2f',correlation)], 'FontSize', 20);


for i = 1:width(volume_68_CN{:,:})
    effsizes_volume_68(:,i) = computeCohen_d(volume_68_FEP{:,i},volume_68_CN{:,i});
end

figure;
scatter(mean(volume_68_CN{:,:}),mean(centiles_68_CN{:,index_1_centiles:index_2_centiles}))
xlabel('Volume CN')
ylabel('Centiles CN')
lsline
correlation = corr(mean(volume_68_CN{:,:})',mean(centiles_68_CN{:,index_1_centiles:index_2_centiles})');
annotation('textbox', 'String', ['r = ',sprintf('%.2f',correlation)], 'FontSize', 20);

figure;
scatter(mean(volume_68_FEP{:,:}),mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles}))
xlabel('Volume FEP')
ylabel('Centiles FEP')
lsline
correlation = corr(mean(volume_68_FEP{:,:})',mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})');
annotation('textbox', 'String', ['r = ',sprintf('%.2f',correlation)], 'FontSize', 20);


%% SUBJECT LEVEL DEGREE
% Volume peak (in years)
peak = readtable([location,'volumes\Table_2_2.csv'],"ReadRowNames",true);
peak = peak(sort(peak.Properties.RowNames),:);

% Average of hemispheres
for i = 1:34
    mean_centiles_CN_68_34(:,i) = mean([centiles_68_CN{:,i + index_1_centiles - 1},centiles_68_CN{:,i + 34 + index_1_centiles - 1}],2);
end
for i = 1:34
    mean_centiles_FEP_68_34(:,i) = mean([centiles_68_FEP{:,i + index_1_centiles - 1},centiles_68_FEP{:,i + 34 + index_1_centiles - 1}],2);
end
for i = 1:34
    mean_degree_CN_68_34(:,i) = mean([degree_68_CN{:,i},degree_68_CN{:,i + 34}],2);
end
for i = 1:34
    mean_degree_FEP_68_34(:,i) = mean([degree_68_FEP{:,i},degree_68_FEP{:,i + 34}],2);
end
for i = 1:34
    mean_volume_CN_68_34(:,i) = mean([volume_68_CN{:,i},volume_68_CN{:,i + 34}],2);
end
for i = 1:34
    mean_volume_FEP_68_34(:,i) = mean([volume_68_FEP{:,i},volume_68_FEP{:,i + 34}],2);
end



for i = 1:width(mean_degree_FEP_68_34)
    effsizes_degree_68_34(:,i) = computeCohen_d(mean_degree_FEP_68_34(:,i),mean_degree_CN_68_34(:,i));
end
for iperm = 1:nperm
    [selection_dx_1,selection_dx_2] = mix_dx(array2table(mean_degree_FEP_68_34),array2table(mean_degree_CN_68_34));
    for ir = 1:34
        effsize_degree_perm_68_34(ir) = computeCohen_d(selection_dx_1(:,ir),selection_dx_2(:,ir));
    end
    effsizes_degree_perm_68_34(iperm,:) = effsize_degree_perm_68_34;
end

pval_eff_degree_68_34 = sum(abs(effsizes_degree_perm_68_34) > abs(effsizes_degree_68_34)) / nperm;
pval_eff_degree_68_34 = mafdr(pval_eff_degree_68_34,'BHFDR',true);
effsizes_degree_significant_68_34 = effsizes_degree_68_34;
effsizes_degree_significant_68_34(pval_eff_degree_68_34 >= 0.05) = NaN;

effsizes_degree_significant_68_34 = array2table(effsizes_degree_significant_68_34,"VariableNames",centiles_old.Properties.VariableNames(index_1_centiles:end));
effsizes_degree_significant_68_34 = effsizes_degree_significant_68_34(:,~isnan(effsizes_degree_significant_68_34{:,:}));

degree_CN_68_34 = mean(mean_degree_CN_68_34);

for i = 1:width(mean_volume_FEP_68_34)
    effsizes_volume_68_34(:,i) = computeCohen_d(mean_volume_FEP_68_34(:,i),mean_volume_CN_68_34(:,i));
end

figure;
scatter(degree_CN_68_34',peak.VolumePeak_inYears_,'filled')
xlabel('Degree CN')
ylabel('Volume peak (in years)')
set(gca, 'FontSize',14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(degree_CN_68_34',peak.VolumePeak_inYears_);
r_par = partialcorr(degree_CN_68_34',peak.VolumePeak_inYears_, mean(mean_volume_CN_68_34)');
for iperm = 1:nperm
    degree_68_34_CN_perm = degree_CN_68_34(:,X_rot_ordered_34(:,iperm));
    r_perm(iperm) = corr(degree_68_34_CN_perm',peak.VolumePeak_inYears_);
    r_par_perm(iperm) = partialcorr(degree_68_34_CN_perm',peak.VolumePeak_inYears_,mean(mean_volume_CN_68_34)');

end
pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' P_{spin} =', sprintf('%.3f',pval_parcorr)];
annotation('textbox','String', textString, 'FontSize', 14,'Position',[0.17 0.15 0.29 0.1]);


figure;
scatter(effsizes_degree_68_34',peak.VolumePeak_inYears_,'filled')
xlabel('Effect sizes of degree')
ylabel('Volume peak (in years)')
set(gca, 'FontSize',14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(effsizes_degree_68_34',peak.VolumePeak_inYears_);
r_par = partialcorr(effsizes_degree_68_34',peak.VolumePeak_inYears_, mean(mean_volume_CN_68_34)');
for iperm = 1:nperm
    effsizes_degree_perm = effsizes_degree_68_34(:,X_rot_ordered_34(:,iperm));
    r_perm(iperm) = corr(effsizes_degree_perm',peak.VolumePeak_inYears_);

%     volume_CN_perm = volume_CN(:,1+X_rot_ordered(:,iperm));
%     peak_perm = peak.VolumePeak_inYears_(X_rot_ordered(:,iperm));
%     r_par_perm(iperm) = partialcorr(effsizes_degree',peak_perm,mean(volume_CN{:,2:end})');
%     r_par_perm(iperm) = partialcorr(effsizes_degree_perm',peak.VolumePeak_inYears_,mean(volume_CN_perm{:,:})');
    r_par_perm(iperm) = partialcorr(effsizes_degree_perm',peak.VolumePeak_inYears_,mean(mean_volume_CN_68_34)');

end

pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' P_{spin} =', sprintf('%.3f',pval_parcorr)];
annotation('textbox','String', textString, 'FontSize', 14,'Position',[0.17 0.15 0.29 0.1]);
lsline

figure;
scatter(mean(mean_centiles_FEP_68_34)',peak.VolumePeak_inYears_,'filled')
xlabel('Centiles FEP')
ylabel('Volume peak (in years)')
set(gca, 'FontSize',14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(mean(mean_centiles_FEP_68_34)',peak.VolumePeak_inYears_);
r_par = partialcorr(peak.VolumePeak_inYears_,mean(mean_centiles_FEP_68_34)',mean(mean_volume_CN_68_34)');
for iperm = 1:nperm
    centiles_FEP_perm = mean_centiles_FEP_68_34(:,X_rot_ordered_34(:,iperm));
    r_perm(iperm) = corr(mean(centiles_FEP_perm)',peak.VolumePeak_inYears_);
    
%     volume_CN_perm = volume_CN(:,1+X_rot_ordered(:,iperm));
%     peak_perm = peak.VolumePeak_inYears_(X_rot_ordered(:,iperm));
%     r_par_perm(iperm) = partialcorr(mean(centiles_FEP_perm{:,:})',peak.VolumePeak_inYears_,mean(volume_CN_perm{:,:})');
%     r_par_perm(iperm) = partialcorr(mean(centiles_FEP{:,:})',peak_perm,mean(volume_CN{:,2:end})');
    r_par_perm(iperm) = partialcorr(mean(centiles_FEP_perm)',peak.VolumePeak_inYears_,mean(mean_volume_CN_68_34)');

end
pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' P_{spin} =', sprintf('%.3f',pval_parcorr)];
annotation('textbox','String', textString, 'FontSize', 14,'Position',[0.17 0.15 0.29 0.1]);
lsline


% Velocity peak (in years)
figure;
scatter(degree_CN_68_34',peak.VelocityPeak_inYears_,'filled')
xlabel('Degree CN')
ylabel('Velocity peak (in years)')
set(gca, 'FontSize',14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(degree_CN_68_34',peak.VelocityPeak_inYears_);
r_par = partialcorr(degree_CN_68_34',peak.VelocityPeak_inYears_, mean(mean_volume_CN_68_34)');
for iperm = 1:nperm
    degree_68_34_CN_perm = degree_CN_68_34(:,X_rot_ordered_34(:,iperm));
    r_perm(iperm) = corr(degree_68_34_CN_perm',peak.VelocityPeak_inYears_);
    r_par_perm(iperm) = partialcorr(degree_68_34_CN_perm',peak.VelocityPeak_inYears_,mean(mean_volume_CN_68_34)');

end
pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' P_{spin} =', sprintf('%.3f',pval_parcorr)];
annotation('textbox','String', textString, 'FontSize', 14,'Position',[0.17 0.15 0.29 0.1]);
lsline


figure;
scatter(effsizes_degree_68_34',peak.VelocityPeak_inYears_,'filled')
xlabel('Effect sizes of degree')
ylabel('Velocity peak (in years)')
set(gca, 'FontSize',14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(effsizes_degree_68_34',peak.VelocityPeak_inYears_);
r_par = partialcorr(effsizes_degree_68_34',peak.VelocityPeak_inYears_, mean(mean_volume_CN_68_34)');
for iperm = 1:nperm
    effsizes_degree_perm = effsizes_degree_68_34(:,X_rot_ordered_34(:,iperm));
    r_perm(iperm) = corr(effsizes_degree_perm',peak.VelocityPeak_inYears_);

    r_par_perm(iperm) = partialcorr(effsizes_degree_perm',peak.VelocityPeak_inYears_,mean(mean_volume_CN_68_34)');
end
pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' P_{spin} =', sprintf('%.3f',pval_parcorr)];
annotation('textbox','String', textString, 'FontSize', 14,'Position',[0.17 0.15 0.29 0.1]);
lsline


figure;
scatter(mean(mean_centiles_FEP_68_34)',peak.VelocityPeak_inYears_,'filled')
xlabel('Centiles FEP')
ylabel('Velocity peak (in years)')
set(gca, 'FontSize',14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
r = corr(mean(mean_centiles_FEP_68_34)',peak.VelocityPeak_inYears_);
r_par = partialcorr(peak.VelocityPeak_inYears_,mean(mean_centiles_FEP_68_34)',mean(mean_volume_CN_68_34)');
for iperm = 1:nperm
    centiles_FEP_perm = mean_centiles_FEP_68_34(:,X_rot_ordered_34(:,iperm));
    r_perm(iperm) = corr(mean(centiles_FEP_perm)',peak.VolumePeak_inYears_);

    r_par_perm(iperm) = partialcorr(mean(centiles_FEP_perm)',peak.VolumePeak_inYears_,mean(mean_volume_CN_68_34)');
end
pval_corr = sum(abs(r_perm) > abs(r))/nperm;
pval_parcorr = sum(abs(r_par_perm) > abs(r_par))/nperm;
textString = ['r = ', sprintf('%.2f',r), '; ', '\color{red}',' P_{spin} =', sprintf('%.3f',pval_corr), '\color{black}','; ', 'r_{parcorr} = ', sprintf('%.2f',r_par),'; ','\color{red}',' P_{spin} =', sprintf('%.3f',pval_parcorr)];
annotation('textbox','String', textString, 'FontSize', 14,'Position',[0.17 0.15 0.29 0.1]);
lsline


%% PSYCHOSIS CO-VULNERABILITY (edges)
load('vulner_effsizes_centiles.mat')

mean_edge_CN = zeros(width(degree_68_FEP));
for i = 1:length(edge_CN_68)
    mean_edge_CN = mean_edge_CN + table2array(edge_CN_68{i,2});
end
mean_edge_CN_68 = mean_edge_CN/length(edge_CN_68);

for i = 1:size(mean_edge_CN_68)
    mean_edge_CN_68(i, i) = NaN;
end

mean_edge_CN_68_34 = (mean_edge_CN_68(1:34,1:34)+mean_edge_CN_68(35:end,35:end))/2;
mean_edge_CN_68_34_col = mean_edge_CN_68_34(find(triu(ones(size(mean_edge_CN_68_34)),1)));

for i = 1:nperm
    mean_edge_CN_perm = mean_edge_CN_68_34(X_rot_ordered(1:34,i),X_rot_ordered(1:34,i));
    mean_edge_CN_perm_col = mean_edge_CN_perm(find(triu(ones(size(mean_edge_CN_perm)),1)));
    r_perm_vulner(i) = corr(mean_edge_CN_perm_col,vulner);
end
correlation_CN_vulner = corr(mean_edge_CN_68_34_col, vulner);
pval_vulner = sum(abs(r_perm_vulner) > abs(correlation_CN_vulner)) / nperm;

figure;
scatter(mean_edge_CN_68_34_col,vulner,'filled')
ylabel('Psychosis co-vulnerability correlation')
xlabel('Edges CN')
set(gca, 'FontSize', 14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1,'HandleVisibility', 'off');
textString = ['r = ', sprintf('%.2f',correlation_CN_vulner), '; ', '\color{red}', ' P_{spin} = ', sprintf('%.2e',pval_vulner)];
annotation('textbox', [0.53 0.21 0.1 0.1], 'String', textString, 'FontSize', 12);


%% Sensorimotor-association-axis ranking
SA_ranking = readtable([location,'\sensorimotor-association_axis_ranking_DK.csv']);
volumes = readtable([location,'\volumes\CorticalMeasuresENIGMA_GrayAvg.csv'],'ReadRowNames',true);
[volumes_names,idx] = sort(volumes(:,1:68).Properties.VariableNames);
SA_ranking = SA_ranking(idx',:);
volume_68_CN = volumes(centiles_68_CN.Properties.RowNames,idx);
volume_68_FEP = volumes(centiles_68_FEP.Properties.RowNames,idx);

r_SA_centiles = corr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',SA_ranking{:,:});
r_par_SA_centiles = partialcorr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',SA_ranking{:,:},mean(volume_68_FEP{:,:})');

r_SA_effsizes_degree = corr(effsizes_degree_68',SA_ranking{:,:});
r_par_SA_effsizes_degree = partialcorr(effsizes_degree_68',SA_ranking{:,:},mean(volume_68_CN{:,:})');


nperm = 10000;
for iperm = 1:nperm
    SA_ranking_perm = SA_ranking(X_rot_ordered(:,iperm),:);
    r_SA_centiles_perm(iperm) = corr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',SA_ranking_perm{:,:});
    r_par_SA_centiles_perm(iperm) = partialcorr(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles})',SA_ranking_perm{:,:},mean(volume_68_FEP{:,:})');

    r_SA_effsizes_degree_perm(iperm) = corr(effsizes_degree_68',SA_ranking_perm{:,:});
    r_par_SA_effsizes_degree_perm(iperm) = partialcorr(effsizes_degree_68',SA_ranking_perm{:,:},mean(volume_68_CN{:,:})');

end

p_SA_centiles = sum(abs(r_SA_centiles_perm) > abs(r_SA_centiles))/nperm;
p_par_SA_centiles = sum(abs(r_par_SA_centiles_perm) > abs(r_par_SA_centiles))/nperm;

p_SA_effsizes_degree = sum(abs(r_SA_effsizes_degree_perm) > abs(r_SA_effsizes_degree))/nperm;
p_par_SA_effsizes_degree = sum(abs(r_par_SA_effsizes_degree_perm) > abs(r_par_SA_effsizes_degree))/nperm;


figure;
scatter(mean(centiles_68_FEP{:,index_1_centiles:index_2_centiles}),SA_ranking{:,:},'filled')
lsline
xlabel('FEP centiles')
ylabel('SA-ranking')
if p_SA_centiles < 0.05
    textString = ['r = ', sprintf('%.2f',r_SA_centiles),'; ','\color{red}','P_{spin} = ',sprintf('%.2f',p_SA_centiles),'; ','\color{black}','r_{parcorr} = ', sprintf('%.2f',r_par_SA_centiles),'; ','\color{red}','P_{spin parcor} = ',sprintf('%.3f',p_par_SA_centiles)];
else
    textString = ['r = ', sprintf('%.2f',r_SA_centiles),'; P_{spin} = ',sprintf('%.2f',p_SA_centiles),'; r_{parcorr} = ', sprintf('%.2f',r_par_SA_centiles),'; P_{spin parcor} = ',sprintf('%.2f',p_par_SA_centiles)];
end    
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.2 0.11 0.65 0.08]);
set(gca, 'FontSize', 14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1,'HandleVisibility', 'off');

figure;
scatter(effsizes_degree_68,SA_ranking{:,:},'filled')
lsline
xlabel('Effect size of MIND degree')
ylabel('SA-ranking')
if p_SA_effsizes_degree < 0.05
    textString = ['r = ', sprintf('%.2f',r_SA_effsizes_degree),'; ','\color{red}','P_{spin} = ',sprintf('%.2f',p_SA_effsizes_degree),'; ','\color{black}','r_{parcorr} = ', sprintf('%.2f',r_par_SA_effsizes_degree),'; ','\color{red}','P_{spin parcor} = ',sprintf('%.3f',p_par_SA_effsizes_degree)];
else
    textString = ['r = ', sprintf('%.2f',r_SA_effsizes_degree),'; P_{spin} = ',sprintf('%.2f',p_SA_effsizes_degree),'; r_{parcorr} = ', sprintf('%.2f',r_par_SA_effsizes_degree),'; P_{spin parcor} = ',sprintf('%.2f',p_par_SA_effsizes_degree)];
end    
annotation('textbox','String', textString, 'FontSize', 10,'Position',[0.2 0.11 0.65 0.08]);
set(gca, 'FontSize', 14)
set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1,'HandleVisibility', 'off');