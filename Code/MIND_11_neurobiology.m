%% Script to find associations between MIND and neurobiological features 

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
location = 'D:';
location = 'C:\Users\usuario';
location_centiles = [location,'\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\Centiles\'];

molecular_map_ordered = readtable([location,'\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\Data\Molecular_maps\molecular_map_ordered.csv'],"ReadRowNames",true);
molecular_names_ordered = readtable([location,'\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Molecular\Predictive_Maps\molecular_names_ordered_Nat.xlsx'],ReadRowNames=true);

X_rotated = csvread([location,'\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Molecular\parcellations\rotation\perm_sphere_10000_DK.csv'])';
X_rot_ordered = X_rotated(:,:);

load("edge_CN_68.mat")
load("effsizes_edges_68.mat")

% Change as desired
MIND = 'HC';
MIND = 'effect size';

mean_edge_CN = zeros(height(molecular_map_ordered));
for i = 1:length(edge_CN_68)
    mean_edge_CN = mean_edge_CN + table2array(edge_CN_68{i,2});
end
mean_edge_CN_68 = mean_edge_CN/length(edge_CN_68);

for i = 1:size(mean_edge_CN_68)
    mean_edge_CN_68(i, i) = NaN;
end

mean_edge_CN_col = mean_edge_CN_68(find(triu(ones(size(mean_edge_CN_68)),1)));
effsizes_edges_68_col = effsizes_edges_68(find(triu(ones(size(effsizes_edges_68)),1)));
types = unique(molecular_names_ordered.Type,'stable');
figure;
nperm = 10000;
for k = 1:length(types)
    molecular_names_ordered_type_1 = molecular_names_ordered.var_name(molecular_names_ordered.Type_numeric == k);
    molecular_map_ordered_type_1 = molecular_map_ordered(:,molecular_names_ordered_type_1);
    Y_type_1 = corr(molecular_map_ordered_type_1{:,:}');
    Y_type_1 = Y_type_1(find(triu(ones(size(Y_type_1)),1)));
    correlation_CN_neuro_type_1 = corr(mean_edge_CN_col, Y_type_1,'type','Spearman');
    correlation_effsizes_edges_neuro_type_1 = corr(effsizes_edges_68_col, Y_type_1,'type','Spearman');

    % spin test
    for i = 1:nperm
        mean_edge_CN_perm = mean_edge_CN_68(X_rot_ordered(:,i),X_rot_ordered(:,i));
        mean_edge_CN_perm_col = mean_edge_CN_perm(find(triu(ones(size(mean_edge_CN_perm)),1)));
        r_perm_CN_neuro_type_spin(i) = corr(mean_edge_CN_perm_col,Y_type_1,'type','Spearman');

        effsizes_edges_68_perm = effsizes_edges_68(X_rot_ordered(:,i),X_rot_ordered(:,i));
        effsizes_edges_68_perm_col = effsizes_edges_68_perm(find(triu(ones(size(effsizes_edges_68_perm)),1)));
        r_perm_effsizes_edges_neuro_type_spin(i) = corr(effsizes_edges_68_perm_col,Y_type_1,'type','Spearman');
        
    end
    pval_CN_neuro_type_spin = sum(abs(r_perm_CN_neuro_type_spin) > abs(correlation_CN_neuro_type_1)) / nperm;
    pval_effsizes_edges_neuro_type_spin = sum(abs(r_perm_effsizes_edges_neuro_type_spin) > abs(correlation_effsizes_edges_neuro_type_1)) / nperm;    
    

    subplot(2,3,k)
    if strcmp(MIND,'HC')
        [counts, centers] = hist3([mean_edge_CN_col Y_type_1], 'Nbins', [20 20]);
        density = interp2(centers{1}, centers{2}, counts', mean_edge_CN_col,Y_type_1);

        scatter(mean_edge_CN_col,Y_type_1,15,density,'filled')
        xlabel('HC edge')
        if pval_CN_neuro_type_spin < 0.05
            textString = ['r = ', sprintf('%.2f',correlation_CN_neuro_type_1), '; ', '\color{red}', ' P_{spin} = ',sprintf('%.2f',pval_CN_neuro_type_spin)];
        else
            textString = ['r = ', sprintf('%.2f',correlation_CN_neuro_type_1), '; ', ' P_{spin} = ',sprintf('%.2f',pval_CN_neuro_type_spin)];
        end
        set(gca, 'FontSize', 14,'CLim',[0 30])

   
    elseif strcmp(MIND,'effect size')
        [counts, centers] = hist3([effsizes_edges_68_col Y_type_1], 'Nbins', [20 20]);
        density = interp2(centers{1}, centers{2}, counts', effsizes_edges_68_col,Y_type_1);
    
        scatter(effsizes_edges_68_col,Y_type_1,15,density,'filled')
        xlabel('Effect size of edge')
        if pval_effsizes_edges_neuro_type_spin < 0.05
            textString = ['r = ', sprintf('%.2f',correlation_effsizes_edges_neuro_type_1), '; ', '\color{red}', ' P_{spin} = ',sprintf('%.2f',pval_effsizes_edges_neuro_type_spin)];
        else
            textString = ['r = ', sprintf('%.2f',correlation_effsizes_edges_neuro_type_1), '; ', ' P_{spin} = ',sprintf('%.2f',pval_effsizes_edges_neuro_type_spin)];

        end
        set(gca, 'FontSize', 14,'CLim',[0 25])

    end
    ylabel(sprintf('Neurobiological similarity'))
    title(types{k})
    set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1,'HandleVisibility', 'off');
    colorbar

    pos = get(gca, 'Position'); 
    x_annot = pos(1) + pos(3)/1.7; 
    y_annot = pos(2);   
    annotation('textbox', [x_annot y_annot 0.1 0.1], 'String', textString, 'FontSize', 12);

   
    set(gca, 'FontSize', 14)
    set(lsline, 'Color', [0 0.4470 0.7410], 'LineWidth', 1,'HandleVisibility', 'off');
    colorbar
end
