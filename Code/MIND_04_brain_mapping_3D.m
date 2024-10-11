
%% Script to plot a Brain Net from MIND effect sizes

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
cd(location)
mkdir(location,'edges\plots')

effsizes_degree = readtable('degree\effsizes_degree_68.csv','ReadRowNames',true);
effsizes_degree_significant = readtable('degree\effsizes_degree_68_significant.csv','ReadRowNames',true);
node_file = readtable([location,'Desikan-Killiany68_Nat.txt']);

load('edges\effsizes_edges_68.mat')
load("edges\pval_effsizes_edges_corrected.mat")

symptoms = {'','_cognition','_BPRS','_SAPS','_SANS'};
threshold = 0.05;
for i = 1:length(symptoms)
    if i == 1
        high_lows = {''};
    else
        high_lows = {'_low','_high'};

    end

    for j = 1:length(high_lows)
       
        % Nodes (degrees)
        effsizes_degree = readtable(['degree\effsizes_degree_68',high_lows{j},symptoms{i},'.csv'],'ReadRowNames',true);
        if strcmp(high_lows,'_low_high_88')
            effsizes_degree_significant = readtable(['degree\effsizes_degree_68',high_lows{j},symptoms{i},'_ttest_significant.csv'],'ReadRowNames',true);
        else
            effsizes_degree_significant = readtable(['degree\effsizes_degree_68',high_lows{j},symptoms{i},'_significant.csv'],'ReadRowNames',true);
        end
        node_file{:,'Var4'} = effsizes_degree{:,:}';
        node_file{:,'Var5'} = ~isnan(effsizes_degree_significant{:,:})';
        save(['degree\node',high_lows{j},symptoms{i},'.node'], 'node_file', '-ascii');

        fileID = fopen(['degree\node',high_lows{j},symptoms{i},'.node'], 'w');

        % Escribir encabezado
        fprintf(fileID, '#Desikan-Killiany Atlas 68 nodes\n');  % Ajusta el encabezado según tu estructura de datos
        
        % Escribir datos
        for k = 1:size(node_file, 1)
            fprintf(fileID, '%d %d %d %d %d %s\n', node_file{k,"Var1"}, node_file{k,"Var2"},node_file{k,"Var3"},node_file{k,"Var4"},node_file{k,"Var5"},node_file{k,"Var6"}{:});
        end
        
        fclose(fileID); % Cerrar el archivo

        % Edges
        load(['edges\pval_effsizes_edges_corrected',high_lows{j},symptoms{i},'.mat'])
        if strcmp(high_lows{j},'_low') || strcmp(high_lows{j},'_high') || strcmp(high_lows{j},'')
            pval_effsizes_edges_68 = eval(['pval_effsizes_edges_corrected',high_lows{j}]);
            connectivity_threshold = double(pval_effsizes_edges_68 < threshold);
            load(['edges\effsizes_edges_68',high_lows{j},symptoms{i},'.mat'])
            effsizes_edges_68 = eval(['effsizes_edges_68',high_lows{j},symptoms{i}]);
        else
            pval_effsizes_edges_68 = eval('pval_effsizes_edges_corrected_low_high');
            connectivity_threshold = double(pval_effsizes_edges_68 < threshold);
            load(['edges\effsizes_edges_68',high_lows{j},symptoms{i},'.mat'])
            effsizes_edges_68 = eval(['effsizes_edges_68_low_high',symptoms{i}]);            
        end
        
        for l = 1:size(effsizes_edges_68)
            effsizes_edges_68(l, l) = 0;
        end
        effsizes_edges = effsizes_edges_68.*connectivity_threshold;

        save(['edges\effsizes_edges',high_lows{j},symptoms{i},'.edge'], 'effsizes_edges', '-ascii');

        lim_nodes(2*i-mod(j,2),:) = [min(node_file{:,"Var4"}) max(node_file{:,"Var4"})];
        lim_edges(2*i-mod(j,2),:) = [min(min(effsizes_edges)) max(max(effsizes_edges))];
        
      % Parameters (modify BrainNet_MapCfg lines 140-142
        EC.msh.alpha = 0.14;  % surface opacity
        
        EC.lot.view_direction = 1;  % 1 sagital view 
                                    % 2 axial view
                                    % 3 coronal view
    
        EC.lot.view = 1;  % 1 single view
                        % 2 full view
                        % 3 lateral and medial view
                        % 4 lateral, medial and ventral view
                        % 5 lateral, medial and dorsal view

      % Initialize
%         H_BrainNet = BrainNet_MapCfg('BrainMesh_ICBM152.nv','Desikan-Killiany68.node',['connectivity_',high_lows{j},'_',symptoms{i},'.edge']);

        H_BrainNet = BrainNet_MapCfg('BrainMesh_ICBM152.nv',['degree\node',high_lows{j},symptoms{i},'.node'],['edges\effsizes_edges',high_lows{j},symptoms{i},'.edge'],['edges\plots\effsizes',high_lows{j},symptoms{i},'.png']);

    end 
end

