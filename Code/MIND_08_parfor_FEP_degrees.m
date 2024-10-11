%% Script to run in parallel multiple PCA-CCA of FEP degrees

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

% Change as desired
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\2. MIND\Data\';
mkdir(location,'PCA-CCA')
directory = [location,'PCA-CCA\'];

 % Load dataset (Comment as desired)
 
% data = readtable([location, 'degree\effsizes_degree_68.csv'],"ReadRowNames",true);
% variable_title = {'effsizes_degrees'};
%  
data = readtable([location, 'degree\degree_68_CN.csv'],"ReadRowNames",true);
variable_title = {'CN_degrees'};

% symptoms = {'cognition','BPRS','SAPS','SANS'};
% high_lows = {'low_high_88'};
% k=1;
% for high_low = high_lows
%     for symptom = symptoms
%         data_symptom = readtable([location, 'degree\effsizes_degree_68_',high_low{:},'_',symptom{:},'.csv'],"ReadRowNames",true);
%         data(k,:) = data_symptom;
%         variable_title{k,:} = ['effsizes_degrees_',high_low{:},'_',symptom{:}];
%         k = k+1;
%     end 
% end

name = {};

for i = 1:size(variable_title,1)
    % Create subfolders
    if ~exist(fullfile(directory,variable_title{i}),'dir')
       mkdir([directory,variable_title{i}])
    end
    if ~exist(fullfile([directory,variable_title{i}],'data'),'dir')
        mkdir([directory,variable_title{i}],'data')
    end
    if exist(fullfile([directory,variable_title{i}],'Copy_of_framework'),'dir')
        if exist(fullfile([directory,variable_title{i}],'framework'),'dir')
            rmdir(fullfile([directory,variable_title{i}],'framework'),'s')
        end
        movefile(fullfile([directory,variable_title{i}],'Copy_of_framework\'),fullfile([directory,variable_title{i}],'framework\'))
        pause(5)
    end
    
    [molecular_names, weights_real, significant_weights, Y_pred, correl_real, molecular_maps_labels_ordered] = MIND_09_CCA_cent_var_FEP_degrees(name,variable_title{i},data(i,:),location,'real');
    weights_total = weights_real;
    significant_weights_total = significant_weights;
    Y_pred_total = Y_pred;
    save('correl_real.mat',"correl_real")
    save('weights_real.mat',"weights_real")
    
    
    weights_total = array2table(weights_total,"RowNames","FEP","VariableNames",molecular_names);
    significant_weights_total = array2table(significant_weights_total,"RowNames","FEP","VariableNames",molecular_names);
    i_significative = any(weights_total{:,:}~=0,2); % same for weight_total
    weights_degrees = weights_total(i_significative,:);
    significant_weights_degrees = significant_weights_total(i_significative,:);
    
    if ~all(Y_pred==0)
        if size(Y_pred,1) == 1
            Y_pred_table_degrees = array2table(Y_pred_total(any(Y_pred_total~=0,2),:),'RowNames',"FEP",'VariableNames',molecular_maps_labels_ordered);
        else
            Y_pred_table_degrees = array2table(Y_pred_total(any(Y_pred_total~=0,2),:),'RowNames',molecular_maps_labels_ordered,'VariableNames',molecular_maps_labels_ordered);
        end
    end
    molecular_names_table = readtable([location,'molecular_names.xlsx'],ReadVariableNames=true,ReadRowNames=true);
    molecular_types_orig = molecular_names_table{:,2};
    weights_degrees =  rows2vars(weights_degrees);
    weights_degrees.molecular_types(1:46) = molecular_types_orig;
    weights_degrees = sortrows(weights_degrees,{'molecular_types','OriginalVariableNames'});
    weights_degrees =  rows2vars(weights_degrees(:,1:end-1),"VariableNamesSource","OriginalVariableNames");
    weights_degrees.Properties.RowNames = weights_degrees{:,1};
    weights_degrees = weights_degrees(:,2:end);
    
    weights_degrees.Properties.VariableNames(7) = {'α_4β_2'};
    weights_degrees.Properties.VariableNames(1) = {'5-HT_{1A}'};
    weights_degrees.Properties.VariableNames(2) = {'5-HT_{1B}'};
    weights_degrees.Properties.VariableNames(3) = {'5-HT_{2A}'};
    weights_degrees.Properties.VariableNames(4) = {'5-HT_4'};
    weights_degrees.Properties.VariableNames(5) = {'5-HT_6'};
    weights_degrees.Properties.VariableNames(6) = {'5-HTT'};
    weights_degrees.Properties.VariableNames(8) = insertBefore(weights_degrees.Properties.VariableNames(8), 3, '_');
    weights_degrees.Properties.VariableNames([9,10,13,14]) = insertBefore(weights_degrees.Properties.VariableNames([9,10,13,14]), 2, '_');
    weights_degrees.Properties.VariableNames(12) = {'GABA'};
    weights_degrees.Properties.VariableNames(19) = insertBefore(weights_degrees.Properties.VariableNames(19), 6, '_');
    weights_degrees.Properties.VariableNames([23,24]) = replace(weights_degrees.Properties.VariableNames([23,24]),'_',' ');
    weights_degrees.Properties.VariableNames(27) = {'Layer I'};
    weights_degrees.Properties.VariableNames(28) = {'Layer II'};
    weights_degrees.Properties.VariableNames(29) = {'Layer III'};
    weights_degrees.Properties.VariableNames(30) = {'Layer IV'};
    weights_degrees.Properties.VariableNames(31) = {'Layer V'};
    weights_degrees.Properties.VariableNames(32) = {'Layer VI'};
    weights_degrees.Properties.VariableNames(33) = {'Gene PC1'};
    weights_degrees.Properties.VariableNames(34) = {'Myelin'};
    weights_degrees.Properties.VariableNames(35) = {'Neurotransmitter PC1'};
    weights_degrees.Properties.VariableNames(36) = {'Synapse density'};
    weights_degrees.Properties.VariableNames(37) = {'Thickness'};
    weights_degrees.Properties.VariableNames(38) = {'Developmental exp.'};
    weights_degrees.Properties.VariableNames(39) = {'Evolutionary exp.'};
    weights_degrees.Properties.VariableNames(40) = {'Scaling NIH'};
    weights_degrees.Properties.VariableNames(41) = {'Scaling PNC'};
    weights_degrees.Properties.VariableNames(42) = {'CBF'};
    weights_degrees.Properties.VariableNames(43) = {'CBV'};
    weights_degrees.Properties.VariableNames(44) = {'CMRO_2'};
    weights_degrees.Properties.VariableNames(45) = {'CMRGlu'};
    weights_degrees.Properties.VariableNames(46) = {'Glycolytic index'};
    
    if strcmp(gca().YLabel.String,'Weight')
        save([directory,'weights_degrees.mat'],"weights_degrees")
    elseif strcmp(gca().YLabel.String,'Loading')
        loadings_degrees = weights_degrees;
        save([directory,'loadings_degrees.mat'],"loadings_degrees")
    end
    
    writetable(Y_pred_table_degrees,[directory,'Y_pred_',variable_title{i},'.csv'],'WriteRowNames',true)
end
ax = gca;
ax.Position(1) = ax.Position(1) - 0.05; %izda
ax.Position(2) = ax.Position(2) - 0.2; %arriba
ax.Position(3) = ax.Position(3) + 0.15; %dcha
ax.Position(4) = ax.Position(4) + 0.2; %abajo
