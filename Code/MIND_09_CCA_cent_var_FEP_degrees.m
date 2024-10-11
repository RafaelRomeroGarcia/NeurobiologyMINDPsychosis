
function [molecular_names, x_weights, significant_weights, Y_pred, correl, molecular_maps_labels_ordered] = MIND_09_CCA_cent_var_FEP_degrees(name,variable_title,data,location,type)
    warning('off','all')
    cd([location,'PCA-CCA\',variable_title,'\'])

    
    X_rotated = csvread([location,'perm_sphere_10000_DK.csv'])';
    X_rot = X_rotated;
    
    response = data{:,:};
   
    
    % Load molecular maps
    molecular_maps = readtable([location,'all_microsc_DesikanKilliany68.csv'],ReadVariableNames=true);
    molecular_names_table = readtable([location,'molecular_names.xlsx'],ReadVariableNames=true);
    molecular_names = molecular_names_table{:,1};
    
    % Important! List of regions may follow a different order in each file. 
    % We need to set the same order for all of them, let´s say alphabetical
    y_labels = data.Properties.VariableNames;
    [~, y_labels_order] = sort(y_labels);
    molecular_maps_labels = strcat(molecular_maps.Var1, '_', molecular_maps.Var2);
    [molecular_maps_labels_ordered, molecular_maps_labels_order] = sort(molecular_maps_labels);
    
    X_rot_ordered = X_rot(molecular_maps_labels_order,:)';
    
    molecular_map_ordered = molecular_maps{molecular_maps_labels_order,3:end};    
    writetable(array2table(molecular_map_ordered,"RowNames",sort(molecular_maps_labels),"VariableNames",molecular_names),[location,'molecular_map_ordered.csv'],"WriteRowNames",true)
    
    
    if length(response) > length(molecular_map_ordered)
        response_ordered = response;
        for j = 1:width(molecular_map_ordered)
            for i = 1:height(molecular_map_ordered)
                for k = 1:height(molecular_map_ordered)
                    molecular_map_edges(i,k) = (molecular_map_ordered(i,j) + molecular_map_ordered(k,j))/2;
                end
            end
            X(:,j) = molecular_map_edges(find(triu(ones(size(molecular_map_edges)),1)));
        end
        
        response_ordered_matrix = zeros(68);
        response_ordered_matrix(find(triu(ones([68,68]),1))) = response;
        response_ordered_matrix = response_ordered_matrix + response_ordered_matrix';
        if ~exist([location,'Data\PCA_CCA\perm_sphere_edges_10000_DK.csv'])
            for iperm = 1:length(X_rot_ordered)
                for i = 1:width(X_rot_ordered)
                    for j = i:width(X_rot_ordered)
                        [idx] = find(response_ordered_matrix == response_ordered_matrix(i,X_rot_ordered(iperm,j)));
                        X_rot_ordered_edges(i,j) = idx(1);
                    end
                end
                X_rot_ordered_edges = X_rot_ordered_edges + X_rot_ordered_edges';
                X_rot_ordered_edges(1:68+1:68^2) = diag(X_rot_ordered_edges)/2;
                X_rot_ordered_edges_col = X_rot_ordered_edges(find(triu(ones(size(X_rot_ordered_edges)),1)));
                X_rot_ordered_edges_vec(iperm,:) = X_rot_ordered_edges_col;
               
            end
        else
            X_rot_ordered = table2array(readtable([location,'Data\PCA_CCA\perm_sphere_edges_10000_DK.csv']));
        end
    else
        response_ordered = response(:,y_labels_order);
        X = molecular_map_ordered;
    end
    Y = mean(response_ordered,1)';
    
    %----- Analysis
    
    % Set path for analysis
    set_path;
    
    % Project folder
    cfg.dir.project = pwd();
    
    % Machine settings
    cfg.machine.name = 'cca';
    
    cfg.machine.metric = {'correl' 'trexvarx' 'trexvary'};
    
    cfg.machine.param.name = {'VARx', 'VARy'}; % explained variance by the PCA components
   
    cfg.machine.param.VARx = 0.6:0.1:0.9; % variance of data kept in the principal components during the SVD step of CCA, PCA-CCA or RCCA
    % % note that if variance is not sufficiently large, only very few (even 0 or 1) variables might be only kept
    cfg.machine.param.VARy = 1;
    
    cfg.machine.svd.varx = 1; % variance of X kept during the SVD step of CCA, PCA-CCA or RCCA. Default is 1 for CCA and 0.99 for RCCA. 
    % Note that if variance is not sufficiently large, only very few (even 0 or 1) variables might be only kept
    cfg.machine.svd.vary = 1; % variance of Y kept during the SVD step of CCA, PCA-CCA or RCCA. Default is 1 for CCA and 0.99 for RCCA. 
    % Note that if variance is not sufficiently large, only very few (even 0 or 1) variables might be only kept, i.e., for models with 1 output variable cfg.machine.svd.vary = 1 should be used
    
    cfg.machine.alignw = 'wX';

    % Framework settings
    cfg.frwork.name = 'permutation'; % In a holdout predictive (machine learning) framework, the data is divided into training and test sets 
    % by randomly subsampling subjects. In a 'permutation' descriptive framework, the data is not splitted, focusing on in-sample statistical evaluation
    
    cfg.frwork.split.nout = 1; % number of outer splits/folds

    cfg.frwork.nlevel = 1;
    
    % Deflation settings
    cfg.defl.name = 'generalized'; % ['generalized', 'pls-projection', 'pls-modeA', 'pls-regression']
    % In the case of iterative methods, once a pair of weights is obtained, the corresponding
    % associative effect is removed from the data (by a process called deflation) and new associations are sought
    
    % Environment settings
    cfg.env.comp = 'local'; %  ['local', 'cluster']
    cfg.env.save.tableHeading = {'set' 'varx' 'correl' 'pval' 'npcax'};
    
    % Number of permutations
    cfg.stat.nperm = 1000;
    cfg.stat.nboot = 1000;


    % Run all VAR
    if ~exist([cfg.dir.project,'\VAR_',variable_title,'.mat'])
        pval_corr = 0.05;
        iter = 0;
        while all(pval_corr >= 0.05)
            if exist([cfg.dir.project,'\framework'])
                rmdir([cfg.dir.project,'\framework'],'s')
            end
            iter = iter + 1;
            cfg.machine.param.VARx = cfg.machine.param.VARx(iter);
            [var_real,pval_spin,~,~,res] = MIND_10_run_CCA(X,Y,cfg,molecular_names_table,variable_title,X_rot_ordered,type);
            pval_spins(iter) = pval_spin;
            cfg.machine.param.VARx = 0.6:0.1:0.9;
            if iter == 1
                pval_corr = pval_spin;
            else
                pval_corr = mafdr(pval_spins,'BHFDR',true); % multipe comparisons correction (Benjamini-Hochberg False Discovery Rate)
            end
            
            if iter == length(cfg.machine.param.VARx)
                break
            end
        end
        
        save(['pval_spins_',variable_title,'.mat'],'pval_spins');
        save(['pval_spins_corr_',variable_title,'.mat'],'pval_corr');
    else
        load(['pval_spins_corr_',variable_title,'.mat'])
        load(['VAR_',variable_title,'.mat'])
    end
   

    if ~exist([cfg.dir.project,'VAR_',variable_title,'.mat'])
        iter = min(find(pval_corr == min(pval_corr)));
        VAR = cfg.machine.param.VARx(iter);
        save(['VAR_',variable_title,'.mat'],'VAR');
        if sum(pval_corr < 0.05) > 1
            rmdir([cfg.dir.project,'\framework'],'s')
        end
%     else
%         load(['VAR_',variable_title,'.mat'])

    end
    
    if any(pval_corr < 0.05)
        cfg.machine.param.VARx = VAR;
    else
        cfg.machine.param.VARx = cfg.machine.param.VARx(1);
    end

    [var_real,pval_spin,weights,significant_weights,res,Y_pred] = MIND_10_run_CCA(X,Y,cfg,molecular_names_table,variable_title,X_rot_ordered,type);
    x_weights = weights';

    if any(pval_corr < 0.05)
        significant_weights = significant_weights';
        Y_pred = Y_pred';
        correl = corr(Y_pred',Y);
    else
        VAR = nan;
        save(['VAR_',variable_title,'.mat'],'VAR');
        significant_weights = zeros(1,length(molecular_names));
        Y_pred = zeros(1,length(Y));
        correl = corr(Y_pred',Y);
    end
    if length(Y_pred) > height(molecular_maps)
        Y_pred_edges = zeros(height(molecular_maps));
        Y_pred_edges(find(triu(ones(size(molecular_map_edges)),1))) = Y_pred;
        Y_pred_edges = Y_pred_edges + Y_pred_edges';
        Y_pred = Y_pred_edges;
    end
    
 end
