function [var,pval,weights,significant_weights,res,Y_pred] = MIND_10_run_CCA(X,Y,cfg,molecular_names_table,dx_title,X_rot_ordered,type)
    save([cfg.dir.project,'\data\X.mat'],'X')
    save([cfg.dir.project,'\data\Y.mat'],'Y')

    % Update cfg with defaults
    cfg = cfg_defaults(cfg);
    
    % Run analysis
    main(cfg,X_rot_ordered);
    
    % Clean up analysis files to save disc space
    cleanup_files(cfg);

    %----- Visualization

    % Set path for plotting
    set_path('plot');
    
    % Load res
    res.dir.frwork = cfg.dir.frwork;
    res.frwork.level = 1;
    res.gen.selectfile = 'none';
    res = res_defaults(res, 'load');

    res.gen.weight.type = 'correlation'; % for plotting Loadings instead of weights
%     res.behav.weight.norm = 'minmax'; % Normalize weight {'std','zscore'}
    res.gen.weight.subset.type = 'significant'; % Filter weights by keeping only significant ones 
        % {'top' (Filter weights by keeping top ones based on absolute value) 'minmax' (Filter weights by keeping top positive and negative ones)} 
    res.gen.weight.stat.ci = 'sd'; % {'pi'}
    res.behav.weight.errorbar.side = 'one';

    var = load([res.dir.res,'\model',res.env.fileend,'.mat'],'trexvarx'); 
    res_pval = load([res.dir.res,'\res',res.env.fileend,'.mat'],'res');
    pval = res_pval.res.stat.pval;

%     y_weights = load([res.dir.res,'\model',res.env.fileend,'.mat'],'wY'); 
%     y_weights = y_weights.wY;
%     
    
    if length(cfg.machine.param.(cfg.machine.param.name{1})) == 1 && exist([cfg.dir.project,'\VAR_',dx_title,'.mat'])
        cfg_machine_param = cfg.machine.param.(cfg.machine.param.name{1});
    
        for i = 1:length(cfg_machine_param)
            res.frwork.split.best = i; % NGS     
    
            % Plot data projections
%             plot_proj(res, {'X' 'Y'}, res.frwork.level, 'osplit', ...
%                 res.frwork.split.best, 'none', '2d_group', ...
%                 'gen.figure.ext', '.svg', ...
%                 'gen.figure.Position', [0 0 500 400], ...
%                 'gen.axes.Position', [0.1798 0.1560 0.7252 0.7690], ...
%                 'gen.axes.XLim', [-1 1], 'gen.axes.YLim', [-1 1], ...
%                 'gen.axes.FontSize', 22, 'gen.legend.FontSize', 22, ...
%                 'gen.legend.Location', 'best', ...
%                 'proj.scatter.SizeData', 120, ...
%                 'proj.scatter.MarkerFaceColor', [0.3 0.3 0.9;0.9 0.3 0.3], ...
%                 'proj.scatter.MarkerEdgeColor', 'k', 'proj.lsline', 'on', ...
%                 'proj.xlabel', 'Modality 1 latent variable', ...
%                 'proj.ylabel', 'Modality 2 latent variable');
%             title(strrep([dx_title,' ', num2str(cfg_machine_param(i))],'_',' '))

            [~, Y_pred] = plot_proj(res, {'X' 'Y'}, res.frwork.level, 'osplit', ...
                res.frwork.split.best, 'none', '2d',type, 'gen.axes.FontSize', 20, ...
                'gen.legend.FontSize', 20,'gen.legend.Location', 'NorthWest', ...
                'proj.scatter.SizeData', 120, 'proj.scatter.MarkerEdgeColor', 'k', ...
                'proj.xlabel', 'Modality X latent variable', 'proj.ylabel', 'Output variable');
            title(strrep([dx_title,' ', num2str(cfg_machine_param(i))],'_',' '))
% %         

%         % Plot X weights as stem plot
%             plot_weight(res, 'X', 'simul', res.frwork.split.best, 'stem', molecular_names_table, ...
%                 'gen.figure.ext', '.svg', ...
%                 'gen.figure.Position', [0 0 500 400], ...
%                 'gen.axes.Position', [0.1798 0.1560 0.7252 0.7690], ...
%                 'gen.axes.YLim', [-1.1 1.2], ...
%                 'gen.axes.YTick', [-1:0.5:1.2], ...
%                 'gen.axes.FontSize', 22, 'gen.legend.FontSize', 22, ...
%                 'gen.legend.Location', 'NorthEast', ...
%                 'simul.xlabel', 'Modality 1 variables', ...
%                 'simul.ylabel', 'Weight', 'simul.weight.norm', 'minmax');
%             title(strrep([dx_title,' ', num2str(cfg_machine_param(i))],'_',' '))
%             xticks(1:length(X))
%             set(gca, 'FontSize', 11);

            % Plot modality X weights as stem plot
            [weights, significant_weights] = plot_weight(res, 'X', 'behav', res.frwork.split.best, 'behav_vert', type, ...
                'gen.figure.Position', [100 200 1200 500], ...
                'gen.axes.FontSize', 30, 'gen.legend.FontSize', 20, ...
                'behav.weight.errorbar.ci', 'sd', 'behav.weight.errorbar.sd', 1, ...
                'behav.label.order', 'weight', 'behav.weight.display', 'subset-only');

%            weights = plot_weight(res, 'X', 'simul', res.frwork.split.best, 'stem', molecular_names_table, ...
%                 'gen.axes.FontSize',20, 'gen.legend.FontSize', 20, ...
%                 'simul.weight.norm','minmax','simul.xlabel','Modality X variables', ...
%                 'simul.weight.display', 'subset-all', 'simul.weight.subset','abs(stat.Z) >= 2');
%            title(strrep([dx_title,' ', num2str(cfg_machine_param(i))],'_',' '));
%                  
        end   
    end 
    
    if ~exist("weights","var")
        weights = zeros(length(X),1);
    end
    if ~exist("significant_weights","var")
        significant_weights = zeros(length(X),1);
    end
    if ~exist("Y_pred","var")
        Y_pred = zeros(length(Y),1);
    end
end