function [selection_dx_1,selection_dx_2,ratio_1,ratio_2] = mix_dx(centiles_dx_1,centiles_dx_2)
    
    % Script to create groups by randomly mixing patients with different dx from different studies (ALSPAC, PAFIP).

    if height(centiles_dx_1) <= height(centiles_dx_2)
    else
        temp = centiles_dx_1;
        centiles_dx_1 = centiles_dx_2;
        centiles_dx_2 = temp;
    end

    centiles_dx_1 = table2array(centiles_dx_1(:,sort(centiles_dx_1.Properties.VariableNames))); % IMPORTANT
    centiles_dx_2 = table2array(centiles_dx_2(:,sort(centiles_dx_2.Properties.VariableNames)));

    ratio = size(centiles_dx_1,1)/(size(centiles_dx_1,1)+size(centiles_dx_2,1));


%     ratio_range = ratio-ratio*0.1:ratio+ratio*0.1;
%     ratio_range((ratio_range > 1)) = 0.9;
%     ratio_range((ratio_range < 0)) = 0.1;

    random_ratio = ratio*0.9 + rand(1,1) * (ratio*0.2);  %ratio_range(randi(length(ratio_range)));
    random_ratio = max([0 random_ratio]);
    random_ratio = min([1 random_ratio]);

    num_filas_matriz_1 = round(random_ratio*size(centiles_dx_1,1));
    num_filas_matriz_2 = size(centiles_dx_1,1) - num_filas_matriz_1;
  
%     dx_1_index_1 = randsample(size(centiles_dx_1,1),randi([1,size(centiles_dx_1,1)]));
    dx_1_index_1 = randsample(size(centiles_dx_1,1),num_filas_matriz_1);
    dx_1_index_2 = setdiff((1:size(centiles_dx_1,1))',dx_1_index_1);

%     dx_2_index_1 = randsample(size(centiles_dx_2,1),size(dx_1_index_2,1));
    dx_2_index_1 = randsample(size(centiles_dx_2,1),num_filas_matriz_2);
    dx_2_index_2 = setdiff((1:size(centiles_dx_2,1))',dx_2_index_1);
            
    selection_dx_1 = [centiles_dx_1(dx_1_index_1,:);centiles_dx_2(dx_2_index_1,:)];
    selection_dx_2 = [centiles_dx_1(dx_1_index_2,:);centiles_dx_2(dx_2_index_2,:)];

    ratio_1 = numel(dx_1_index_1)/size(selection_dx_1,1);
    ratio_2 = numel(dx_1_index_2)/size(selection_dx_2,1);
end


