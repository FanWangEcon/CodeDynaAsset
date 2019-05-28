%% FF_AZ_VF_POST handles post ff_az_vf tabling, graphing, gen mat
% While invoking ff_az_vf, if bl_post == true, proceed with below
function [result_map] = ff_az_vf_post(varargin)

%% Parameters Defaults
params_len = length(varargin);
bl_input_override = 0;
if (params_len == 6)
    bl_input_override = varargin{6};
end

if (bl_input_override)
    % if invoked from outside overrid fully
    [param_map, support_map, armt_map, func_map, result_map, ~] = varargin{:};
    
    params_group = values(result_map, {'mt_val', 'mt_pol_a'});
    [mt_val, mt_pol_a] = params_group{:};
    params_group = values(result_map, {'ar_val_diff_norm', 'ar_pol_diff_norm', 'mt_pol_perc_change'});
    [ar_val_diff_norm, ar_pol_diff_norm, mt_pol_perc_change] = params_group{:};
    
    % Get Parameters
    params_group = values(param_map, {'it_z_n', 'it_a_n'});
    [it_z_n, it_a_n] = params_group{:};
    params_group = values(armt_map, {'ar_a', 'ar_z'});
    [ar_a, ar_z] = params_group{:};
    
else
    clear all;
    close all;
    
    % internal invoke for testing
    it_param_set = 2;
    bl_input_override = true;
    
    % Get Parameters
    [param_map, support_map] = ffs_az_set_default_param(it_param_set);
    [armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
    
    % Generate Default val and policy matrixes
    params_group = values(param_map, {'it_maxiter_val', 'it_z_n', 'it_a_n'});
    [it_maxiter_val, it_z_n, it_a_n] = params_group{:};
    params_group = values(armt_map, {'ar_a', 'ar_z'});
    [ar_a, ar_z] = params_group{:};
    f_util_standin = func_map('f_util_standin');
    
    % Set Defaults
    mt_val = f_util_standin(ar_z, ar_a');
    mt_pol_a = zeros(size(mt_val)) + ar_a'*0.9;
    it_iter_max = min(it_maxiter_val, 50);
    ar_val_diff_norm = rand([it_iter_max, 1]);
    ar_pol_diff_norm = rand([it_iter_max, 1]);
    mt_pol_perc_change = rand([it_iter_max, it_z_n]);
    
    % Set Results Map
    result_map = containers.Map('KeyType','char', 'ValueType','any');
    result_map('mt_val') = mt_val;
    result_map('mt_pol_a') = mt_pol_a;
end

%% Parse Parameter
% support_map
params_group = values(support_map, {'bl_display_final', 'it_display_final_rowmax', 'it_display_final_colmax'});
[bl_display_final, it_display_final_rowmax, it_display_final_colmax] = params_group{:};
params_group = values(support_map, {'bl_graph', 'bl_graph_onebyones'});
[bl_graph] = params_group{:};
params_group = values(support_map, {'bl_mat', 'st_mat_path', 'st_mat_prefix', 'st_mat_name_main', 'st_mat_suffix'});
[bl_mat, st_mat_path, st_mat_prefix, st_mat_name_main, st_mat_suffix] = params_group{:};
% func_map
params_group = values(func_map, {'f_inc', 'f_cons', 'f_coh'});
[f_inc, f_cons, f_coh] = params_group{:};

%% Consumption and Income Matrix
mt_cons = f_cons(ar_z, ar_a', mt_pol_a);
mt_incm = f_inc(ar_z, ar_a');
mt_coh = f_coh(ar_z, ar_a');
result_map('mt_cons') = mt_cons;
result_map('mt_incm') = mt_incm;
result_map('mt_coh') = mt_coh;

%% Saving workspace to mat
if (bl_mat)
    mkdir(support_map('st_mat_path'));
    st_file_name = [st_mat_prefix st_mat_name_main st_mat_suffix];
    save(strcat(st_mat_path, st_file_name));
end

%% Table Generation
if (bl_display_final)
    
    % Display Value Function Iteration Step by Step REsults
    it_iter_max = length(ar_val_diff_norm);
    if (it_iter_max >= it_display_final_rowmax)
        ar_it_rows_iter = (1:1:round(it_display_final_rowmax/2));
        ar_it_rows_iter = [ar_it_rows_iter ((it_iter_max)-round(it_display_final_rowmax/2)+1):1:(it_iter_max)];        
    else
        ar_it_rows_iter = 1:1:it_iter_max;
    end
    tb_valpol_alliter = array2table([ar_val_diff_norm(ar_it_rows_iter)';...
                                     ar_pol_diff_norm(ar_it_rows_iter)';...
                                     mt_pol_perc_change(ar_it_rows_iter, :)']');
                                 
    cl_col_names = ['valgap', 'polgap', strcat('z', string((1:it_z_n)))];
    cl_row_names = strcat('iter=', string(ar_it_rows_iter));
    tb_valpol_alliter.Properties.VariableNames = cl_col_names;
    tb_valpol_alliter.Properties.RowNames = cl_row_names;
    tb_valpol_alliter.Properties.VariableDescriptions{'valgap'} = 'norm(mt_val - mt_val_cur)';
    tb_valpol_alliter.Properties.VariableDescriptions{'polgap'} = 'norm(mt_pol_a - mt_pol_a_cur)';
    tb_valpol_alliter.Properties.VariableDescriptions{'z1'} = 'z1 perc change: sum((mt_pol_a ~= mt_pol_a_cur))/(it_a_n)';
    
    disp('valgap = norm(mt_val - mt_val_cur)');
    disp('polgap = norm(mt_pol_a - mt_pol_a_cur)');
    disp('z1 = z1 perc change: sum((mt_pol_a ~= mt_pol_a_cur))/(it_a_n)');   
    disp(tb_valpol_alliter);
    
    % Display Values by States
    % at most display 11 columns of shocks
    % at most display 50 rows for states
    % display first and last
    if (it_a_n >= it_display_final_rowmax)
        ar_it_rows = (1:1:round(it_display_final_rowmax/2));
        ar_it_rows = [ar_it_rows ((it_a_n)-round(it_display_final_rowmax/2)+1):1:(it_a_n)];
    else
        ar_it_rows = 1:1:it_a_n;
    end
    if (it_z_n >= it_display_final_colmax)
        ar_it_cols = (1:1:round(it_display_final_colmax/2));
        ar_it_cols = [ar_it_cols ((it_z_n)-round(it_display_final_colmax/2)+1):1:(it_z_n)];
    else
        ar_it_cols = 1:1:it_z_n;
    end
    mt_val_print = mt_val(ar_it_rows, ar_it_cols);
    mt_pol_a_print = mt_pol_a(ar_it_rows, ar_it_cols);
    
    % Display Optimal Values
    tb_val = array2table(mt_val_print);
    tb_val.Properties.RowNames = strcat('a', string(ar_it_rows), '=', string(ar_a(ar_it_rows)));
    tb_val.Properties.VariableNames = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z(ar_it_cols))));
    disp('tb_val');
    disp(tb_val);
    
    % Display Optimal Choices
    tb_pol_a = array2table(mt_pol_a_print);
    tb_pol_a.Properties.RowNames = strcat('a', string(ar_it_rows), '=', string(ar_a(ar_it_rows)));
    tb_pol_a.Properties.VariableNames = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z(ar_it_cols))));
    disp('tb_pol_a');
    disp(tb_pol_a);
    
    % Save to result map
    result_map('tb_valpol_alliter') = tb_valpol_alliter;
    result_map('tb_val') = tb_val;
    result_map('tb_pol_a') = tb_pol_a;
end

%% Graphing
if (bl_graph)
    ff_az_vf_post_graph(param_map, support_map, armt_map, result_map);
end

end
