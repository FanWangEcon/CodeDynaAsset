%% Tabulate Value and Policy Iteration Results, Store to Mat, Graph Results (Risky + Safe Asset)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_akz_vf_post(varargin)
%% FF_AKZ_VF_POST post ff_akz_vf graphs, tables, mats.
% Given the solution form ff_akz_vf, graphs, tables, mats. Graphing code is
% separately in ff_akz_vf_post_graph.m. Run this function directly with
% randomly generates matrixes for graph/table testings.
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param armt_map container container with states, choices and shocks
% grids that are inputs for grid based solution algorithm
%
% @param func_map container container with function handles for
% consumption cash-on-hand etc.
%
% @param result_map container contains policy function matrix, value
% function matrix, iteration results
%
% @return result_map container add coh consumption and other matrixes to
% result_map also add table versions of val pol and iter matries
%
% @example
%
%    bl_input_override = true;
%    result_map = containers.Map('KeyType','char', 'ValueType','any');
%    result_map('mt_val') = mt_val;
%    result_map('mt_pol_a') = mt_pol_a;
%    result_map('mt_pol_k') = mt_pol_k;
%    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
%    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
%    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);
%    result_map = ff_akz_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solvepost/html/ff_akz_vf_post_graph.html ff_akz_vf_post_graph>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html ffs_akz_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html ffs_akz_get_funcgrid>
%

%% Default

params_len = length(varargin);
bl_input_override = 0;
if (params_len == 6)
    bl_input_override = varargin{6};
end

if (bl_input_override)
    % if invoked from outside overrid fully
    [param_map, support_map, armt_map, func_map, result_map, ~] = varargin{:};

    params_group = values(result_map, {'mt_val', 'cl_mt_pol_a', 'cl_mt_pol_k'});    
    [mt_val, cl_mt_pol_a, cl_mt_pol_k] = params_group{:};
    [mt_pol_a, mt_pol_k] = deal(cl_mt_pol_a{1}, cl_mt_pol_k{1});
    
    params_group = values(result_map, {'ar_val_diff_norm', 'ar_pol_diff_norm', 'mt_pol_perc_change'});
    [ar_val_diff_norm, ar_pol_diff_norm, mt_pol_perc_change] = params_group{:};

    % Get Parameters
    params_group = values(param_map, {'it_z_n'});
    [it_z_n] = params_group{:};
    params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'mt_coh_wkb', 'it_ameshk_n'});
    [ar_a_meshk, ar_k_mesha, mt_coh_wkb, it_ameshk_n] = params_group{:};

else
    clear all;
    close all;

    % internal invoke for testing
    it_param_set = 4;
    bl_input_override = true;

    % Get Parameters
    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override

    % Generate Default val and policy matrixes
    params_group = values(param_map, {'it_maxiter_val', 'it_z_n'});
    [it_maxiter_val, it_z_n] = params_group{:};
    params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'ar_z', 'mt_coh_wkb', 'it_ameshk_n'});
    [ar_a_meshk, ar_k_mesha, ar_z, mt_coh_wkb, it_ameshk_n] = params_group{:};
    params_group = values(func_map, {'f_util_standin', 'f_cons', 'f_coh'});
    [f_util_standin, f_cons, f_coh] = params_group{:};

    % Set Defaults
    mt_val = f_util_standin(ar_z, ar_a_meshk, ar_k_mesha);
    mt_pol_aksum = mt_coh_wkb.*(cumsum(sort(ar_z))/sum(ar_z)*0.4 + 0.4);
    mt_pol_a = mt_pol_aksum.*(0.7 - cumsum(sort(ar_z))/sum(ar_z)*0.3);
    mt_pol_k = mt_pol_aksum - mt_pol_a;
    it_iter_max = min(it_maxiter_val, 50);
    ar_val_diff_norm = rand([it_iter_max, 1]);
    ar_pol_diff_norm = rand([it_iter_max, 1]);
    mt_pol_perc_change = rand([it_iter_max, it_z_n]);

    % Set Results Map
    result_map = containers.Map('KeyType','char', 'ValueType','any');
    result_map('mt_val') = mt_val;
    result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
    result_map('cl_mt_pol_k') = {mt_pol_k, zeros(1)};
    
end

%% Parse Parameter

% Get Parameters
params_group = values(param_map, {'st_model'});
[st_model] = params_group{:};

% armt_map
if (strcmp(st_model, 'ipwkbzr'))
    params_group = values(armt_map, {'ar_z_r_borr_mesh_wage_w1r2', 'ar_z_wage_mesh_r_borr_w1r2'});
    [ar_z_r_borr_mesh_wage_w1r2, ar_z_wage_mesh_r_borr_w1r2] = params_group{:};
    params_group = values(param_map, {'fl_z_r_borr_n'});
    [fl_z_r_borr_n] = params_group{:};        
else
    params_group = values(armt_map, {'ar_z'});
    [ar_z] = params_group{:};    
end

% support_map
params_group = values(support_map, {'bl_display_final', 'it_display_final_rowmax', 'it_display_final_colmax'});
[bl_display_final, it_display_final_rowmax, it_display_final_colmax] = params_group{:};
params_group = values(support_map, {'bl_graph', 'bl_graph_onebyones'});
[bl_graph] = params_group{:};
params_group = values(support_map, {'bl_mat', 'st_mat_path', 'st_mat_prefix', 'st_mat_name_main', 'st_mat_suffix'});
[bl_mat, st_mat_path, st_mat_prefix, st_mat_name_main, st_mat_suffix] = params_group{:};

%% Generate Consumption and Income Matrix
if (~isKey(result_map, 'cl_mt_cons'))
    f_cons = func_map('f_cons');
    mt_cons = f_cons(mt_coh_wkb, mt_pol_a, mt_pol_k);
    result_map('cl_mt_cons') = {mt_cons, zeros(1)};
end
if (~isKey(result_map, 'cl_mt_coh'))
    f_coh = func_map('f_coh');
    mt_coh = f_coh(ar_z, ar_a_meshk, ar_k_mesha);
    result_map('cl_mt_coh') = {mt_coh, zeros(1)};
end

%% Save Mat

if (bl_mat)
    mkdir(support_map('st_mat_path'));
    st_file_name = [st_mat_prefix st_mat_name_main st_mat_suffix];
    save(strcat(st_mat_path, st_file_name));
end

%% Generate and Save Graphs

if (bl_graph)
    bl_input_override = true;
    ff_akz_vf_post_graph(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
end

%% Display Val Pol Iter Table

if (bl_display_final)
    
    % Columns to Display    
    if (it_z_n >= it_display_final_colmax)
        ar_it_cols = (1:1:round(it_display_final_colmax/2));
        ar_it_cols = [ar_it_cols ((it_z_n)-round(it_display_final_colmax/2)+1):1:(it_z_n)];
    else
        ar_it_cols = 1:1:it_z_n;
    end
    ar_it_cols = unique(ar_it_cols);    

    % Column Z Names
    if (strcmp(st_model, 'ipwkbzr'))
        if (fl_z_r_borr_n == 1)
            ar_st_col_zs = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z_wage_mesh_r_borr_w1r2(ar_it_cols))));
        else
            ar_st_col_zs = matlab.lang.makeValidName(strcat('zi', string(ar_it_cols), ...
                                                            ':zr=', string(ar_z_r_borr_mesh_wage_w1r2(ar_it_cols)), ...
                                                            ';zw=', string(ar_z_wage_mesh_r_borr_w1r2(ar_it_cols))));
        end
    else        
        ar_st_col_zs = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z(ar_it_cols))));        
    end
    
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
    tb_valpol_alliter.Properties.VariableDescriptions{'z1'} = 'z1 perc change: sum((mt_pol_a ~= mt_pol_a_cur))/(it_ameshk_n)';

    disp('valgap = norm(mt_val - mt_val_cur): value function difference across iterations');
    disp('polgap = norm(mt_pol_a - mt_pol_a_cur): policy function difference across iterations');
    disp(['z1 = z1 perc change: (sum((mt_pol_a ~= mt_pol_a_cur))+sum((mt_pol_k ~= mt_pol_k_cur)))/(2*it_ameshk_n):' ...
          'percentage of state space points conditional on shock where the policy function is changing across iterations']);
    disp(tb_valpol_alliter);

    % Display Values by States
    % at most display 11 columns of shocks
    % at most display 50 rows for states
    % display first and last
    if (it_ameshk_n >= it_display_final_rowmax)
        ar_it_rows = (1:1:round(it_display_final_rowmax/2));
        ar_it_rows = [ar_it_rows ((it_ameshk_n)-round(it_display_final_rowmax/2)+1):1:(it_ameshk_n)];
    else
        ar_it_rows = 1:1:it_ameshk_n;
    end
    mt_val_print = mt_val(ar_it_rows, ar_it_cols);
    mt_pol_a_print = mt_pol_a(ar_it_rows, ar_it_cols);
    mt_pol_k_print = mt_pol_k(ar_it_rows, ar_it_cols);
    mt_pol_w_print = mt_pol_a_print + mt_pol_k_print;
    
    % Row Name Define
    ar_st_row_coh = strcat('coh', string(ar_it_rows), ...
                            ':k=', string(ar_a_meshk(ar_it_rows)'), ...
                            ',b=', string(ar_k_mesha(ar_it_rows)'));
                                    
    % Display Optimal Values
    tb_val = array2table(mt_val_print);
    tb_val.Properties.RowNames = ar_st_row_coh;
    tb_val.Properties.VariableNames = ar_st_col_zs;
    disp('tb_val: V(a,z) value at each state space point');
    disp(tb_val);

    % Display Optimal Choices for a
    tb_pol_a = array2table(mt_pol_a_print);
    tb_pol_a.Properties.RowNames = ar_st_row_coh;
    tb_pol_a.Properties.VariableNames = ar_st_col_zs;
    disp('tb_pol_a: optimal safe savings choice for each state space point');
    disp(tb_pol_a);

    % Display Optimal Choices for k
    tb_pol_k = array2table(mt_pol_k_print);
    tb_pol_k.Properties.RowNames = ar_st_row_coh;
    tb_pol_k.Properties.VariableNames = ar_st_col_zs;
    disp('tb_pol_k: optimal risky investment choice for each state space point');
    disp(tb_pol_k);

    % Display Optimal Choices for k+a
    tb_pol_w = array2table(mt_pol_w_print);
    tb_pol_w.Properties.RowNames = ar_st_row_coh;
    tb_pol_w.Properties.VariableNames = ar_st_col_zs;
    disp('tb_pol_w: risky + safe investment choices (first stage choice, choose within risky vs safe)');
    disp(tb_pol_w);

    % Save to result map
    result_map('tb_valpol_alliter') = tb_valpol_alliter;
    result_map('tb_val') = tb_val;
    result_map('tb_pol_a') = tb_pol_a;
    
end

end
