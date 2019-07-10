%% Tabulate Value and Policy Iteration Results, Store to Mat, Graph Results
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_az_vf_post(varargin)
%% FF_AZ_VF_POST post ff_az_vf graphs, tables, mats.
% Given the solution form ff_az_vf, graphs, tables, mats. Graphing code is
% separately in ff_az_vf_post_graph.m. Run this function directly with
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
%    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
%    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
%    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);
%    result_map = ff_az_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_vf_post_graph.html ff_az_vf_post_graph>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html ffs_az_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_get_funcgrid.html ffs_az_get_funcgrid>
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

    params_group = values(result_map, {'mt_val', 'cl_mt_pol_a'});
    [mt_val, cl_mt_pol_a] = params_group{:};
    mt_pol_a = deal(cl_mt_pol_a{1});
    
    params_group = values(result_map, {'ar_val_diff_norm', 'ar_pol_diff_norm', 'mt_pol_perc_change'});
    [ar_val_diff_norm, ar_pol_diff_norm, mt_pol_perc_change] = params_group{:};

    % Get Parameters
    params_group = values(param_map, {'it_z_n', 'it_a_n'});
    [it_z_n, it_a_n] = params_group{:};
    params_group = values(armt_map, {'ar_a'});
    [ar_a] = params_group{:};

else
    clear all;
    close all;

    % internal invoke for testing
    it_param_set = 4;
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
    mt_pol_a = zeros(size(mt_val)) + ar_a'*(cumsum(sort(ar_z))/sum(ar_z)*0.4 + 0.4);
    it_iter_max = min(it_maxiter_val, 50);
    ar_val_diff_norm = rand([it_iter_max, 1]);
    ar_pol_diff_norm = rand([it_iter_max, 1]);
    mt_pol_perc_change = rand([it_iter_max, it_z_n]);

    % Set Results Map
    result_map = containers.Map('KeyType','char', 'ValueType','any');
    result_map('mt_val') = mt_val;
    result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
end

%% Parse Parameter

% param_map
params_group = values(param_map, {'st_model'});
[st_model] = params_group{:};

% armt_map
if (strcmp(st_model, 'az'))
    params_group = values(armt_map, {'ar_z'});
    [ar_z] = params_group{:};
elseif (ismember(st_model, ['abz', 'abz_fibs']))
    params_group = values(armt_map, {'ar_z_r_borr_mesh_wage', 'ar_z_wage_mesh_r_borr'});
    [ar_z_r_borr_mesh_wage, ar_z_wage_mesh_r_borr] = params_group{:};
    params_group = values(param_map, {'fl_z_r_borr_n'});
    [fl_z_r_borr_n] = params_group{:};        
end

% support_map
params_group = values(support_map, {'bl_display_final', 'it_display_final_rowmax', 'it_display_final_colmax'});
[bl_display_final, it_display_final_rowmax, it_display_final_colmax] = params_group{:};
params_group = values(support_map, {'bl_graph', 'bl_graph_onebyones'});
[bl_graph] = params_group{:};
params_group = values(support_map, {'bl_mat', 'st_mat_path', 'st_mat_prefix', 'st_mat_name_main', 'st_mat_suffix'});
[bl_mat, st_mat_path, st_mat_prefix, st_mat_name_main, st_mat_suffix] = params_group{:};

%% Generate Consumption and Income Matrix

if (~isKey(result_map, 'cl_mt_pol_c'))
    f_cons = func_map('f_cons');
    mt_cons = f_cons(ar_z, ar_a', mt_pol_a);
    result_map('cl_mt_pol_c') = {mt_cons, zeros(1)};
end
if (~isKey(result_map, 'cl_mt_coh'))
    f_coh = func_map('f_coh');
    mt_coh = f_coh(ar_z, ar_a');
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
    ff_az_vf_post_graph(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
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
    if (strcmp(st_model, 'az'))
        ar_st_col_zs = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z(ar_it_cols))));
    elseif (ismember(st_model, ['abz', 'abz_fibs']))
        if (fl_z_r_borr_n == 1)
            ar_st_col_zs = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z_wage_mesh_r_borr(ar_it_cols))));
        else
            ar_st_col_zs = matlab.lang.makeValidName(strcat('zi', string(ar_it_cols), ...
                                                            ':zr=', string(ar_z_r_borr_mesh_wage(ar_it_cols)), ...
                                                            ';zw=', string(ar_z_wage_mesh_r_borr(ar_it_cols))));                          
        end
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
                                     mt_pol_perc_change(ar_it_rows_iter, ar_it_cols)']');

    cl_col_names = ['valgap', 'polgap', ar_st_col_zs];
    cl_row_names = strcat('iter=', string(ar_it_rows_iter));
    tb_valpol_alliter.Properties.VariableNames = cl_col_names;
    tb_valpol_alliter.Properties.RowNames = cl_row_names;
    tb_valpol_alliter.Properties.VariableDescriptions{'valgap'} = 'norm(mt_val - mt_val_cur)';
    tb_valpol_alliter.Properties.VariableDescriptions{'polgap'} = 'norm(mt_pol_a - mt_pol_a_cur)';

    disp('valgap = norm(mt_val - mt_val_cur): value function difference across iterations');
    disp('polgap = norm(mt_pol_a - mt_pol_a_cur): policy function difference across iterations');
    disp(['z1 = z1 perc change: sum((mt_pol_a ~= mt_pol_a_cur))/(it_a_n): percentage of state space'...
          ' points conditional on shock where the policy function is changing across iterations']);
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
    ar_it_rows = unique(ar_it_rows);
    
    % Subsetting
    mt_val_print = mt_val(ar_it_rows, ar_it_cols);
    mt_pol_a_print = mt_pol_a(ar_it_rows, ar_it_cols);
        
    % Display Optimal Values
    tb_val = array2table(mt_val_print);
    tb_val.Properties.RowNames = strcat('a', string(ar_it_rows), '=', string(ar_a(ar_it_rows)));
    tb_val.Properties.VariableNames = ar_st_col_zs;
    disp('tb_val: V(a,z) value at each state space point');
    disp(tb_val);

    % Display Optimal Choices
    tb_pol_a = array2table(mt_pol_a_print);
    tb_pol_a.Properties.RowNames = strcat('a', string(ar_it_rows), '=', string(ar_a(ar_it_rows)));
    tb_pol_a.Properties.VariableNames = ar_st_col_zs;
    disp('tb_pol_a: optimal asset choice for each state space point');
    disp(tb_pol_a);

    % Save to result map
    result_map('tb_valpol_alliter') = tb_valpol_alliter;
    result_map('tb_val') = tb_val;
    result_map('tb_pol_a') = tb_pol_a;
end

end
