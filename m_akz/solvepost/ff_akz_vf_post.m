%%
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

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
%    result_map = ff_akz_vf_post(param_map, support_map, armt_map, func_map, result_map,    bl_input_override);
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solvepost/ff_akz_vf_post_graph.m ff_akz_vf_post_graph>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/paramfunc/ffs_akz_set_default_param.m ffs_akz_set_default_param>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/paramfunc/ffs_akz_get_funcgrid.m ffs_akz_get_funcgrid>
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

    params_group = values(result_map, {'mt_val', 'mt_pol_a', 'mt_pol_k'});
    [mt_val, mt_pol_a, mt_pol_k] = params_group{:};
    params_group = values(result_map, {'ar_val_diff_norm', 'ar_pol_diff_norm', 'mt_pol_perc_change'});
    [ar_val_diff_norm, ar_pol_diff_norm, mt_pol_perc_change] = params_group{:};

    % Get Parameters
    params_group = values(param_map, {'it_z_n'});
    [it_z_n] = params_group{:};
    params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'ar_z', 'mt_coh_wkb', 'it_ameshk_n'});
    [ar_a_meshk, ar_k_mesha, ar_z, mt_coh_wkb, it_ameshk_n] = params_group{:};

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
    result_map('mt_pol_a') = mt_pol_a;
    result_map('mt_pol_k') = mt_pol_k;
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

%% Generate Consumption and Income Matrix

mt_cons = f_cons(mt_coh_wkb, mt_pol_a, mt_pol_k);
mt_incm = f_inc(ar_z, ar_a_meshk, ar_k_mesha);
result_map('mt_cons') = mt_cons;
result_map('mt_incm') = mt_incm;

%% Save Mat

if (bl_mat)
    mkdir(support_map('st_mat_path'));
    st_file_name = [st_mat_prefix st_mat_name_main st_mat_suffix];
    save(strcat(st_mat_path, st_file_name));
end

%% Display Val Pol Iter Table

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
    tb_valpol_alliter.Properties.VariableDescriptions{'z1'} = 'z1 perc change: sum((mt_pol_a ~= mt_pol_a_cur))/(it_ameshk_n)';

    disp('valgap = norm(mt_val - mt_val_cur)');
    disp('polgap = norm(mt_pol_a - mt_pol_a_cur)');
    disp('z1 = z1 perc change: sum((mt_pol_a ~= mt_pol_a_cur))/(it_ameshk_n)');
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
    if (it_z_n >= it_display_final_colmax)
        ar_it_cols = (1:1:round(it_display_final_colmax/2));
        ar_it_cols = [ar_it_cols ((it_z_n)-round(it_display_final_colmax/2)+1):1:(it_z_n)];
    else
        ar_it_cols = 1:1:it_z_n;
    end
    mt_val_print = mt_val(ar_it_rows, ar_it_cols);
    mt_pol_a_print = mt_pol_a(ar_it_rows, ar_it_cols);
    mt_pol_k_print = mt_pol_k(ar_it_rows, ar_it_cols);
    mt_pol_w_print = mt_pol_a_print + mt_pol_k_print;

    % Display Optimal Values
    tb_val = array2table(mt_val_print);
    tb_val.Properties.RowNames = strcat('coh', string(ar_it_rows),...
                                        ':k=', string(ar_a_meshk(ar_it_rows)'),...
                                        ',b=', string(ar_k_mesha(ar_it_rows)'));
    tb_val.Properties.VariableNames = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z(ar_it_cols))));
    disp('tb_val');
    disp(tb_val);

    % Display Optimal Choices for a
    tb_pol_a = array2table(mt_pol_a_print);
    tb_pol_a.Properties.RowNames = strcat('coh', string(ar_it_rows),...
                                        ':k=', string(ar_a_meshk(ar_it_rows)'),...
                                        ',b=', string(ar_k_mesha(ar_it_rows)'));
    tb_pol_a.Properties.VariableNames = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z(ar_it_cols))));
    disp('tb_pol_a');
    disp(tb_pol_a);

    % Display Optimal Choices for k
    tb_pol_k = array2table(mt_pol_k_print);
    tb_pol_k.Properties.RowNames = strcat('coh', string(ar_it_rows),...
                                        ':k=', string(ar_a_meshk(ar_it_rows)'),...
                                        ',b=', string(ar_k_mesha(ar_it_rows)'));
    tb_pol_k.Properties.VariableNames = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z(ar_it_cols))));
    disp('tb_pol_k');
    disp(tb_pol_k);

    % Display Optimal Choices for k+a
    tb_pol_w = array2table(mt_pol_w_print);
    tb_pol_w.Properties.RowNames = strcat('coh', string(ar_it_rows),...
                                        ':k=', string(ar_a_meshk(ar_it_rows)'),...
                                        ',b=', string(ar_k_mesha(ar_it_rows)'));
    tb_pol_w.Properties.VariableNames = matlab.lang.makeValidName(strcat('z', string(ar_it_cols), '=', string(ar_z(ar_it_cols))));
    disp('tb_pol_w');
    disp(tb_pol_w);

    % Save to result map
    result_map('tb_valpol_alliter') = tb_valpol_alliter;
    result_map('tb_val') = tb_val;
    result_map('tb_pol_a') = tb_pol_a;
end

%% Generate and Save Graphs

if (bl_graph)
    bl_input_override = true;
    ff_akz_vf_post_graph(param_map, support_map, armt_map, result_map, bl_input_override);
end

end
