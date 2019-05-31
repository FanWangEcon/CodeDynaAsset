%% 
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository> 
% Table of Content.*

function [armt_map, func_map] = ffs_akz_get_funcgrid(varargin)
%% FFS_AKZ_GET_FUNCGRID get funcs, params, states choices shocks grids
% centralized gateway for retrieving parameters, and solution grids and
% functions. 
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param bl_input_override boolean if true varargin contained param_map and
% support_map fully overrides local default. Local default is not invoked.
% This could be important for speed if this function is getting invoked
% within certain loops. Default is 0.
%
% @return armt_map container container with states, choices and shocks
% grids that are inputs for grid based solution algorithm
%
% @return func_map container container with function handles for
% consumption cash-on-hand etc.
%
% @example
%
%    it_param_set = 2;
%    bl_input_override = true;
%    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);
%    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override);
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/paramfunc/ffs_akz_set_functions.m ffs_akz_set_functions>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/tools/ffto_gen_tauchen_jhl.m ffto_gen_tauchen_jhl>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/tools/fft_gen_grid_loglin.m fft_gen_grid_loglin>
%

%% Default

bl_input_override = 0;
if (length(varargin) == 3)
    bl_input_override = varargin{3};
end
if (bl_input_override)
    % override when called from outside
    [param_map, support_map, ~] = varargin{:};
else
    % default internal run
    [param_map, support_map] = ffs_akz_set_default_param();
    support_map('bl_graph_funcgrids') = true;
    support_map('bl_display_funcgrids') = true;
    default_maps = {param_map, support_map};

    % numvarargs is the number of varagin inputted
    [default_maps{1:length(varargin)}] = varargin{:};
    param_map = [param_map; default_maps{1}];
    support_map = [support_map; default_maps{2}];
end

%% Parse Parameters

params_group = values(param_map, {'it_z_n', 'fl_z_mu', 'fl_z_rho', 'fl_z_sig'});
[it_z_n, fl_z_mu, fl_z_rho, fl_z_sig] = params_group{:};

params_group = values(param_map, {'fl_b_bd', 'fl_a_min', 'fl_a_max', 'bl_loglin', 'fl_loglin_threshold', 'it_a_n'});
[fl_b_bd, fl_a_min, fl_a_max, bl_loglin, fl_loglin_threshold, it_a_n] = params_group{:};

params_group = values(param_map, {'fl_k_min', 'fl_k_max', 'bl_k_loglin', 'fl_k_loglin_threshold', 'it_k_n'});
[fl_k_min, fl_k_max, bl_k_loglin, fl_k_loglin_threshold, it_k_n] = params_group{:};

params_group = values(param_map, {'fl_crra', 'fl_c_min'});
[fl_crra, fl_c_min] = params_group{:};

params_group = values(param_map, {'fl_Amean', 'fl_alpha', 'fl_delta'});
[fl_Amean, fl_alpha, fl_delta] = params_group{:};

params_group = values(param_map, {'fl_r_save', 'fl_r_borr', 'fl_w'});
[fl_r_save, fl_r_borr, fl_w] = params_group{:};

params_group = values(support_map, {'bl_graph_funcgrids', 'bl_display_funcgrids'});
[bl_graph_funcgrids, bl_display_funcgrids] = params_group{:};

%% Get Asset and Choice Grid
% Compared to the az model, we now need to generate the grid for a as well
% as k, and mesh them together for the coh grid. 

if (bl_loglin)
    % C:\Users\fan\M4Econ\asset\grid\ff_grid_loglin.m
    ar_a = fft_gen_grid_loglin(it_a_n, fl_a_max, fl_b_bd, fl_loglin_threshold);
else
    ar_a = linspace(fl_b_bd, fl_a_max, it_a_n);
    ar_a = [0 ar_a];
    ar_a = sort(unique(ar_a));
end

if (bl_k_loglin)
    % C:\Users\fan\M4Econ\asset\grid\ff_grid_loglin.m
    ar_k = fft_gen_grid_loglin(it_k_n, fl_k_max, fl_k_min, fl_k_loglin_threshold);
else
    ar_k = linspace(fl_k_min, fl_k_max, it_k_n);
    ar_k = [0 ar_k];
    ar_k = sort(unique(ar_k));
end

[mt_a, mt_k] = ndgrid(ar_a, ar_k);
ar_a_meshk = mt_a(:);
ar_k_mesha = mt_k(:);

%% Get Shock Grids

[~, mt_z_trans, ar_stationary, ar_z] = ffto_gen_tauchen_jhl(fl_z_mu,fl_z_rho,fl_z_sig,it_z_n);

%% Get Equations

[f_util_log, f_util_crra, f_util_standin, f_prod, f_inc, f_coh, f_cons] = ...
    ffs_akz_set_functions(fl_crra, fl_c_min, fl_Amean, fl_alpha, fl_delta, fl_r_save, fl_r_borr, fl_w);

%% Generate Cash-on-Hand/State Grid
% The endogenous state variable is cash-on-hand, it has it_z_n*it_a_n
% number of points, covering all reachable points when ar_a is the choice
% vector and ar_z is the shock vector. requires inputs from get Asset and
% choice grids, get shock grids, and get equations above.

mt_coh = f_coh(ar_z, ar_a_meshk, ar_k_mesha);

%% Store

armt_map = containers.Map('KeyType','char', 'ValueType','any');
armt_map('ar_a') = ar_a;
armt_map('ar_k') = ar_k;
armt_map('mt_z_trans') = mt_z_trans;
armt_map('ar_stationary') = ar_stationary;
armt_map('ar_z') = ar_z;
armt_map('ar_a_meshk') = ar_a_meshk;
armt_map('ar_k_mesha') = ar_k_mesha;
armt_map('it_ameshk_n') = length(ar_a_meshk);
armt_map('mt_coh') = mt_coh;

func_map = containers.Map('KeyType','char', 'ValueType','any');
func_map('f_util_log') = f_util_log;
func_map('f_util_crra') = f_util_crra;
func_map('f_util_standin') = f_util_standin;
func_map('f_prod') = f_prod;
func_map('f_inc') = f_inc;
func_map('f_coh') = f_coh;
func_map('f_cons') = f_cons;

%% Display

if (bl_display_funcgrids)

    disp('ar_z');
    disp(size(ar_z));
    disp(ar_z);

    disp('mt_z_trans');
    disp(size(mt_z_trans));
    disp(mt_z_trans');    
    
    disp('ar_a_meshk');
    disp(size(ar_a_meshk));
    summary(array2table(ar_a_meshk));
    disp(ar_a_meshk');    

    disp('ar_k_mesha');
    disp(size(ar_k_mesha));
    summary(array2table(ar_k_mesha));
    disp(ar_k_mesha');    
    
    disp('mt_coh');
    disp(size(mt_coh));
    summary(array2table(mt_coh));    
    disp(mt_coh);
    
    param_map_keys = keys(func_map);
    param_map_vals = values(func_map);
    for i = 1:length(func_map)
        st_display = strjoin(['pos =' num2str(i) '; key =' string(param_map_keys{i}) '; val =' func2str(param_map_vals{i})]);
        disp(st_display);
    end
    
end

end
