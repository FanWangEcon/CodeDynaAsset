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

params_group = values(param_map, {'fl_b_bd', 'fl_w_min', 'fl_w_max', 'it_w_n', 'it_coh_n'});
[fl_b_bd, fl_w_min, fl_w_max, it_w_n, it_coh_n] = params_group{:};

params_group = values(param_map, {'fl_k_min', 'fl_k_max', 'it_ak_n'});
[fl_k_min, fl_k_max, it_ak_n] = params_group{:};

params_group = values(param_map, {'fl_crra', 'fl_c_min', 'it_c_interp_grid_gap'});
[fl_crra, fl_c_min, it_c_interp_grid_gap] = params_group{:};

params_group = values(param_map, {'fl_Amean', 'fl_alpha', 'fl_delta'});
[fl_Amean, fl_alpha, fl_delta] = params_group{:};

params_group = values(param_map, {'fl_r_save', 'fl_r_borr', 'fl_w'});
[fl_r_save, fl_r_borr, fl_w] = params_group{:};

params_group = values(support_map, {'bl_graph_funcgrids', 'bl_display_funcgrids'});
[bl_graph_funcgrids, bl_display_funcgrids] = params_group{:};

%% Get Asset and Choice Grid
% This generate triangular choice structure. Household choose total
% aggregate savings, and within that how much to put into risky capital and
% how much to put into safe assets. See
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html
% ffs_akz_set_default_param> for details.

ar_w = linspace(fl_w_min, fl_w_max, it_w_n);
ar_k = linspace(fl_k_min, fl_k_max, it_ak_n);

mt_a = ar_w - ar_k';
mt_k = ar_w - mt_a;

mt_bl_constrained = (mt_a < fl_b_bd); % this generates NAN, some NAN have -Inf Util
mt_a_wth_na = mt_a;
mt_k_wth_na = mt_k;
mt_a_wth_na(mt_bl_constrained) = nan;
mt_k_wth_na(mt_bl_constrained) = nan;

ar_a_mw_wth_na = mt_a_wth_na(:);
ar_k_mw_wth_na = mt_k_wth_na(:);

ar_a_meshk = ar_a_mw_wth_na(~isnan(ar_a_mw_wth_na));
ar_k_mesha = ar_k_mw_wth_na(~isnan(ar_k_mw_wth_na));

%% Get Shock Grids

[~, mt_z_trans, ar_stationary, ar_z] = ffto_gen_tauchen_jhl(fl_z_mu,fl_z_rho,fl_z_sig,it_z_n);

%% Get Equations

[f_util_log, f_util_crra, f_util_standin, f_prod, f_inc, f_coh, f_cons] = ...
    ffs_akz_set_functions(fl_crra, fl_c_min, fl_Amean, fl_alpha, fl_delta, fl_r_save, fl_r_borr, fl_w);

%% Generate Cash-on-Hand/State Matrix
% The endogenous state variable is cash-on-hand, it has it_z_n*it_a_n
% number of points, covering all reachable points when ar_a is the choice
% vector and ar_z is the shock vector. requires inputs from get Asset and
% choice grids, get shock grids, and get equations above.

mt_coh = f_coh(ar_z, ar_a_meshk, ar_k_mesha);


%% IWKZ Interpolation Cash-on-hand Interpolation Grid
% For the iwkz problems, we solve the problem along a grid of cash-on-hand
% values, the interpolate to find v(k',b',z) at (k',b') choices. 

fl_max_mt_coh = max(max(mt_coh));
fl_min_mt_coh = min(min(mt_coh));
ar_interp_coh_grid = linspace(fl_min_mt_coh, fl_max_mt_coh, it_coh_n);

%% IWKZ Interpolation Consumption Grid
% We also interpolate over consumption to speed the program up. We only
% solve for u(c) at this grid for the iwkz problmes, and then interpolate
% other c values. 

fl_c_max = max(max(mt_coh)) - fl_b_bd;
it_interp_c_grid_n = (fl_c_max-fl_c_min)/(it_c_interp_grid_gap);
ar_interp_c_grid = linspace(fl_c_min, fl_c_max, it_interp_c_grid_n);

%% Store

armt_map = containers.Map('KeyType','char', 'ValueType','any');
armt_map('ar_w') = ar_w;
armt_map('ar_k') = ar_k;
armt_map('mt_z_trans') = mt_z_trans;
armt_map('ar_stationary') = ar_stationary;
armt_map('ar_z') = ar_z;
armt_map('mt_k_wth_na') = mt_k_wth_na;
armt_map('ar_a_mw_wth_na') = ar_a_mw_wth_na;
armt_map('ar_k_mw_wth_na') = ar_k_mw_wth_na;
armt_map('ar_a_meshk') = ar_a_meshk;
armt_map('ar_k_mesha') = ar_k_mesha;
armt_map('it_ameshk_n') = length(ar_a_meshk);
armt_map('mt_coh') = mt_coh;
armt_map('ar_interp_coh_grid') = ar_interp_coh_grid;
armt_map('ar_interp_c_grid') = ar_interp_c_grid;

func_map = containers.Map('KeyType','char', 'ValueType','any');
func_map('f_util_log') = f_util_log;
func_map('f_util_crra') = f_util_crra;
func_map('f_util_standin') = f_util_standin;
func_map('f_prod') = f_prod;
func_map('f_inc') = f_inc;
func_map('f_coh') = f_coh;
func_map('f_cons') = f_cons;

%% Graph

if (bl_graph_funcgrids)

    figure('PaperPosition', [0 0 7 4]);
    hold on;

    chart = plot(mt_a_wth_na, mt_k_wth_na, 'blue');

    clr = jet(numel(chart));
    for m = 1:numel(chart)
       set(chart(m),'Color',clr(m,:))
    end
    if (length(ar_w) <= 100) 
        scatter(ar_a_meshk, ar_k_mesha, 5, 'filled');
    end
    xline(0);
    yline(0);

    title('Choice Grids Conditional on k+a=w')
    ylabel('Capital Choice')
    xlabel({'Borrowing or Saving'})
    legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
    legendCell = cellstr(num2str(ar_w', 'k+a=%3.2f'));
    legend(chart(legend2plot), legendCell(legend2plot), 'Location','northeast');

    grid on;

end

%% Display

if (bl_display_funcgrids)

    disp('ar_z');
    disp(size(ar_z));
    disp(ar_z);

    disp('ar_w');
    disp(size(ar_w));
    disp(ar_w');    
    
    disp('mt_z_trans');
    disp(size(mt_z_trans));
    disp(mt_z_trans');    
    
    disp('ar_interp_coh_grid');
    disp(size(ar_interp_coh_grid));
    summary(array2table(ar_interp_coh_grid'));
    
    disp('ar_interp_c_grid');
    disp(size(ar_interp_c_grid));
    summary(array2table(ar_interp_c_grid'));    
        
    disp('mt_a_wth_na');
    disp(size(mt_a_wth_na));
%     summary(array2table(mt_a_wth_na));
    disp(mt_a_wth_na);    

    disp('ar_a_meshk');
    disp(size(ar_a_meshk));
    summary(array2table(ar_a_meshk));
    disp(ar_a_meshk');    

    disp('mt_k_wth_na');
    disp(size(mt_k_wth_na));
%     summary(array2table(mt_k_wth_na));
    disp(mt_k_wth_na);    
    
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
