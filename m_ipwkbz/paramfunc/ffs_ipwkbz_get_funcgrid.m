%% Generate States/Choices/Shocks Grids, get Functions (Interpolated + Percentage + Risky + Safe Asset + Save + Borrow)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [armt_map, func_map] = ffs_ipwkbz_get_funcgrid(varargin)
%% FFS_IPWKBZ_GET_FUNCGRID get funcs, params, states choices shocks grids
% centralized gateway for retrieving parameters, and solution grids and
% functions. Similar to
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html
% ffs_akz_get_funcgrid> function. This code deals with problems with
% savings and borrowing.
%
% This file is a continuation of
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkz/paramfunc/html/ffs_ipwkz_get_funcgrid.html
% ffs_ipwkz_get_funcgrid>. A significant change here is the inclusion of
% borrowing shocks, which changes significantly the implementation of the
% two-stage solution structure. Specifically, for the coh matrix here,
% there is another dimension which is the borrowing interest rate shock
% states, which expands the number of cash-on-hand points given the same
% productivity shock z. So compared to <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkz/paramfunc/html/ffs_ipwkz_get_funcgrid.html
% ffs_ipwkz_get_funcgrid>, the same number of columns for the
% _mt_coh_wkb_full_ matirx, but more rows.
%
% Note that for borrowing, we can not start at the min(coh(k,w-k,z))
% reacheable given the w and k choice grids. That would be:
% min(coh(k,w-k,z)) < min(w). At min(coh), there is a single point, and
% also for all values between min(coh) and min(w), all percentage w_perc
% choices are below min(w), requiring extrapolation. It is in some sense
% safer to extrapolate by starting solution at min(coh) = min(w). Now, for
% the lowest levels of coh(k, w-k, z), which are < min(w), we now require
% extrapolation. But this extrapolation is straight forward, when coh <
% min(w), households have to default, the utility from default is the same
% regardless of the level of coh at the time of default. So nearest
% extrapolation is fully correct.
%
% Note that the first stage w grid is based on cash-on-hand level reached
% by the coh(k,w-k,z) possible choice and shock combinations. This
% coh(k,w-k,z) > max(w), which also means that at max(coh) grid, the w_perc
% choices at higher points require extrapolation. Extrapolation is based on
% nearest extrapolation.
%
% Note that even when w = 0, as long as interest rate is low, only the
% lowest level of borrowing is invalid.
%
% All discussion of of the word _wage_ refers to productivity shock.
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
% grids that are inputs for grid based solution algorithm. Contains these:
%
% *armt_map* base arrays:
%
% # ar_interp_c_grid: 1 by I^c
% # ar_interp_coh_grid: 1 by I^{coh}
% # ar_w_level: 1 by I^{W=k+b}
% # ar_w_perc: 1 by P^{W=k+b}
% # ar_ak_perc: 1 by P^{k and b}
%
% *armt_map* 1st stage level coh on hand related arrays:
%
% # mt_interp_coh_grid_mesh_z_wage: I^{coh} by M^w
% # mt_z_wage_mesh_interp_coh_grid: I^{coh} by M^w
% # mt_interp_coh_grid_mesh_w_perc: I^{coh} by P^{LAM=k+b}
% # mt_w_perc_mesh_interp_coh_grid: I^{coh} by P^{LAM=k+b}
%
% *armt_map* 2nd stage reachable coh(k(w), a(w,k), z', r)
%
% # mt_coh_wkb: (I^k x I^w x M^r) by (M^z)
% # mt_z_wage_mesh_coh_wkb: (I^k x I^w x M^r) by (M^z)
%
% *armt_map* 2nd stage additional arrays
%
% # mt_k: (I^w) by (P^{k and b})
% # ar_a_meshk: 1 by (I^w x P^{k and b})
% # ar_k_mesha: 1 by (I^w x P^{k and b})
% # ar_aplusk_mesh: 1 by (I^w x P^{k and b})
% # it_ameshk_n: scalar
%
% *armt_map* Shock Grids Arrays and Mesh
%
% # ar_z_r_borr: 1 by (M^r)
% # ar_z_r_borr_prob: 1 by (M^r)
% # ar_z_wage: 1 by (M^z)
% # ar_z_wage_prob: 1 by (M^z)
% # ar_z_r_borr_mesh_wage_w1r2: 1 by (M^z x M^r)
% # ar_z_wage_mesh_r_borr_w1r2: 1 by (M^z x M^r)
% # ar_z_r_borr_mesh_wage_r1w2: 1 by (M^r x M^z)
% # ar_z_wage_mesh_r_borr_r1w2: 1 by (M^r x M^z)
%
% @return func_map container container with function handles for
% consumption cash-on-hand etc.
%
% @example
%
%    it_param_set = 2;
%    bl_input_override = true;
%    [param_map, support_map] = ffs_ipwkbz_set_default_param(it_param_set);
%    [armt_map, func_map] = ffs_ipwkbz_get_funcgrid(param_map, support_map, bl_input_override);
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkbz/paramfunc/ffs_ipwkbz_set_functions.m ffs_ipwkbz_set_functions>
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

    close all;

    % default internal run
    [param_map, support_map] = ffs_ipwkbz_set_default_param(4);

    support_map('bl_graph_funcgrids') = true;
    support_map('bl_graph_funcgrids_detail') = true;
    bl_display_funcgrids = true;
    support_map('bl_display_funcgrids') = bl_display_funcgrids;

    st_param_which = 'medium';

    if (ismember(st_param_which, ['default']))

        param_map('it_ak_perc_n') = 250;

    elseif (ismember(st_param_which, ['medium']))

        % to be able to visually see choice grid points
        param_map('fl_b_bd') = -20; % borrow bound, = 0 if save only
        param_map('fl_default_aprime') = 0;
        param_map('bl_default') = 0; % if borrowing is default allowed

        param_map('fl_w_min') = param_map('fl_b_bd');
        param_map('it_w_perc_n') = 25;
        param_map('it_ak_perc_n') = 45;

        param_map('fl_w_interp_grid_gap') = 2;
        param_map('fl_coh_interp_grid_gap') = 2;

%         param_map('fl_z_r_borr_min') = 0.025;
%         param_map('fl_z_r_borr_max') = 0.95;
%         param_map('fl_z_r_borr_n') = 3;

        param_map('fl_z_r_borr_min') = 0.025;
        param_map('fl_z_r_borr_max') = 0.95;
        param_map('fl_z_r_borr_n') = 2;

    elseif (strcmp(st_param_which, 'small'))

        param_map('fl_z_r_borr_n') = 2;
        param_map('it_z_wage_n') = 3;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

        param_map('fl_b_bd') = -20; % borrow bound, = 0 if save only
        param_map('fl_default_aprime') = 0;
        param_map('bl_default') = 0; % if borrowing is default allowed

        param_map('fl_w_min') = param_map('fl_b_bd');
        param_map('it_w_perc_n') = 5;
        param_map('it_ak_perc_n') = 6;

        param_map('fl_w_interp_grid_gap') = 3;
        param_map('fl_coh_interp_grid_gap') = 3;

        param_map('fl_z_r_borr_min') = 0.025;
        param_map('fl_z_r_borr_max') = 0.95;
        param_map('fl_z_r_borr_n') = 3;

    end

    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

    default_maps = {param_map, support_map};

    % numvarargs is the number of varagin inputted
    [default_maps{1:length(varargin)}] = varargin{:};
    param_map = [param_map; default_maps{1}];
    support_map = [support_map; default_maps{2}];

    % Display Parameters
    if (bl_display_funcgrids)
        fft_container_map_display(param_map);
        fft_container_map_display(support_map);
    end

end

%% Parse Parameters 1a

params_group = values(param_map, {'fl_b_bd', 'fl_w_min', 'fl_w_max'});
[fl_b_bd, fl_w_min, fl_w_max] = params_group{:};

params_group = values(param_map, {'fl_crra', 'fl_c_min'});
[fl_crra, fl_c_min] = params_group{:};

params_group = values(param_map, {'fl_Amean', 'fl_alpha', 'fl_delta'});
[fl_Amean, fl_alpha, fl_delta] = params_group{:};

params_group = values(param_map, {'fl_r_save', 'fl_w'});
[fl_r_save, fl_w] = params_group{:};

%% Parse Parameters 1b

params_group = values(param_map, {...
    'it_w_perc_n', 'it_ak_perc_n',...
    'it_c_interp_grid_gap', 'fl_w_interp_grid_gap', 'fl_coh_interp_grid_gap'});
[it_w_perc_n, it_ak_perc_n,...
    it_c_interp_grid_gap, fl_w_interp_grid_gap, fl_coh_interp_grid_gap] = params_group{:};

%% Parse Parameters 2

% param_map shock income
params_group = values(param_map, {'it_z_wage_n', 'fl_z_wage_mu', 'fl_z_wage_rho', 'fl_z_wage_sig'});
[it_z_wage_n, fl_z_wage_mu, fl_z_wage_rho, fl_z_wage_sig] = params_group{:};

% param_map shock borrowing interest
params_group = values(param_map, {'st_z_r_borr_drv_ele_type', 'st_z_r_borr_drv_prb_type', 'fl_z_r_borr_poiss_mean', ...
    'fl_z_r_borr_max', 'fl_z_r_borr_min', 'fl_z_r_borr_n'});
[st_z_r_borr_drv_ele_type, st_z_r_borr_drv_prb_type, fl_z_r_borr_poiss_mean, ...
    fl_z_r_borr_max, fl_z_r_borr_min, fl_z_r_borr_n] = params_group{:};

% param_map shock income
params_group = values(param_map, {'it_z_n'});
[it_z_n] = params_group{:};

%% Parse Parameters 3

params_group = values(support_map, {'bl_graph_funcgrids', 'bl_graph_funcgrids_detail', 'bl_display_funcgrids'});
[bl_graph_funcgrids, bl_graph_funcgrids_detail, bl_display_funcgrids] = params_group{:};

%% Generate Asset and Choice Grid for 2nd stage Problem
% This generate triangular choice structure. Household choose total
% aggregate savings, and within that how much to put into risky capital and
% how much to put into safe assets, in percentages. See
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbz/paramfunc/html/ffs_ipwkbz_set_default_param.html
% ffs_ipwkbz_set_default_param> for details.

% percentage grid for 1st stage choice problem, level grid for 2nd stage
% solving optimal k given w and z.
ar_w_perc = linspace(0.001, 0.999, it_w_perc_n);
it_w_interp_n = ((fl_w_max-fl_w_min)/(fl_w_interp_grid_gap));
ar_w_level_full = fft_array_add_zero(linspace(fl_w_min, fl_w_max, it_w_interp_n), true);
ar_w_level = ar_w_level_full;
it_w_interp_n = length(ar_w_level_full);

% max k given w, need to consider the possibility of borrowing.
ar_k_max = ar_w_level_full - fl_b_bd;

% k percentage choice grid
ar_ak_perc = linspace(0.001, 0.999, it_ak_perc_n);

% 2nd stage percentage choice matrixes
% (ar_k_max') is it_w_interp_n by 1, and (ar_ak_perc) is 1 by it_ak_perc_n
% mt_k is a it_w_interp_n by it_ak_perc_n matrix of choice points of k'
% conditional on w, each column is a different w, each row for each col a
% different k' value.
mt_k = (ar_k_max'*ar_ak_perc)';
mt_a = (ar_w_level_full - mt_k);

% can not have choice that are beyond feasible bound given the percentage
% structure here.
mt_bl_constrained = (mt_a < fl_b_bd);
if (sum(mt_bl_constrained) > 0 )
    error('at %s second stage choice points, percentage choice exceed bounds, can not happen',...
        num2str(sum(mt_bl_constrained)));
end

ar_a_meshk_full = mt_a(:);
ar_k_mesha_full = mt_k(:);

ar_a_meshk = ar_a_meshk_full;
ar_k_mesha = ar_k_mesha_full;

%% Get Shock: Income Shock (ar1)

[~, mt_z_wage_trans, ar_z_wage_prob, ar_z_wage] = ffto_gen_tauchen_jhl(fl_z_wage_mu,fl_z_wage_rho,fl_z_wage_sig,it_z_wage_n);

%% Get Shock: Interest Rate Shock (iid)

% get borrowing grid and probabilities
param_dsv_map = containers.Map('KeyType','char', 'ValueType','any');
param_dsv_map('st_drv_ele_type') = st_z_r_borr_drv_ele_type;
param_dsv_map('st_drv_prb_type') = st_z_r_borr_drv_prb_type;
param_dsv_map('fl_poiss_mean') = fl_z_r_borr_poiss_mean;
param_dsv_map('fl_max') = fl_z_r_borr_max;
param_dsv_map('fl_min') = fl_z_r_borr_min;
param_dsv_map('fl_n') = fl_z_r_borr_n;
[ar_z_r_borr, ar_z_r_borr_prob] = fft_gen_discrete_var(param_dsv_map, true);

% iid transition matrix
mt_z_r_borr_prob_trans = repmat(ar_z_r_borr_prob, [length(ar_z_r_borr_prob), 1]);

%% Get Shock: Mesh Shocks Together
% R is outter, W is Inner

% Kronecker product to get full transition matrix for the two shocks
mt_z_trans = kron(mt_z_r_borr_prob_trans, mt_z_wage_trans);

% mesh the shock vectors
[mt_z_wage_mesh_r_borr_w1r2, mt_z_r_borr_mesh_wage_w1r2] = ndgrid(ar_z_wage, ar_z_r_borr);
ar_z_wage_mesh_r_borr_w1r2 = mt_z_wage_mesh_r_borr_w1r2(:)';
ar_z_r_borr_mesh_wage_w1r2 = mt_z_r_borr_mesh_wage_w1r2(:)';

% mesh the shock vectors
[mt_z_r_borr_mesh_wage_r1w2, mt_z_wage_mesh_r_borr_r1w2] = ndgrid(ar_z_r_borr, ar_z_wage);
ar_z_wage_mesh_r_borr_r1w2 = mt_z_wage_mesh_r_borr_r1w2(:)';
ar_z_r_borr_mesh_wage_r1w2 = mt_z_r_borr_mesh_wage_r1w2(:)';

if (bl_display_funcgrids)

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Borrow R Shock: ar_z_r_borr_mesh_wage_w1r2');
    disp('Prod/Wage Shock: ar_z_wage_mesh_r_borr_w1r2');
    disp('show which shock is inner and which is outter');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');

    tb_two_shocks = array2table([ar_z_r_borr_mesh_wage_w1r2;...
                                 ar_z_wage_mesh_r_borr_w1r2]');
    cl_col_names = ["Borrow R Shock (Meshed)", "Wage R Shock (Meshed)"];
    cl_row_names = strcat('zi=', string((1:it_z_n)));
    tb_two_shocks.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_two_shocks.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);

    it_row_display = it_z_wage_n*2;

    disp(size(tb_two_shocks));
    disp(head(tb_two_shocks, it_row_display));
    disp(tail(tb_two_shocks, it_row_display));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Borrow R Shock: ar_z_wage_mesh_r_borr_r1w2');
    disp('Prod/Wage Shock: ar_z_r_borr_mesh_wage_r1w2');
    disp('show which shock is inner and which is outter');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');

    tb_two_shocks = array2table([ar_z_wage_mesh_r_borr_r1w2;...
                                 ar_z_r_borr_mesh_wage_r1w2]');
    cl_col_names = ["Borrow R Shock (Meshed)", "Wage R Shock (Meshed)"];
    cl_row_names = strcat('zi=', string((1:length(ar_z_r_borr_mesh_wage_r1w2))));
    tb_two_shocks.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_two_shocks.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);

    it_row_display = fl_z_r_borr_n*2;

    disp(size(tb_two_shocks));
    disp(head(tb_two_shocks, it_row_display));
    disp(tail(tb_two_shocks, it_row_display));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Borrow Rate Transition Matrix: mt_z_r_borr_prob_trans');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    it_col_n_keep = 15;
    it_row_n_keep = 15;
    [it_row_n, it_col_n] = size(mt_z_r_borr_prob_trans);
    [ar_it_cols, ar_it_rows] = fft_row_col_subset(it_col_n, it_col_n_keep, it_row_n, it_row_n_keep);
    cl_st_full_rowscols = cellstr([num2str(ar_z_r_borr', 'r%3.2f')]);
    tb_z_r_borr_prob_trans = array2table(round(mt_z_r_borr_prob_trans(ar_it_rows, ar_it_cols), 6));
    cl_col_names = strcat('zi=', num2str(ar_it_cols'), ':', cl_st_full_rowscols(ar_it_cols));
    cl_row_names = strcat('zi=', num2str(ar_it_rows'), ':', cl_st_full_rowscols(ar_it_rows));
    tb_z_r_borr_prob_trans.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_z_r_borr_prob_trans.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);

    disp(size(tb_z_r_borr_prob_trans));
    disp(tb_z_r_borr_prob_trans);

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Wage Prod Shock Transition Matrix: mt_z_r_borr_prob_trans');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    it_col_n_keep = 15;
    it_row_n_keep = 15;
    [it_row_n, it_col_n] = size(mt_z_wage_trans);
    [ar_it_cols, ar_it_rows] = fft_row_col_subset(it_col_n, it_col_n_keep, it_row_n, it_row_n_keep);
    cl_st_full_rowscols = cellstr([num2str(ar_z_wage', 'w%3.2f')]);
    tb_z_wage_trans = array2table(round(mt_z_wage_trans(ar_it_rows, ar_it_cols),6));
    cl_col_names = strcat('zi=', num2str(ar_it_cols'), ':', cl_st_full_rowscols(ar_it_cols));
    cl_row_names = strcat('zi=', num2str(ar_it_rows'), ':', cl_st_full_rowscols(ar_it_rows));
    tb_z_wage_trans.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_z_wage_trans.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);

    disp(size(tb_z_wage_trans));
    disp(tb_z_wage_trans);

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Full Transition Matrix: mt_z_trans');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    it_col_n_keep = it_z_wage_n*3;
    it_row_n_keep = it_z_wage_n*3;
    [it_row_n, it_col_n] = size(mt_z_trans);
    [ar_it_cols, ar_it_rows] = fft_row_col_subset(it_col_n, it_col_n_keep, it_row_n, it_row_n_keep);
    cl_st_full_rowscols = cellstr([num2str(ar_z_r_borr_mesh_wage_w1r2', 'r%3.2f;'), ...
                                   num2str(ar_z_wage_mesh_r_borr_w1r2', 'w%3.2f')]);
    tb_mt_z_trans = array2table(round(mt_z_trans(ar_it_rows, ar_it_cols),6));
    cl_col_names = strcat('i', num2str(ar_it_cols'), ':', cl_st_full_rowscols(ar_it_cols));
    cl_row_names = strcat('i', num2str(ar_it_rows'), ':', cl_st_full_rowscols(ar_it_rows));
    tb_mt_z_trans.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_mt_z_trans.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);

    disp(size(tb_mt_z_trans));
    disp(tb_mt_z_trans);

end

%% Get Equations

[f_util_log, f_util_crra, f_util_standin, f_util_standin_coh, f_prod, f_inc, f_coh, f_cons] = ...
    ffs_ipwkbz_set_functions(fl_crra, fl_c_min, fl_b_bd, fl_Amean, fl_alpha, fl_delta, fl_r_save, fl_w);

%% Generate Cash-on-Hand/State Matrix
% The endogenous state variable is cash-on-hand, it has it_z_n*it_a_n
% number of points, covering all reachable points when ar_a is the choice
% vector and ar_z is the shock vector. requires inputs from get Asset and
% choice grids, get shock grids, and get equations above.
%
% # mt_coh_wkb_full: this is the (I^k x I^w) by (M^r x M^z) matrix, where
% rows = it_w_interp_n*it_ak_perc_n, and cols = fl_z_r_borr_n*it_z_wage_n.
% # mt_coh_wkb_full: this is the (I^k x I^w x M^r) by (M^z) matrix, where
% rows = it_w_interp_n*it_ak_perc_n*fl_z_r_borr_n, and cols = it_z_wage_n.
%
mt_coh_wkb_full = f_coh(ar_z_r_borr_mesh_wage_r1w2, ar_z_wage_mesh_r_borr_r1w2, ...
                        ar_a_meshk_full, ar_k_mesha_full);

it_coh_wkb_reshape_rows = it_w_interp_n*it_ak_perc_n*fl_z_r_borr_n;
it_coh_wkb_reshape_cols = it_z_wage_n;
mt_coh_wkb_full = reshape(mt_coh_wkb_full, [it_coh_wkb_reshape_rows, it_coh_wkb_reshape_cols]);

% Generate Aggregate Variables
ar_aplusk_mesh = ar_a_meshk_full + ar_k_mesha_full;

if (bl_display_funcgrids || bl_graph_funcgrids)

    % Genereate Table
    tab_ak_choices = array2table([ar_aplusk_mesh, ar_k_mesha_full, ar_a_meshk_full]);
    cl_col_names = {'ar_aplusk_mesh', 'ar_k_mesha_full', 'ar_a_meshk_full'};
    tab_ak_choices.Properties.VariableNames = cl_col_names;

    % Label Table Variables
    tab_ak_choices.Properties.VariableDescriptions{'ar_aplusk_mesh'} = ...
        '*ar_aplusk_mesha*: ar_aplusk_mesha = ar_a_meshk_full + ar_k_mesha_full;';
    tab_ak_choices.Properties.VariableDescriptions{'ar_a_meshk_full'} = ...
        '*ar_a_meshk_full*:';
    tab_ak_choices.Properties.VariableDescriptions{'ar_k_mesha_full'} = ...
        '*ar_k_mesha_full*:';

    cl_var_desc = tab_ak_choices.Properties.VariableDescriptions;
    for it_var_name = 1:length(cl_var_desc)
        disp(cl_var_desc{it_var_name});
    end

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('tab_ak_choices');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    it_rows_toshow = length(ar_w_level)*2;
    disp(size(tab_ak_choices));
    disp(head(array2table(tab_ak_choices), it_rows_toshow));
    disp(tail(array2table(tab_ak_choices), it_rows_toshow));

    % Generate Shock Full Mat to See coh_wkb_full zr and zw
    mt_z_r_borr_mwge_makfull = zeros([length(ar_a_meshk_full), it_z_n]) + ar_z_r_borr_mesh_wage_r1w2;
    mt_z_wage_mesh_r_makfull = zeros([length(ar_a_meshk_full), it_z_n]) + ar_z_wage_mesh_r_borr_r1w2;
    % use the same reshape command as above
    mt_z_r_borr_mwge_makfull = reshape(mt_z_r_borr_mwge_makfull, [it_coh_wkb_reshape_rows, it_coh_wkb_reshape_cols]);
    mt_z_wage_mesh_r_makfull = reshape(mt_z_wage_mesh_r_makfull, [it_coh_wkb_reshape_rows, it_coh_wkb_reshape_cols]);

    % Generate W and K Choices Full Mat to see which w and k
    mt_aplusk_mesh = zeros([length(ar_a_meshk_full), it_z_n]) + ar_aplusk_mesh;
    mt_a_meshk_full = zeros([length(ar_a_meshk_full), it_z_n]) + ar_a_meshk_full;
    mt_k_mesha_full = zeros([length(ar_a_meshk_full), it_z_n]) + ar_k_mesha_full;
    % use the same reshape command as above
    mt_aplusk_mesh = reshape(mt_aplusk_mesh, [it_coh_wkb_reshape_rows, it_coh_wkb_reshape_cols]);
    mt_a_meshk_full = reshape(mt_a_meshk_full, [it_coh_wkb_reshape_rows, it_coh_wkb_reshape_cols]);
    mt_k_mesha_full = reshape(mt_k_mesha_full, [it_coh_wkb_reshape_rows, it_coh_wkb_reshape_cols]);

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_coh_wkb_full');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_coh_wkb_full));
    disp(head(array2table(mt_coh_wkb_full), it_rows_toshow));
    disp(tail(array2table(mt_coh_wkb_full), it_rows_toshow));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Shock Borrow R: mt_z_r_borr_mwge_makfull');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_z_r_borr_mwge_makfull));
    disp(head(array2table(mt_z_r_borr_mwge_makfull), it_rows_toshow));
    disp(tail(array2table(mt_z_r_borr_mwge_makfull), it_rows_toshow));
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Shock Productivity: mt_z_wage_mesh_r_makfull');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_z_wage_mesh_r_makfull));
    disp(head(array2table(mt_z_wage_mesh_r_makfull), it_rows_toshow));
    disp(tail(array2table(mt_z_wage_mesh_r_makfull), it_rows_toshow));
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('W = Aprime + Kprime; mt_aplusk_mesh');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_aplusk_mesh));
    disp(head(array2table(mt_aplusk_mesh), it_rows_toshow));
    disp(tail(array2table(mt_aplusk_mesh), it_rows_toshow));
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Kprime: mt_k_mesha_full');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_k_mesha_full));
    disp(head(array2table(mt_k_mesha_full), it_rows_toshow));
    disp(tail(array2table(mt_k_mesha_full), it_rows_toshow));
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Aprime: mt_a_meshk_full');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_a_meshk_full));
    disp(head(array2table(mt_a_meshk_full), it_rows_toshow));
    disp(tail(array2table(mt_a_meshk_full), it_rows_toshow));

end

%% Check if COH is within Borrowing Bounds
% some coh levels are below borrowing bound, can not borrow enough to pay
% debt

mt_bl_coh_wkb_invalid = (mt_coh_wkb_full < fl_b_bd);

% (k,a) invalid if coh(k,a,z) < bd for any z
ar_bl_wkb_invalid = max(mt_bl_coh_wkb_invalid,[], 2);
mt_bl_wkb_invalid = reshape(ar_bl_wkb_invalid, [it_ak_perc_n, it_w_interp_n*fl_z_r_borr_n]);

% find the first w_level choice where some k(w) percent choices are valid?
ar_bl_w_level_invalid =  min(mt_bl_wkb_invalid, [], 1);

% w choices can not be lower than fl_w_level_min_valid. If w choices are
% lower, given the current borrowing interest rate as well as the minimum
% income level in the future, and the maximum borrowing level available
% next period, and given the shock distribution, there exists some state in
% the future when the household when making this choice will be unable to
% borrow sufficiently to maintain positive consumption.
ar_w_level_full_dup = repmat(ar_w_level_full, [1,fl_z_r_borr_n]);
fl_w_level_min_valid = min(ar_w_level_full_dup(~ar_bl_w_level_invalid));

%% Update Valid 2nd stage choice matrix
% ar_w_level = linspace(fl_w_level_min_valid, fl_w_max, it_w_interp_n);
% ar_k_max = ar_w_level - fl_b_bd;
% mt_k = (ar_k_max'*ar_ak_perc)';
% mt_a = (ar_w_level - mt_k);
% ar_a_meshk = mt_a(:);
% ar_k_mesha = mt_k(:);

%% Select only Valid (k(w), a) choices

% mt_coh_wkb = mt_coh_wkb_full(~ar_bl_wkb_invalid, :);
mt_coh_wkb = mt_coh_wkb_full;
mt_z_r_borr_mesh_coh_wkb = repmat(ar_z_r_borr, [size(mt_coh_wkb,1),1]);
mt_z_wage_mesh_coh_wkb = repmat(ar_z_wage, [size(mt_coh_wkb,1),1]);

%% Generate 1st Stage States: Interpolation Cash-on-hand Interpolation Grid
% For the iwkz problems, we solve the problem along a grid of cash-on-hand
% values, the interpolate to find v(k',b',z) at (k',b') choices. Crucially,
% we have to coh matrxies

fl_max_mt_coh = max(max(mt_coh_wkb));

% This is savings only condition
% fl_min_mt_coh = min(min(mt_coh_wkb));

% This could be condition if no defaults are allowed
% fl_min_mt_coh = fl_w_level_min_valid;

% This is borrowing with default or not condition
fl_min_mt_coh = fl_b_bd;

it_coh_interp_n = (fl_max_mt_coh-fl_min_mt_coh)/(fl_coh_interp_grid_gap);
ar_interp_coh_grid = fft_array_add_zero(linspace(fl_min_mt_coh, fl_max_mt_coh, it_coh_interp_n), true);
mt_interp_coh_grid_mesh_w_perc = repmat(ar_interp_coh_grid, [it_w_perc_n, 1]);

[mt_interp_coh_grid_mesh_z_wage, mt_z_wage_mesh_interp_coh_grid] = ndgrid(ar_interp_coh_grid, ar_z_wage);

mt_interp_coh_grid_mesh_z = repmat(ar_interp_coh_grid', [1, it_z_n]);

%% Generate 1st Stage Choices: Interpolation Cash-on-hand Interpolation Grid
% previously, our ar_w was the first stage choice grid, the grid was the
% same for all coh levels. Now, for each coh level, there is a different
% ar_w. ar_interp_coh_grid is (1 by ar_interp_coh_grid) and ar_w_perc is (
% 1 by it_w_perc_n). Conditional on z, each choice matrix is (it_w_perc_n
% by ar_interp_coh_grid). Here we are pre-computing the choice matrix. This
% could be a large matrix if the choice grid is large. This is the matrix
% of aggregate savings choices
%

if (fl_min_mt_coh < 0)
    % borrowing bound is below zero
    mt_w_perc_mesh_interp_coh_grid = ((ar_interp_coh_grid-fl_min_mt_coh)'*ar_w_perc)' + fl_min_mt_coh;
else
    % savings only
    mt_w_perc_mesh_interp_coh_grid = ((ar_interp_coh_grid)'*ar_w_perc)';
end


%% Generate Interpolation Consumption Grid
% We also interpolate over consumption to speed the program up. We only
% solve for u(c) at this grid for the iwkz problmes, and then interpolate
% other c values.

fl_c_max = max(max(mt_coh_wkb_full)) - fl_b_bd;
it_interp_c_grid_n = (fl_c_max-fl_c_min)/(it_c_interp_grid_gap);
ar_interp_c_grid = linspace(fl_c_min, fl_c_max, it_interp_c_grid_n);

%% Initialize armt_map to store, state, choice, shock matrixes

armt_map = containers.Map('KeyType','char', 'ValueType','any');
armtdesc_map = containers.Map('KeyType','char', 'ValueType','any');

%% Store armt_map (1): base arrays
% Dimensions of Various Grids: I for level grid, M for shock grid, P for
% percent grid. Dimensions are:
%
% # ar_interp_c_grid: 1 by I^c
% # ar_interp_coh_grid: 1 by I^{coh}
% # ar_w_level: 1 by I^{W=k+b}
% # ar_w_perc: 1 by P^{W=k+b}
% # ar_ak_perc: 1 by P^{k and b}
%
% more descriptions:
%
% # ar_interp_c_grid: 1 by I^c, 1st stage consumption interpolation
% # ar_interp_coh_grid: 1 by I^{coh}, 1st stage value function V(coh,z)
% # ar_w_level: 1 by I^{W=k+b}, 2nd stage k*(w,z) w grid. 2nd stage, level
% of w over which we solve the optimal percentage k' choices. Need to
% generate interpolant based on this so that we know optimal k* given
% ar_w_perc(coh) in the 1st stage
% # ar_w_perc: 1 by P^{W=k+b}, 1st stage w \in {w_perc(coh)} choice set.
% 1st stage, percentage w choice given coh, at each coh level the number of
% choice points is the same for this problem with
% percentage grid points.
% # ar_ak_perc: 1 by P^{k and b}, 2nd stage k \in {ask_perc(w,z)} set
%

armt_map('ar_interp_c_grid') = ar_interp_c_grid;
armt_map('ar_interp_coh_grid') = ar_interp_coh_grid;
armt_map('ar_w_level') = ar_w_level;
armt_map('ar_w_perc') = ar_w_perc;
armt_map('ar_ak_perc') = ar_ak_perc;

%% Store armt_map (2): 1st stage level coh on hand related arrays
% Dimensions of Various Grids: I for level grid, M for shock grid, P for
% percent grid. Dimensions are:
%
% # mt_interp_coh_grid_mesh_z_wage: I^{coh} by M^w
% # mt_z_wage_mesh_interp_coh_grid: I^{coh} by M^w
% # mt_interp_coh_grid_mesh_w_perc: I^{coh} by P^{LAM=k+b}
% # mt_w_perc_mesh_interp_coh_grid: I^{coh} by P^{LAM=k+b}
%
% more descriptions:
%
% # *mt_w_perc_mesh_interp_coh_grid* 1st stage, generate w(coh, percent),
% meaning the level of w given coh and the percentage grid of ar_w_perc.
% Mesh this with the coh grid, Rows here correspond to percentage of w
% choices, columns correspond to cash-on-hand. The columns of cash-on-hand
% is determined by ar_interp_coh_grid, because we solve the 1st stage
% problem at that coh grid.
%

armt_map('mt_interp_coh_grid_mesh_z_wage') = mt_interp_coh_grid_mesh_z_wage;
armt_map('mt_z_wage_mesh_interp_coh_grid') = mt_z_wage_mesh_interp_coh_grid;

armt_map('mt_interp_coh_grid_mesh_w_perc') = mt_interp_coh_grid_mesh_w_perc;
armt_map('mt_w_perc_mesh_interp_coh_grid') = mt_w_perc_mesh_interp_coh_grid;

armt_map('mt_interp_coh_grid_mesh_z') = mt_interp_coh_grid_mesh_z;

%% Store armt_map (3): 2nd stage reachable coh(k(w), a(w,k), z', r)
% Dimensions of Various Grids: I for level grid, M for shock grid, P for
% percent grid. These are grids for 1st stage solution
%
% # mt_coh_wkb: (I^k x I^w x M^r) by (M^z)
% # mt_z_wage_mesh_coh_wkb: (I^k x I^w x M^r) by (M^z)
%

armt_map('mt_coh_wkb') = mt_coh_wkb_full;
armt_map('mt_z_wage_mesh_coh_wkb') = mt_z_wage_mesh_coh_wkb;
% armt_map('mt_z_r_borr_mesh_coh_wkb') = mt_z_r_borr_mesh_coh_wkb;

%% Store armt_map (4): 2nd stage additional arrays
% Dimensions of Various Grids: I for level grid, M for shock grid, P for
% percent grid. These are grids for 1st stage solution
%
% # mt_k: (I^w) by (P^{k and b})
% # ar_a_meshk: 1 by (I^w x P^{k and b})
% # ar_k_mesha: 1 by (I^w x P^{k and b})
% # ar_aplusk_mesh: 1 by (I^w x P^{k and b})
% # it_ameshk_n: scalar
%

armt_map('mt_k') = mt_k;
armt_map('ar_a_meshk') = ar_a_meshk;
armt_map('ar_k_mesha') = ar_k_mesha;
armt_map('ar_aplusk_mesh') = ar_aplusk_mesh;
armt_map('it_ameshk_n') = length(ar_a_meshk);

%% Store armt_map (5): Shock Grids Arrays and Mesh
% Dimensions of Various Grids: I for level grid, M for shock grid, P for
% percent grid. These are grids for 1st stage solution
%
% # ar_z_r_borr: 1 by (M^r)
% # ar_z_r_borr_prob: 1 by (M^r)
% # ar_z_wage: 1 by (M^z)
% # ar_z_wage_prob: 1 by (M^z)
% # ar_z_r_borr_mesh_wage_w1r2: 1 by (M^z x M^r)
% # ar_z_wage_mesh_r_borr_w1r2: 1 by (M^z x M^r)
% # ar_z_r_borr_mesh_wage_r1w2: 1 by (M^r x M^z)
% # ar_z_wage_mesh_r_borr_r1w2: 1 by (M^r x M^z)
%

armt_map('ar_z_r_borr') = ar_z_r_borr;
armt_map('ar_z_r_borr_prob') = ar_z_r_borr_prob;

armt_map('ar_z_wage') = ar_z_wage;
armt_map('ar_z_wage_prob') = ar_z_wage_prob;

armt_map('ar_z_r_borr_mesh_wage_w1r2') = ar_z_r_borr_mesh_wage_w1r2;
armt_map('ar_z_wage_mesh_r_borr_w1r2') = ar_z_wage_mesh_r_borr_w1r2;
armt_map('ar_z_r_borr_mesh_wage_r1w2') = ar_z_r_borr_mesh_wage_r1w2;
armt_map('ar_z_wage_mesh_r_borr_r1w2') = ar_z_wage_mesh_r_borr_r1w2;

armt_map('mt_z_trans') = mt_z_trans;

%% Store Function Map
func_map = containers.Map('KeyType','char', 'ValueType','any');
func_map('f_util_log') = f_util_log;
func_map('f_util_crra') = f_util_crra;
func_map('f_util_standin') = f_util_standin;
func_map('f_util_standin_coh') = f_util_standin_coh;
func_map('f_prod') = f_prod;
func_map('f_inc') = f_inc;
func_map('f_coh') = f_coh;
func_map('f_cons') = f_cons;

%% Graph

if (bl_graph_funcgrids)

    %% Generate Limited Legends
    % 8 graph points, 2 levels of borrow rates, and 4 levels of rbr rates
    ar_it_z_r_borr = ([1 round((fl_z_r_borr_n)/2) (fl_z_r_borr_n)]);
    ar_it_z_wage = ([1 round((it_z_wage_n)/2) (it_z_wage_n)]);

    % combine by index
    mt_it_z_graph = ar_it_z_wage' + it_z_wage_n*(ar_it_z_r_borr-1);
    ar_it_z_graph = mt_it_z_graph(:)';
    ar_it_z_graph_zwage = ([1 round((it_z_wage_n)/4) 2*round((it_z_wage_n)/4) 3*round((it_z_wage_n)/4) (it_z_wage_n)]);

    % legends index final
    cl_st_legendCell = cellstr([num2str(ar_z_r_borr_mesh_wage_w1r2', 'zr=%3.2f;'), ...
                                num2str(ar_z_wage_mesh_r_borr_w1r2', 'zw=%3.2f')]);

    % legends index final full mat wage only
    cl_st_legendCell_zwage = cellstr([num2str(ar_z_wage', 'zw=%3.2f')]);


    %% Graph 1: a and k choice grid graphs
    % compare the figure here to the same figure in
    % <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html
    % ffs_akz_get_funcgrid>. there the grid points are on an even grid,
    % half of the grid points have NA. for the grid here, the grid points
    % get denser as we get closer to low w = k'+b' levels. This is what is
    % different visually about percentage points based choice grid for the
    % 2nd stage problem.

    figure('PaperPosition', [0 0 7 4]);
    hold on;

    chart = plot(mt_a, mt_k, 'blue');
    clr = jet(numel(chart));
    for m = 1:numel(chart)
        set(chart(m),'Color',clr(m,:))
    end
%     if (length(ar_w_level_full) <= 100)
        scatter(ar_a_meshk, ar_k_mesha, 3, 'filled', ...
            'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
%     end
    if (length(ar_w_level_full) <= 100)
        gf_invalid_scatter = scatter(mt_a_meshk_full(ar_bl_wkb_invalid),...
                                     mt_k_mesha_full(ar_bl_wkb_invalid),...
                20, 'O', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black');
    end

    xline(0);
    yline(0);

    title('Risky K Percentage Grids Given w=k+a (2nd Stage)')
    ylabel('Capital Choice (mt\_k)')
    xlabel({'Borrowing (<0) or Saving (>0) (mt\_a)'...
        'Each Diagonal Line a Different w=k+a level'...
        'Percentage for Risky K along Each Diagonal'})

    legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
    legendCell = cellstr(num2str(ar_w_level', 'k+a=%3.2f'));

    if (length(ar_w_level_full) <= 100)
        chart(length(chart)+1) = gf_invalid_scatter;
        legendCell{length(legendCell) + 1} = 'Invalid: COH(a,b,z)<bar(b) some z';
        legend(chart([legend2plot length(legendCell)]), legendCell([legend2plot length(legendCell)]), 'Location', 'northeast');
    else
        legend(chart([legend2plot]), legendCell([legend2plot]), 'Location', 'northeast');
    end

    grid on;

    %% Graph 2: coh by shock
    % compare the figure here to the same figure in
    % <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html
    % ffs_akz_get_funcgrid>. there the grid points are on an even grid.
    % Visually, one could see that the blue/red line segments here are
    % always the same length, but in the ffs_akz_get_funcgrid figure, they
    % are increasingly longer as we move towards the right. They are even
    % because the number of percentage points available is constant
    % regardless of w = k' + b' levels. But previously, the number of grid
    % points available is increasing as w increases since choice grid is
    % based on levels.
    %

    figure('PaperPosition', [0 0 7 4]);
    chart = plot(0:1:(size(mt_coh_wkb_full,1)-1), mt_coh_wkb_full);
    clr = jet(numel(chart));
    for m = 1:numel(chart)
        set(chart(m),'Color',clr(m,:))
    end

    % zero lines
    xline(0);
    yline(0);

    % invalid points separating lines
    yline_borrbound = yline(fl_b_bd);
    yline_borrbound.HandleVisibility = 'on';
    yline_borrbound.LineStyle = '--';
    yline_borrbound.Color = 'blue';
    yline_borrbound.LineWidth = 2.5;

    title('Cash-on-Hand given w(k+b),k,z');
    ylabel('Cash-on-Hand (mt\_coh\_wkb\_full)');
    xlabel({'Index of Cash-on-Hand Discrete Point (0:1:(size(mt\_coh\_wkb\_full,1)-1))'...
            'Super-Segment: borrow r; Sub-Segment: w=k+b; within seg increasing k'...
            'For each w and z, coh maximizing k is different'});

    cl_st_legendCell_here = cl_st_legendCell_zwage;

    cl_st_legendCell_here{length(cl_st_legendCell_here) + 1} = 'borrow-constraint';
    chart(length(chart)+1) = yline_borrbound;
    legend(chart([ar_it_z_graph_zwage length(cl_st_legendCell_here)]), ...
                  cl_st_legendCell_here([ar_it_z_graph_zwage length(cl_st_legendCell_here)]), 'Location', 'southeast');

    grid on;

    %% Graph 3: 1st State Aggregate Savings Choices by COH interpolation grids

    figure('PaperPosition', [0 0 7 4]);
    hold on;

    chart = plot(ar_interp_coh_grid, mt_w_perc_mesh_interp_coh_grid');

    clr = jet(numel(chart));
    for m = 1:numel(chart)
        set(chart(m),'Color',clr(m,:))
    end
    if (length(ar_interp_coh_grid) <= 100)
        [~, mt_interp_coh_grid_mesh_w_perc] = ndgrid(ar_w_perc, ar_interp_coh_grid);
        scatter(mt_interp_coh_grid_mesh_w_perc(:), mt_w_perc_mesh_interp_coh_grid(:), 3, 'filled', ...
            'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    end

    % invalid points separating lines
    yline_borrbound = yline(fl_w_level_min_valid);
    yline_borrbound.HandleVisibility = 'on';
    yline_borrbound.LineStyle = '--';
    yline_borrbound.Color = 'red';
    yline_borrbound.LineWidth = 2.5;

    xline0 = xline(0);
    xline0.HandleVisibility = 'off';
    yline0 = yline(0);
    yline0.HandleVisibility = 'off';

    title({'Aggregate Savings Percentage Grids (1st Stage)' ...
           'y=mt\_w\_by\_interp\_coh\_interp\_grid, and, y=ar\_interp\_coh\_grid'});
    ylabel('1st Stage Aggregate Savings Choices');
    xlabel({'Cash-on-Hand Levels (Interpolation Points)'...
            'w(coh)>min-agg-save, coh(k(w),w-k)>=bar(b)'});

    legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
    legendCell = cellstr(num2str(ar_w_perc', 'ar w perc=%3.2f'));
    legendCell{length(legendCell) + 1} = 'min-agg-save';
    chart(length(chart)+1) = yline_borrbound;
    legend(chart([legend2plot length(legendCell)]), legendCell([legend2plot length(legendCell)]), 'Location', 'northwest');

    grid on;

end


%% Graph Details, Generally do Not Run
if (bl_graph_funcgrids_detail)

    %% Graph 1: 2nd stage coh reached by k' b' choices by index

    figure('PaperPosition', [0 0 7 4]);
    ar_coh_kpzgrid_unique = unique(sort(mt_coh_wkb_full(:)));
    scatter(1:length(ar_coh_kpzgrid_unique), ar_coh_kpzgrid_unique);
    xline(0);
    yline(0);
    title('Cash-on-Hand given w(k+b),k,z');
    ylabel('Cash-on-Hand (y=ar\_coh\_kpzgrid\_unique)');
    xlabel({'Index of Cash-on-Hand Discrete Point' 'x = 1:length(ar\_coh\_kpzgrid\_unique)'});
    grid on;

    %% Graph 2: 2nd stage coh reached by k' b' choices  by coh

    figure('PaperPosition', [0 0 7 4]);
    ar_coh_kpzgrid_unique = unique(sort(mt_coh_wkb_full(:)));
    scatter(ar_coh_kpzgrid_unique, ar_coh_kpzgrid_unique, '.');
    xline(0);
    yline(0);
    title('Cash-on-Hand given w(k+b),k,z; See Clearly Sparsity Density of Grid across Z');
    ylabel('Cash-on-Hand (y = ar\_coh\_kpzgrid\_unique)');
    xlabel({'Cash-on-Hand' 'x = ar\_coh\_kpzgrid\_unique'});
    grid on;
end

%% Display

if (bl_display_funcgrids)

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('ar_z_wage');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(ar_z_wage));
    disp(ar_z_wage);

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('ar_w_level_full');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(ar_w_level_full));
    disp(ar_w_level_full);

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_w_perc_mesh_interp_coh_grid');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_w_perc_mesh_interp_coh_grid));
    disp(head(array2table(mt_w_perc_mesh_interp_coh_grid), 10));
    disp(tail(array2table(mt_w_perc_mesh_interp_coh_grid), 10));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_z_trans');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_z_trans));
    disp(head(array2table(mt_z_trans), 10));
    disp(tail(array2table(mt_z_trans), 10));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('ar_interp_coh_grid');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    summary(array2table(ar_interp_coh_grid'));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('ar_interp_c_grid');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    summary(array2table(ar_interp_c_grid'));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_interp_coh_grid_mesh_z');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_interp_coh_grid_mesh_z_wage));
    disp(head(array2table(mt_interp_coh_grid_mesh_z_wage), 10));
    disp(tail(array2table(mt_interp_coh_grid_mesh_z_wage), 10));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_a');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_a));
    disp(head(array2table(mt_a), 10));
    disp(tail(array2table(mt_a), 10));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('ar_a_meshk');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    summary(array2table(ar_a_meshk));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_k');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_k));
    disp(head(array2table(mt_k), 10));
    disp(tail(array2table(mt_k), 10));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('ar_k_mesha');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    summary(array2table(ar_k_mesha));

    param_map_keys = keys(func_map);
    param_map_vals = values(func_map);
    for i = 1:length(func_map)
        st_display = strjoin(['pos =' num2str(i) '; key =' string(param_map_keys{i}) '; val =' func2str(param_map_vals{i})]);
        disp(st_display);
    end

end

end
