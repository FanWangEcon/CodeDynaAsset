%%
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

function [armt_map, func_map] = ffs_fibs_get_funcgrid(varargin)
%% FFS_FIBS_GET_FUNCGRID get funcs, params, states choices shocks grids
% centralized gateway for retrieving parameters, and solution grids and
% functions. Similar to
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html
% ffs_akz_get_funcgrid> function. This code deals with problems with
% savings and borrowing.
%
% The graphs below show the difference between percentage choice grid and
% level choice grid. See comments by graphs below for explanations of
% differences between the choice grids here and choice grids in the
% ffs_akz_get_funcgrid function.
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
%    [param_map, support_map] = ffs_fibs_set_default_param(it_param_set);
%    [armt_map, func_map] = ffs_fibs_get_funcgrid(param_map, support_map, bl_input_override);
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_fibs/paramfunc/ffs_fibs_set_functions.m ffs_fibs_set_functions>
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
    [param_map, support_map] = ffs_fibs_set_default_param();
    support_map('bl_graph_funcgrids') = true;
    support_map('bl_graph_funcgrids_detail') = true;
    support_map('bl_display_funcgrids') = true;

    % to be able to visually see choice grid points
    param_map('fl_b_bd') = -20; % borrow bound, = 0 if save only
    param_map('fl_default_aprime') = 0;
    param_map('bl_default') = 0; % if borrowing is default allowed

    param_map('fl_w_min') = param_map('fl_b_bd');
    param_map('it_w_perc_n') = 25;
    param_map('it_ak_perc_n') = 45;
    param_map('fl_w_interp_grid_gap') = 2;
    param_map('fl_coh_interp_grid_gap') = 2;
    param_map('fl_coh_interp2nd_grid_gap') = 5;    
    default_maps = {param_map, support_map};

    % numvarargs is the number of varagin inputted
    [default_maps{1:length(varargin)}] = varargin{:};
    param_map = [param_map; default_maps{1}];
    support_map = [support_map; default_maps{2}];
end

%% Parse Parameters

params_group = values(param_map, {'it_z_n', 'fl_z_mu', 'fl_z_rho', 'fl_z_sig'});
[it_z_n, fl_z_mu, fl_z_rho, fl_z_sig] = params_group{:};

params_group = values(param_map, {'fl_nan_replace', 'fl_b_bd', 'fl_w_min', 'fl_w_max', ...
    'it_w_perc_n', 'fl_w_interp_grid_gap', 'fl_coh_interp_grid_gap'});
[fl_nan_replace, fl_b_bd, fl_w_min, fl_w_max, ...
    it_w_perc_n, fl_w_interp_grid_gap, fl_coh_interp_grid_gap] = params_group{:};

params_group = values(param_map, {'fl_k_min', 'fl_k_max', 'it_ak_perc_n'});
[fl_k_min, fl_k_max, it_ak_perc_n] = params_group{:};

params_group = values(param_map, {'fl_crra', 'fl_c_min', 'it_c_interp_grid_gap', 'fl_coh_interp2nd_grid_gap'});
[fl_crra, fl_c_min, it_c_interp_grid_gap, fl_coh_interp2nd_grid_gap] = params_group{:};

params_group = values(param_map, {'fl_Amean', 'fl_alpha', 'fl_delta'});
[fl_Amean, fl_alpha, fl_delta] = params_group{:};

params_group = values(param_map, {'fl_r_save', 'fl_r_borr', 'fl_w'});
[fl_r_save, fl_r_borr, fl_w] = params_group{:};

params_group = values(support_map, {'bl_graph_funcgrids', 'bl_graph_funcgrids_detail', 'bl_display_funcgrids'});
[bl_graph_funcgrids, bl_graph_funcgrids_detail, bl_display_funcgrids] = params_group{:};

%% Generate Asset and Choice Grid for 2nd stage Problem
% This generate triangular choice structure. Household choose total
% aggregate savings, and within that how much to put into risky capital and
% how much to put into safe assets, in percentages. See
% <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc/html/ffs_fibs_set_default_param.html
% ffs_fibs_set_default_param> for details.

% percentage grid for 1st stage choice problem, level grid for 2nd stage
% solving optimal k given w and z.
ar_w_perc = linspace(0.001, 0.999, it_w_perc_n);
it_w_interp_n = (fl_w_max-fl_w_min)/(fl_w_interp_grid_gap);
ar_w_level_full = linspace(fl_w_min, fl_w_max, it_w_interp_n);
ar_w_level = ar_w_level_full;

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

% ar_a_meshk_full = mt_a(:);
% ar_k_mesha_full = mt_k(:);
% 
% ar_a_meshk = ar_a_meshk_full;
% ar_k_mesha = ar_k_mesha_full;

%% Duplicating mt_a and mt_k by coh_interp2nd_grid

% for borrowing region + 1 point for debt = 0
it_coh_interp2nd_n = (0-fl_b_bd)/(fl_coh_interp2nd_grid_gap);
ar_interp2nd_coh_grid = linspace(fl_b_bd, 0, it_coh_interp2nd_n);

% mesh these for interpolation inside VFI
[mt_w_level_by_interp2nd_coh_grid, mt_interp2nd_coh_grid_mesh_w_level] = ndgrid(ar_w_level, ar_interp2nd_coh_grid);

% expand the mt_a andmt_k matrixes
mt_a = repmat(mt_a, [1, it_coh_interp2nd_n]);
mt_k = repmat(mt_k, [1, it_coh_interp2nd_n]);

%% Array version of Choice Grids

ar_a_meshk_full = mt_a(:);
ar_k_mesha_full = mt_k(:);

ar_a_meshk = ar_a_meshk_full;
ar_k_mesha = ar_k_mesha_full;

%% Get Shock Grids

[~, mt_z_trans, ar_stationary, ar_z] = ffto_gen_tauchen_jhl(fl_z_mu,fl_z_rho,fl_z_sig,it_z_n);

%% Get Equations

[f_util_log, f_util_crra, f_util_standin, f_prod, f_inc, f_coh, f_cons] = ...
    ffs_fibs_set_functions(fl_crra, fl_c_min, fl_b_bd, fl_Amean, fl_alpha, fl_delta, fl_r_save, fl_r_borr, fl_w);

%% Generate Cash-on-Hand/State Matrix
% The endogenous state variable is cash-on-hand, it has it_z_n*it_a_n
% number of points, covering all reachable points when ar_a is the choice
% vector and ar_z is the shock vector. requires inputs from get Asset and
% choice grids, get shock grids, and get equations above.

mt_coh_wkb_full = f_coh(ar_z, ar_a_meshk_full, ar_k_mesha_full);

%% Check if COH is within Borrowing Bounds
% some coh levels are below borrowing bound, can not borrow enough to pay
% debt

mt_bl_coh_wkb_invalid = (mt_coh_wkb_full < fl_b_bd);

% (k,a) invalid if coh(k,a,z) < bd for any z
ar_bl_wkb_invalid = max(mt_bl_coh_wkb_invalid,[], 2);
mt_bl_wkb_invalid = reshape(ar_bl_wkb_invalid, size(mt_a));

% find the first w_level choice where some k(w) percent choices are valid?
ar_bl_w_level_invalid =  reshape(min(mt_bl_wkb_invalid, [], 1), size(mt_w_level_by_interp2nd_coh_grid));

% w choices can not be lower than fl_w_level_min_valid. If w choices are
% lower, given the current borrowing interest rate as well as the minimum
% income level in the future, and the maximum borrowing level available
% next period, and given the shock distribution, there exists some state in
% the future when the household when making this choice will be unable to
% borrow sufficiently to maintain positive consumption.
fl_w_level_min_valid = min(ar_w_level_full(~ar_bl_w_level_invalid(:,1)));

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
mt_z_mesh_coh_wkb = repmat(ar_z, [size(mt_coh_wkb,1),1]);

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
ar_interp_coh_grid = linspace(fl_min_mt_coh, fl_max_mt_coh, it_coh_interp_n);
[mt_interp_coh_grid_mesh_z, mt_z_mesh_coh_interp_grid] = ndgrid(ar_interp_coh_grid, ar_z);
mt_interp_coh_grid_mesh_w_perc = repmat(ar_interp_coh_grid, [it_w_perc_n, 1]);

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
    mt_w_by_interp_coh_interp_grid = ((ar_interp_coh_grid-fl_min_mt_coh)'*ar_w_perc)' + fl_min_mt_coh;
else
    % savings only
    mt_w_by_interp_coh_interp_grid = ((ar_interp_coh_grid)'*ar_w_perc)';
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

%% Store armt_map (1): 2nd Stage Problem Arrays and Matrixes
armt_map('ar_ak_perc') = ar_ak_perc;
armt_map('mt_k') = mt_k;
armt_map('ar_a_meshk') = ar_a_meshk;
armt_map('ar_k_mesha') = ar_k_mesha;
armt_map('it_ameshk_n') = length(ar_a_meshk);
armt_map('mt_coh_wkb') = mt_coh_wkb_full;
armt_map('mt_z_mesh_coh_wkb') = mt_z_mesh_coh_wkb;

%% Store armt_map (2): First Stage Aggregate Savings
% w = k' + b', w is aggregate Savings%
%
% # *ar_w_perc* 1st stage, percentage w choice given coh, at each coh
% level the number of choice points is the same for this problem with
% percentage grid points.
% # *ar_w_level* 2nd stage, level of w over which we solve the optimal
% percentage k' choices. Need to generate interpolant based on this so that
% we know optimal k* given ar_w_perc(coh) in the 1st stage
% # *mt_w_by_interp_coh_interp_grid* 1st stage, generate w(coh, percent),
% meaning the level of w given coh and the percentage grid of ar_w_perc.
% Mesh this with the coh grid, Rows here correspond to percentage of w
% choices, columns correspond to cash-on-hand. The columns of cash-on-hand
% is determined by ar_interp_coh_grid, because we solve the 1st stage
% problem at that coh grid.
%

armt_map('ar_w_perc') = ar_w_perc;
armt_map('ar_w_level') = ar_w_level;
armt_map('mt_w_by_interp_coh_interp_grid') = mt_w_by_interp_coh_interp_grid;
armt_map('mt_interp_coh_grid_mesh_w_perc') = mt_interp_coh_grid_mesh_w_perc;

%% Store armt_map (3): First Stage Consumption and Cash-on-Hand Grids

armt_map('ar_interp_c_grid') = ar_interp_c_grid;
armt_map('ar_interp_coh_grid') = ar_interp_coh_grid;
armt_map('mt_interp_coh_grid_mesh_z') = mt_interp_coh_grid_mesh_z;
armt_map('mt_z_mesh_coh_interp_grid') = mt_z_mesh_coh_interp_grid;

%% Store armt_map (4): Shock Grids
armt_map('mt_z_trans') = mt_z_trans;
armt_map('ar_stationary') = ar_stationary;
armt_map('ar_z') = ar_z;

%% Store Function Map
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
    mt_a_plot = mt_a(:, 1:length(ar_w_level));
    mt_k_plot = mt_k(:, 1:length(ar_w_level));
    chart = plot(mt_a_plot, mt_k_plot, 'blue');
    clr = jet(numel(chart));
    for m = 1:numel(chart)
        set(chart(m),'Color',clr(m,:))
    end
    if (length(ar_w_level_full) <= 100)
        scatter(ar_a_meshk, ar_k_mesha, 3, 'filled', ...
            'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    end
    if (length(ar_w_level_full) <= 100)
        gf_invalid_scatter = scatter(ar_a_meshk_full(ar_bl_wkb_invalid),...
                                     ar_k_mesha_full(ar_bl_wkb_invalid),...
                20, 'O', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black');
    end

    xline(0);
    yline(0);

    title('Risky K Percentage Grids Given w=k+a (2nd Stage)')
    ylabel('Capital Choice')
    xlabel({'Borrowing (<0) or Saving (>0)'...
        'Each Diagonal Line a Different w=k+a level'...
        'Percentage for Risky K along Each Diagonal'})

    legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
    legendCell = cellstr(num2str(ar_w_level', 'k+a=%3.2f'));

    chart(length(chart)+1) = gf_invalid_scatter;
    legendCell{length(legendCell) + 1} = 'Invalid: COH(a,b,z)<bar(b) some z';
    legend(chart([legend2plot length(legendCell)]), legendCell([legend2plot length(legendCell)]), 'Location', 'northeast');

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
    ylabel('Cash-on-Hand');
    xlabel({'Index of Cash-on-Hand Discrete Point'...
        'Each Segment is a w=k+b; within segment increasing k'...
        'For each w and z, coh maximizing k is different'});

    legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
    legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));

    legendCell{length(legendCell) + 1} = 'borrow-constraint';
    chart(length(chart)+1) = yline_borrbound;
    legend(chart([legend2plot length(legendCell)]), legendCell([legend2plot length(legendCell)]), 'Location', 'southeast');

    grid on;

    %% Graph 3: 1st State Aggregate Savings Choices by COH interpolation grids

    figure('PaperPosition', [0 0 7 4]);
    hold on;

    chart = plot(ar_interp_coh_grid, mt_w_by_interp_coh_interp_grid');

    clr = jet(numel(chart));
    for m = 1:numel(chart)
        set(chart(m),'Color',clr(m,:))
    end
    if (length(ar_interp_coh_grid) <= 100)
        [~, mt_interp_coh_grid_mesh_w_perc] = ndgrid(ar_w_perc, ar_interp_coh_grid);
        scatter(mt_interp_coh_grid_mesh_w_perc(:), mt_w_by_interp_coh_interp_grid(:), 3, 'filled', ...
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

    title('Aggregate Savings Percentage Grids (1st Stage)');
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
    ylabel('Cash-on-Hand');
    xlabel({'Index of Cash-on-Hand Discrete Point'});
    grid on;

    %% Graph 2: 2nd stage coh reached by k' b' choices  by coh

    figure('PaperPosition', [0 0 7 4]);
    ar_coh_kpzgrid_unique = unique(sort(mt_coh_wkb_full(:)));
    scatter(ar_coh_kpzgrid_unique, ar_coh_kpzgrid_unique, '.');
    xline(0);
    yline(0);
    title('Cash-on-Hand given w(k+b),k,z; See Clearly Sparsity Density of Grid across Z');
    ylabel('Cash-on-Hand');
    xlabel({'Cash-on-Hand'});
    grid on;
end

%% Display

if (bl_display_funcgrids)

    disp('ar_z');
    disp(size(ar_z));
    disp(ar_z);

    disp('ar_w_level');
    disp(size(ar_w_level_full));
    disp(ar_w_level_full');

    disp('mt_w_by_interp_coh_interp_grid');
    disp(size(mt_w_by_interp_coh_interp_grid));

    disp('mt_z_trans');
    disp(size(mt_z_trans));
    disp(mt_z_trans');

    disp('ar_interp_coh_grid');
    disp(size(ar_interp_coh_grid));
    summary(array2table(ar_interp_coh_grid'));

    disp('ar_interp_c_grid');
    disp(size(ar_interp_c_grid));
    summary(array2table(ar_interp_c_grid'));

    disp('mt_interp_coh_grid_mesh_z');
    disp(size(mt_interp_coh_grid_mesh_z));
    summary(array2table(mt_interp_coh_grid_mesh_z));

    disp('mt_a');
    disp(size(mt_a));
    %     summary(array2table(mt_a));
    %     disp(mt_a);

    disp('ar_a_meshk');
    disp(size(ar_a_meshk));
    summary(array2table(ar_a_meshk));
    %     disp(ar_a_meshk');

    disp('mt_k');
    disp(size(mt_k));
    %     summary(array2table(mt_k));
    %     disp(mt_k);

    disp('ar_k_mesha');
    disp(size(ar_k_mesha));
    summary(array2table(ar_k_mesha));
    %     disp(ar_k_mesha');

    disp('mt_coh');
    disp(size(mt_coh_wkb_full));
    summary(array2table(mt_coh_wkb_full));
    %     disp(mt_coh_wkb);

    disp('mt_coh_wkb_full')
    disp(mt_coh_wkb_full)

    param_map_keys = keys(func_map);
    param_map_vals = values(func_map);
    for i = 1:length(func_map)
        st_display = strjoin(['pos =' num2str(i) '; key =' string(param_map_keys{i}) '; val =' func2str(param_map_vals{i})]);
        disp(st_display);
    end

end

end
