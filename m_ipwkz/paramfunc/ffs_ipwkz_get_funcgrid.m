%%
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

function [armt_map, func_map] = ffs_ipwkz_get_funcgrid(varargin)
%% FFS_IPWKZ_GET_FUNCGRID get funcs, params, states choices shocks grids
% centralized gateway for retrieving parameters, and solution grids and
% functions. Similar to
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html
% ffs_akz_get_funcgrid> function.
%
% The graphs below show the difference between percentage choice grid and
% level choice grid. See comments by graphs below for explanations of
% differences between the choice grids here and choice grids in the
% ffs_akz_get_funcgrid function. 
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
%    [param_map, support_map] = ffs_ipwkz_set_default_param(it_param_set);
%    [armt_map, func_map] = ffs_ipwkz_get_funcgrid(param_map, support_map, bl_input_override);
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkz/paramfunc/ffs_ipwkz_set_functions.m ffs_ipwkz_set_functions>
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
    [param_map, support_map] = ffs_ipwkz_set_default_param();
    support_map('bl_graph_funcgrids') = true;
    support_map('bl_display_funcgrids') = true;
    
    % to be able to visually see choice grid points 
    param_map('fl_b_bd') = -10;
    param_map('fl_w_min') = param_map('fl_b_bd');
    param_map('it_w_perc_n') = 25;
    param_map('it_ak_perc_n') = 45;
    param_map('fl_w_interp_grid_gap') = 2; 
    param_map('fl_coh_interp_grid_gap') = 2;
    default_maps = {param_map, support_map};

    % numvarargs is the number of varagin inputted
    [default_maps{1:length(varargin)}] = varargin{:};
    param_map = [param_map; default_maps{1}];
    support_map = [support_map; default_maps{2}];
end

%% Parse Parameters

params_group = values(param_map, {'it_z_n', 'fl_z_mu', 'fl_z_rho', 'fl_z_sig'});
[it_z_n, fl_z_mu, fl_z_rho, fl_z_sig] = params_group{:};

params_group = values(param_map, {'fl_b_bd', 'fl_w_min', 'fl_w_max', ...
    'it_w_perc_n', 'fl_w_interp_grid_gap', 'fl_coh_interp_grid_gap'});
[fl_b_bd, fl_w_min, fl_w_max, ...
    it_w_perc_n, fl_w_interp_grid_gap, fl_coh_interp_grid_gap] = params_group{:};

params_group = values(param_map, {'fl_k_min', 'fl_k_max', 'it_ak_perc_n'});
[fl_k_min, fl_k_max, it_ak_perc_n] = params_group{:};

params_group = values(param_map, {'fl_crra', 'fl_c_min', 'it_c_interp_grid_gap'});
[fl_crra, fl_c_min, it_c_interp_grid_gap] = params_group{:};

params_group = values(param_map, {'fl_Amean', 'fl_alpha', 'fl_delta'});
[fl_Amean, fl_alpha, fl_delta] = params_group{:};

params_group = values(param_map, {'fl_r_save', 'fl_r_borr', 'fl_w'});
[fl_r_save, fl_r_borr, fl_w] = params_group{:};

params_group = values(support_map, {'bl_graph_funcgrids', 'bl_display_funcgrids'});
[bl_graph_funcgrids, bl_display_funcgrids] = params_group{:};

%% Generate Asset and Choice Grid for 2nd stage Problem
% This generate triangular choice structure. Household choose total
% aggregate savings, and within that how much to put into risky capital and
% how much to put into safe assets, in percentages. See
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkz/paramfunc/html/ffs_ipwkz_set_default_param.html
% ffs_ipwkz_set_default_param> for details.

% percentage grid for 1st stage choice problem, level grid for 2nd stage
% solving optimal k given w and z.
ar_w_perc = linspace(0, 1, it_w_perc_n);
it_w_interp_n = (fl_w_max-fl_w_min)/(fl_w_interp_grid_gap);
ar_w_level = linspace(fl_w_min, fl_w_max, it_w_interp_n);

% max k given w, need to consider the possibility of borrowing.
ar_k_max = ar_w_level - fl_b_bd;

% k percentage choice grid
ar_ak_perc = linspace(0, 1, it_ak_perc_n);

% 2nd stage percentage choice matrixes
% (ar_k_max') is it_w_interp_n by 1, and (ar_ak_perc) is 1 by it_ak_perc_n
% mt_k is a it_w_interp_n by it_ak_perc_n matrix of choice points of k'
% conditional on w, each column is a different w, each row for each col a
% different k' value. 
mt_k = (ar_k_max'*ar_ak_perc)';
mt_a = (ar_w_level - mt_k);

% can not have choice that are beyond feasible bound given the percentage
% structure here. 
mt_bl_constrained = (mt_a < fl_b_bd);
if (sum(mt_bl_constrained) > 0 )
    error('at %s second stage choice points, percentage choice exceed bounds, can not happen',...
          num2str(sum(mt_bl_constrained)));
end

ar_a_meshk = mt_a(:);
ar_k_mesha = mt_k(:);

%% Get Shock Grids

[~, mt_z_trans, ar_stationary, ar_z] = ffto_gen_tauchen_jhl(fl_z_mu,fl_z_rho,fl_z_sig,it_z_n);

%% Get Equations

[f_util_log, f_util_crra, f_util_standin, f_prod, f_inc, f_coh, f_cons] = ...
    ffs_ipwkz_set_functions(fl_crra, fl_c_min, fl_Amean, fl_alpha, fl_delta, fl_r_save, fl_r_borr, fl_w);

%% Generate Cash-on-Hand/State Matrix
% The endogenous state variable is cash-on-hand, it has it_z_n*it_a_n
% number of points, covering all reachable points when ar_a is the choice
% vector and ar_z is the shock vector. requires inputs from get Asset and
% choice grids, get shock grids, and get equations above.

mt_coh_wkb = f_coh(ar_z, ar_a_meshk, ar_k_mesha);
mt_z_mesh_coh_wkb = repmat(ar_z, [size(mt_coh_wkb,1),1]);

%% Generate 1st Stage States: Interpolation Cash-on-hand Interpolation Grid
% For the iwkz problems, we solve the problem along a grid of cash-on-hand
% values, the interpolate to find v(k',b',z) at (k',b') choices. Crucially,
% we have to coh matrxies

fl_max_mt_coh = max(max(mt_coh_wkb));
fl_min_mt_coh = min(min(mt_coh_wkb));
it_coh_interp_n = (fl_max_mt_coh-fl_min_mt_coh)/(fl_coh_interp_grid_gap);
ar_interp_coh_grid = linspace(fl_min_mt_coh, fl_max_mt_coh, it_coh_interp_n);
[mt_interp_coh_grid_mesh_z, mt_z_mesh_coh_interp_grid] = ndgrid(ar_interp_coh_grid, ar_z);

%% Generate 1st Stage Choices: Interpolation Cash-on-hand Interpolation Grid
% previously, our ar_w was the first stage choice grid, the grid was the
% same for all coh levels. Now, for each coh level, there is a different
% ar_w. ar_interp_coh_grid is (1 by ar_interp_coh_grid) and ar_w_perc is (
% 1 by it_w_perc_n). Conditional on z, each choice matrix is (it_w_perc_n
% by ar_interp_coh_grid). Here we are pre-computing the choice matrix. This
% could be a large matrix if the choice grid is large. This is the matrix
% of aggregate savings choices
%

mt_w_by_interp_coh_interp_grid = ((ar_interp_coh_grid-fl_min_mt_coh)'*ar_w_perc)' + fl_min_mt_coh;

%% Generate Interpolation Consumption Grid
% We also interpolate over consumption to speed the program up. We only
% solve for u(c) at this grid for the iwkz problmes, and then interpolate
% other c values.

fl_c_max = max(max(mt_coh_wkb)) - fl_b_bd;
it_interp_c_grid_n = (fl_c_max-fl_c_min)/(it_c_interp_grid_gap);
ar_interp_c_grid = linspace(fl_c_min, fl_c_max, it_interp_c_grid_n);

%% Store

armt_map = containers.Map('KeyType','char', 'ValueType','any');

% 2nd stage choices
armt_map('ar_ak_perc') = ar_ak_perc;
armt_map('mt_k') = mt_k;
armt_map('ar_a_meshk') = ar_a_meshk;
armt_map('ar_k_mesha') = ar_k_mesha;
armt_map('it_ameshk_n') = length(ar_a_meshk);
armt_map('mt_coh_wkb') = mt_coh_wkb;
armt_map('mt_z_mesh_coh_wkb') = mt_z_mesh_coh_wkb;

% 1st stage aggregate savings
armt_map('ar_w_perc') = ar_w_perc;
armt_map('ar_w_level') = ar_w_level;
armt_map('mt_w_by_interp_coh_interp_grid') = mt_w_by_interp_coh_interp_grid;

% 1st stage consumption coh
armt_map('ar_interp_c_grid') = ar_interp_c_grid;
armt_map('ar_interp_coh_grid') = ar_interp_coh_grid;
armt_map('mt_interp_coh_grid_mesh_z') = mt_interp_coh_grid_mesh_z;
armt_map('mt_z_mesh_coh_interp_grid') = mt_z_mesh_coh_interp_grid;

% Shocks
armt_map('mt_z_trans') = mt_z_trans;
armt_map('ar_stationary') = ar_stationary;
armt_map('ar_z') = ar_z;

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

    chart = plot(mt_a, mt_k, 'blue');

    clr = jet(numel(chart));
    for m = 1:numel(chart)
       set(chart(m),'Color',clr(m,:))
    end
    if (length(ar_w_level) <= 100)
        scatter(ar_a_meshk, ar_k_mesha, 3, 'filled', ...
                'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    end
    xline(0);
    yline(0);

    title('Choice Grids Conditional on k+a=w')
    ylabel('Capital Choice')
    xlabel({'Borrowing or Saving'})
    legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
    legendCell = cellstr(num2str(ar_w_level', 'k+a=%3.2f'));
    legend(chart(legend2plot), legendCell(legend2plot), 'Location','northeast');

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
    chart = plot(mt_coh_wkb);
    clr = jet(numel(chart));
    for m = 1:numel(chart)
       set(chart(m),'Color',clr(m,:))
    end
    xline(0);
    yline(0);

    title('Cash-on-Hand given w(k+b),k,z');
    ylabel('Cash-on-Hand');
    xlabel({'Index of Cash-on-Hand Discrete Point'...
        'Each Segment is a w=k+b; within segment increasing k'...
        'For each w and z, coh maximizing k is different'});

    legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
    legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));
    legend(chart(legend2plot), legendCell(legend2plot), 'Location','southeast');

    grid on;
    

    %% Graph 3: 2nd stage coh reached by k' b' choices by index
    
    figure('PaperPosition', [0 0 7 4]);
    ar_coh_kpzgrid_unique = unique(sort(mt_coh_wkb(:)));
    scatter(1:length(ar_coh_kpzgrid_unique), ar_coh_kpzgrid_unique);
    xline(0);
    yline(0);    
    title('Cash-on-Hand given w(k+b),k,z');
    ylabel('Cash-on-Hand');
    xlabel({'Index of Cash-on-Hand Discrete Point'});
    grid on;

    %% Graph 4: 2nd stage coh reached by k' b' choices  by coh
    
    figure('PaperPosition', [0 0 7 4]);
    ar_coh_kpzgrid_unique = unique(sort(mt_coh_wkb(:)));
    scatter(ar_coh_kpzgrid_unique, ar_coh_kpzgrid_unique, '.');
    xline(0);
    yline(0);
    title('Cash-on-Hand given w(k+b),k,z; See Clearly Sparsity Density of Grid across Z');
    ylabel('Cash-on-Hand');
    xlabel({'Cash-on-Hand'});
    grid on;
    
    %% Graph 5: 1st State Aggregate Savings Choices by COH interpolation grids
    
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
    xline0 = xline(0);
    xline0.HandleVisibility = 'off';
    yline0 = yline(0);
    yline0.HandleVisibility = 'off';
    
    title('Aggregate Savings Grid Given Cash-on-Hand');
    ylabel('1st Stage Aggregate Savings Choices');
    xlabel({'Cash-on-Hand Levels (Interpolation Points)'});

    legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
    legendCell = cellstr(num2str(ar_w_perc', 'ar w perc=%3.2f'));
    legend(chart(legend2plot), legendCell(legend2plot), 'Location','northwest');

    grid on;

end

%% Display

if (bl_display_funcgrids)

    disp('ar_z');
    disp(size(ar_z));
    disp(ar_z);

    disp('ar_w_level');
    disp(size(ar_w_level));
    disp(ar_w_level');
    
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
    disp(size(mt_coh_wkb));
    summary(array2table(mt_coh_wkb));
%     disp(mt_coh_wkb);

    param_map_keys = keys(func_map);
    param_map_vals = values(func_map);
    for i = 1:length(func_map)
        st_display = strjoin(['pos =' num2str(i) '; key =' string(param_map_keys{i}) '; val =' func2str(param_map_vals{i})]);
        disp(st_display);
    end

end

end
