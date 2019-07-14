%% 2nd Stage Optimization (Interpolated + Percentage + Risky + Safe Asse + Save + Borr + FIBS)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [mt_ev_condi_z_max, mt_ev_condi_z_max_idx, mt_ev_condi_z_max_kp, mt_ev_condi_z_max_bp] = ff_ipwkbzr_fibs_evf(varargin)
%% FF_IPWKBZ_FIBS_EVF solves the k' vs b' problem given aggregate savings
% This file is based on
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/solve/html/ff_ipwkbzr_evf.html
% ff_ipwkbzr_evf>, see that file for more comments. Compare graphs side by
% side from this file and
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/solve/html/ff_ipwkbzr_evf.html
% ff_ipwkbzr_evf> to see visually the effect of introducing formal and
% informal choices with bridge loan.
%
% In contrast to ff_ipwkbzr_evf.m, here, we need to deal with borrowing and
% savings formal and informal. These will change how the testing matrix is
% constructed. When bridge loan is allowed, we also need to construct the
% output matrixes differently. In ff_ipwkbzr_evf.m, the assumption is that
% coh today does not matter, so to find optimal k* choice, we only need to
% know the aggregate savings level. But now, we need to know the coh level
% as well.
%
% Below two reachable coh matrixes are constructed, one for when aggregate
% savings choice w >= 0, and another for when aggregate savings <= 0. Then
% they are stacked together. And we still have the same outputs as
% ff_ipwkbzr_evf.m. The difference is that while for savings where w >=0,
% each row are w levels for the output matrixes, but for w <=0, each row is
% for w level + coh percentage combinations.
%
% @param mt_val matrix state_n I^2 by shock_n. This is the value
% matrix each row is a feasible reachable state given the choice
% vectors/matrix and each column is a shock state.
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param armt_map container container with states, choices and shocks
% grids that are inputs for grid based solution algorithm
%
% @return mt_ev_condi_z_max matrix *(choice_w_pos_n + choice_w_neg_n x
% coh_perc_n)* by *shock_n* max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w)) conditional
% on z and w, at the optimal k' choice (w=k'+b') what is the expected
% utility? This is the value result from the 2nd stage problem. Note the
% result integrates over z'.
%
% @return mt_ev_condi_z_max_idx matrix *(choice_w_pos_n + choice_w_neg_n x
% coh_perc_n)* by *(shock_n)* this is the argmax from
% max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w)). Given the vector of k' choices,
% which index maximized conditional on z and w integrating over z'.
%
% @return mt_ev_condi_z_max_kp matrix  *(choice_w_pos_n + choice_w_neg_n x
% coh_perc_n)* by *(shock_n)* the k' choice at
% max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w))
%
% @return mt_ev_condi_z_max_bp matrix  *(choice_w_pos_n + choice_w_neg_n x
% coh_perc_n)* by *(shock_n)* the b'=w-k' choice at
% max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w))
%
% @example
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_fibs/m_ipwkbzr_paramfunc/ffs_ipwkbzr_fibs_set_default_param.m ffs_ipwkbzr_fibs_set_default_param>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_fibs/m_ipwkbzr_paramfunc/ffs_ipwkbzr_fibs_get_funcgrid.m ffs_ipwkbzr_fibs_get_funcgrid>
%

%% Default
% If comparing with
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/solve/html/ff_ipwkbzr_evf.html
% ff_ipwkbzr_evf>, note that the borrowing and savings interest rates are
% the same there. Run st_param_which = 'default' to replicate identical
% result as ff_ipwkbzr_evf.m.
%

if (~isempty(varargin))
    
    % override when called from outside
    [mt_val, param_map, support_map, armt_map] = varargin{:};
    
else
    
    close all;
    % Not default parameters, but parameters that generate defaults
    it_param_set = 4;
    [param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);

    support_map('bl_graph_evf') = true;
    bl_display_evf = true;    
    support_map('bl_display_evf') = bl_display_evf;

    st_param_which = 'default';

    if (strcmp(st_param_which, 'default'))

        param_map('it_ak_perc_n') = 250;
        param_map('bl_bridge') = true;
        
    elseif (strcmp(st_param_which, 'small'))

        param_map('fl_z_r_infbr_n') = 2;
        param_map('it_z_wage_n') = 3;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');
        
        param_map('fl_b_bd') = -20; % borrow bound, = 0 if save only
        param_map('fl_default_aprime') = 0;
        param_map('bl_default') = 0; % if borrowing is default allowed

        param_map('fl_w_min') = param_map('fl_b_bd');
        param_map('it_w_perc_n') = 7;
        param_map('it_ak_perc_n') = 7;
        param_map('it_coh_bridge_perc_n') = 3;
        
        param_map('fl_w_interp_grid_gap') = 2;
        param_map('fl_coh_interp_grid_gap') = 2;
                
        param_map('fl_z_r_infbr_min') = 0.025;
        param_map('fl_z_r_infbr_max') = 0.95;
        param_map('fl_z_r_infbr_n') = 3;
        
        param_map('bl_bridge') = true;
        
    elseif (strcmp(st_param_which, 'ff_ipwkbzrr_evf'))
        
        % ff_ipwkbzrr_evf default
        param_map('fl_r_fsv') = 0.0;
        param_map('fl_r_fbr') = 1.000;
        param_map('it_ak_perc_n') = 250;
        param_map('bl_bridge') = false;
        param_map('it_coh_bridge_perc_n') = 1;
        
    elseif (strcmp(st_param_which, 'ff_ipwkbzr_evf'))
        
        % ff_ipwkbzr_evf default
        param_map('fl_r_fsv') = 0.025;
        param_map('fl_z_r_infbr_min') = 0.025;
        param_map('fl_z_r_infbr_max') = 0.025;
        param_map('fl_z_r_infbr_n') = 1;
        param_map('fl_r_fbr') = 0.025;
        param_map('it_ak_perc_n') = 250;

        param_map('bl_bridge') = false;
        param_map('it_coh_bridge_perc_n') = 1;
        
    end
    
    % Dimension Adjustments
    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');       
    param_map('fl_w_interp_grid_gap') = (param_map('fl_w_max')-param_map('fl_b_bd'))/param_map('it_ak_perc_n');    
    
    % Generate Grids
    [armt_map, func_map] = ffs_ipwkbzr_fibs_get_funcgrid(param_map, support_map);

    % Get Defaults
    params_group = values(param_map, {'it_z_n'});
    [it_z_n] = params_group{:};
    params_group = values(armt_map, {'mt_coh_wkb', 'ar_z_r_infbr'});
    [mt_coh_wkb, ar_z_r_infbr] = params_group{:};    
    params_group = values(armt_map, {'ar_ameshk_tnext_with_r', 'ar_k_mesha'});
    [ar_ameshk_tnext_with_r, ar_k_mesha] = params_group{:};
    params_group = values(func_map, {'f_util_standin_coh'});
    [f_util_standin_coh] = params_group{:};

    % mt_coh_wkb is: ((P^{k}_{a>=0} + P^{k}_{a<0} x P^{w frac bridge} ) x I^w x M^r) by (M^z) matrix
    % mt_coh_wkb(:): ((P^{k}_{a>=0} + P^{k}_{a<0} x P^{w frac bridge} ) x I^w x M^r x M^z) by 1
    % ar_z_r_infbr is: (1 x M^r)
    % mt_val: ((P^{k}_{a>=0} + P^{k}_{a<0} x P^{w frac bridge} ) x I^w x M^r x M^z) by (M^r)
    mt_val = f_util_standin_coh(mt_coh_wkb(:), ar_z_r_infbr);
    % mt_val is: (I^k x I^w x M^r) by (M^z x M^r)
    mt_val = reshape(mt_val, [size(mt_coh_wkb, 1), it_z_n]);
        
    % Display Parameters
    if (bl_display_evf)
        fft_container_map_display(param_map);
        fft_container_map_display(support_map);
    end
    
end

%% Parse Parameters

% armt_map
params_group = values(armt_map, {'ar_z_r_infbr_mesh_wage_w1r2', 'ar_z_wage_mesh_r_infbr_w1r2'});
[ar_z_r_infbr_mesh_wage_w1r2, ar_z_wage_mesh_r_infbr_w1r2] = params_group{:};
params_group = values(armt_map, {'mt_z_trans', 'ar_ak_perc', 'ar_w_level', 'ar_k_mesha', 'ar_a_meshk', 'ar_aplusk_mesh'});
[mt_z_trans, ar_ak_perc, ar_w_level, ar_k_mesha, ar_a_meshk, ar_aplusk_mesh] = params_group{:};
params_group = values(armt_map, {'ar_w_level_full'});
[ar_w_level_full] = params_group{:};

% param_map
params_group = values(param_map, {'it_z_n', 'fl_z_r_infbr_n', 'it_z_wage_n'});
[it_z_n, fl_z_r_infbr_n, it_z_wage_n] = params_group{:};
params_group = values(param_map, {'fl_nan_replace', 'fl_b_bd'});
[fl_nan_replace, fl_b_bd] = params_group{:};

% support_map
params_group = values(support_map, {'bl_graph_onebyones','bl_display_evf', 'bl_graph_evf'});
[bl_graph_onebyones, bl_display_evf, bl_graph_evf] = params_group{:};
params_group = values(support_map, {'bl_img_save', 'st_img_path', 'st_img_prefix', 'st_img_name_main', 'st_img_suffix'});
[bl_img_save, st_img_path, st_img_prefix, st_img_name_main, st_img_suffix] = params_group{:};
params_group = values(support_map, {'it_display_summmat_rowmax', 'it_display_summmat_colmax'});
[it_display_summmat_rowmax, it_display_summmat_colmax] = params_group{:};

% append function name
st_func_name = 'ff_ipwkbzr_fibs_evf';
st_img_name_main = [st_func_name st_img_name_main];

%% Integrate *E(V(coh(k',forinf(b',zr)),zw',zr')|zw,zr)*
% Start with E(V(coh(k',forinf(b',zr)),zw',zr')|zw,zr), integrate to find
% EV(k',b';zw,zr).
% 
% Note that mt_ev_condi_z rows are less by length of r rate shock times.
%
% # mt_val = (it_w_interp_n*it_ak_perc_n*length(fl_z_r_infbr_n)) by (it_z_wage_n*length(fl_z_r_infbr_n))
% # mt_ev_condi_z = (it_w_interp_n*it_ak_perc_n) by (it_z_wage_n*length(fl_z_r_infbr_n))
%

% 1. Number of W/B/K Choice Combinations
it_ak_perc_n = length(ar_ak_perc);
it_w_interp_n = length(ar_w_level_full);
it_wak_n = it_w_interp_n*it_ak_perc_n;

% 2. Initialize mt_ev_condi_z = E(V(coh(k',b',zr'),zw',zr')|zw,zr)
% rows = it_wak_n
% cols = it_z_n
mt_ev_condi_z = zeros([it_wak_n, it_z_n]);

for it_z_r_infbr_ctr = 1:1:fl_z_r_infbr_n
    
    % Transition Row Subset: ((M^z) by (M^z x M^r))' for one m^r
    it_mt_z_trans_row_start = it_z_wage_n*(it_z_r_infbr_ctr-1) + 1;
    it_mt_z_trans_row_end = it_mt_z_trans_row_start + it_z_wage_n - 1;    
    mt_z_trans_cur_z_r_infbr = mt_z_trans(it_mt_z_trans_row_start:it_mt_z_trans_row_end, :);

    % Val Segment : ((M^z) by (M^z x M^r))' for one m^r
    it_mt_val_row_start = it_wak_n*(it_z_r_infbr_ctr-1) + 1;
    it_mt_val_row_end = it_mt_val_row_start + it_wak_n - 1;
    mt_val_cur_z_r_infbr = mt_val(it_mt_val_row_start:it_mt_val_row_end, :);
    
    % E(V(coh(k',b',zr'),zw',zr')|zw,zr) for one zr and all zw
    mt_ev_condi_z(:, it_mt_z_trans_row_start:it_mt_z_trans_row_end) = ...
        mt_val_cur_z_r_infbr*mt_z_trans_cur_z_r_infbr';
    
end

%% Reshape *E(V(coh,z'|z,w))* to allow for maxing
% dim(mt_ev_condi_z): *IxJ by M*

mt_ev_condi_z_full = reshape(mt_ev_condi_z, [it_ak_perc_n, it_w_interp_n*it_z_n]);

%% Maximize *max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w))* optimal value and index
% Maximization, find optimal k'/b' combination given z and w=k'+b'

[ar_ev_condi_z_max, ar_ev_condi_z_max_idx] = max(mt_ev_condi_z_full);
mt_ev_condi_z_max = reshape(ar_ev_condi_z_max, [it_w_interp_n, it_z_n]);
mt_ev_condi_z_max_idx = reshape(ar_ev_condi_z_max_idx, [it_w_interp_n, it_z_n]);
if(bl_display_evf)

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_full: J by IxM');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_ev_condi_z_full));
%     disp(head(array2table(mt_ev_condi_z_full), 20));
%     disp(tail(array2table(mt_ev_condi_z_full), 20));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_max: I by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_ev_condi_z_max));
    disp(head(array2table(mt_ev_condi_z_max), 20));
    disp(tail(array2table(mt_ev_condi_z_max), 20));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_max_idx: I by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_ev_condi_z_max_idx));
    disp(head(array2table(mt_ev_condi_z_max_idx), 20));
    disp(tail(array2table(mt_ev_condi_z_max_idx), 20));

end

%% Reindex K' and B' Choices for each State at the Optimal *w'=k'+b'* choice
% The K' and B' Optimal Choices Associated with EV opti
% dim(mt_ev_condi_z_max_kp): *I by M*
ar_add_grid = linspace(0, it_ak_perc_n*(it_w_interp_n-1), it_w_interp_n);
mt_ev_condi_z_max_idx = mt_ev_condi_z_max_idx + ar_add_grid';

if(bl_display_evf)
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_max_idx: I by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_ev_condi_z_max_idx));
    disp(head(array2table(mt_ev_condi_z_max_idx), 20));
    disp(tail(array2table(mt_ev_condi_z_max_idx), 20));
end

mt_ev_condi_z_max_kp = reshape(ar_k_mesha(mt_ev_condi_z_max_idx), [it_w_interp_n, it_z_n]);
mt_ev_condi_z_max_bp = reshape(ar_a_meshk(mt_ev_condi_z_max_idx), [it_w_interp_n, it_z_n]);

if(bl_display_evf)
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_max_kp: I by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_ev_condi_z_max_kp));
    disp(head(array2table(mt_ev_condi_z_max_kp), 20));
    disp(tail(array2table(mt_ev_condi_z_max_kp), 20));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_max_bp: I by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_ev_condi_z_max_bp));
    disp(head(array2table(mt_ev_condi_z_max_bp), 20));
    disp(tail(array2table(mt_ev_condi_z_max_bp), 20));
end

%% Graph

if (bl_graph_evf)

    %% Generate Limited Legends
    % 8 graph points, 2 levels of borrow rates, and 4 levels of rbr rates
    ar_it_z_r_infbr = ([1 round((fl_z_r_infbr_n)/2) (fl_z_r_infbr_n)]);
    ar_it_z_wage = ([1 round((it_z_wage_n)/2) (it_z_wage_n)]);

    % combine by index
    mt_it_z_graph = ar_it_z_wage' + it_z_wage_n*(ar_it_z_r_infbr-1);
    ar_it_z_graph = mt_it_z_graph(:)';

    % legends index final
    cl_st_legendCell = cellstr([num2str(ar_z_r_infbr_mesh_wage_w1r2', 'zr=%3.2f;'), ...
                                num2str(ar_z_wage_mesh_r_infbr_w1r2', 'zw=%3.2f')]);
    
    %% Graph 1, V and EV
    if (~bl_graph_onebyones)
        figure('PaperPosition', [0 0 14 4]);
        hold on;
    end


    for subplot_j=1:1:2

        if (~bl_graph_onebyones)
            hAxis(subplot_j) = subplot(1,2,subplot_j);
        else
            figure('PaperPosition', [0 0 7 4]);
        end

        if (subplot_j==1)
            chart = plot(mt_val);
        end
        if (subplot_j==2)
            chart = plot(mt_ev_condi_z);
        end

        clr = jet(numel(chart));
        for m = 1:numel(chart)
            set(chart(m),'Color',clr(m,:))
        end

        legend(chart(ar_it_z_graph), cl_st_legendCell(ar_it_z_graph), 'Location','southeast');
        
        if (subplot_j==1)
            title('V(coh,zp); w(k+b),k,z');
        end
        if (subplot_j==2)
            title('E_z(V(coh,zp|z))');
        end

        ylabel('Next Period Value');
        xlabel({'Index of Cash-on-Hand Discrete Point'...
            'Each Segment is a w=k+b; within segment increasing k'...
            'EV and V identical if shock is fully persistent'});
        grid on;
        grid minor;
    end

    % Share y axis
    if (~bl_graph_onebyones)
        linkaxes(hAxis,'y');
    end

    % save file
    if (bl_img_save)
        mkdir(support_map('st_img_path'));
        st_file_name = [st_img_prefix st_img_name_main '_vev' st_img_suffix];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end

    %% Graph 2, max(EV)

    if(~bl_graph_onebyones)
        figure('PaperPosition', [0 0 7 4]);
    end

    for sub_j=1:1:1

        if(sub_j==1)
            mt_outcome = mt_ev_condi_z_max;
            st_y_label = 'max_{k''}(E(V(coh(k'',b''=w-k''),z''|z,w))';
        end

        if(~bl_graph_onebyones)
            subplot(1,1,sub_j)
        else
            figure('PaperPosition', [0 0 7 4]);
        end
        hold on;

        clr = jet(length(ar_it_z_graph));
        i_ctr = 0;
        for i = ar_it_z_graph
            i_ctr = i_ctr + 1;
            ar_x = ar_w_level_full;
            ar_y = mt_outcome(:, i);
            scatter(ar_x, ar_y, 5, ...
                'MarkerEdgeColor', clr(i_ctr,:), ...
                'MarkerFaceColor', clr(i_ctr,:));
        end

        grid on;
        grid minor;
        title(['2nd Stage Exp Value at Optimal K given W=K''+B'''])
        ylabel(st_y_label)
        xlabel({'Aggregate Savings'})

        legendCell_here = cl_st_legendCell;
        legendCell_here{length(legendCell_here) + 1} = 'max-agg-save';
        legend(legendCell_here([ar_it_z_graph length(legendCell_here)]), 'Location','southeast');
        
        xline0 = xline(0);
        xline0.HandleVisibility = 'off';
        yline0 = yline(0);
        yline0.HandleVisibility = 'off';

    end

    % save file
    if (bl_img_save)
        mkdir(support_map('st_img_path'));
        st_file_name = [st_img_prefix st_img_name_main '_maxev' st_img_suffix];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end

    %% Graph 3, at max(EV) optimal choice category, color regions, borrow save

    % Borrow Vs Save
    [ar_z_mw, ar_w_mz] = meshgrid(ar_z_wage_mesh_r_infbr_w1r2, ar_w_level_full);
    mt_it_borr_idx = (mt_ev_condi_z_max_bp < 0);
    mt_it_riskyhalf_idx = ((mt_ev_condi_z_max_kp./mt_ev_condi_z_max_bp) > 0.5);
    mt_it_kzero_idx = (mt_ev_condi_z_max_kp == 0);
    mt_it_isnan_idx = (isnan(mt_ev_condi_z_max_kp));

    figure('PaperPosition', [0 0 7 4]);
    % States: ar_w, ar_z
    % Choices: mt_ev_condi_z_max_kp, mt_ev_condi_z_max_bp
    hold on;
    it_sca_size = 10;
    chart_br = scatter(ar_w_mz(mt_it_borr_idx),...
        ar_z_mw(mt_it_borr_idx),...
        it_sca_size, 'blue', 'filled');
    %     legend([chart_br], {'Borrow'}, 'Location','northeast');
    chart_khalf = scatter(ar_w_mz(~mt_it_borr_idx & mt_it_riskyhalf_idx),...
        ar_z_mw(~mt_it_borr_idx & mt_it_riskyhalf_idx),...
        it_sca_size, 'black', 'filled');
    %     legend([chart_khalf], {'Save >0.5 K'}, 'Location','northeast');
    chart_sv = scatter(ar_w_mz(~mt_it_borr_idx & ~mt_it_riskyhalf_idx),...
        ar_z_mw(~mt_it_borr_idx & ~mt_it_riskyhalf_idx),...
        it_sca_size, 'red', 'filled');
    %     legend([chart_sv], {'Save <0.5 K'}, 'Location','northeast');
    chart_invalid = scatter(ar_w_mz(mt_it_kzero_idx | mt_it_isnan_idx),...
        ar_z_mw(mt_it_kzero_idx | mt_it_isnan_idx),...
        it_sca_size, 'yellow', 'filled');
    legend([chart_br, chart_khalf, chart_sv, chart_invalid], ...
        {'Borrow','Save >0.5 K','Save <0.5 K', 'k=0 or k=nan'}, 'Location','northeast');
    title('Borrow and Save Regions')
    ylabel('Shocks')
    xlabel({'Total Savings w=k+b'})
    grid on;

    % save file
    if (bl_img_save)
        mkdir(support_map('st_img_path'));
        st_file_name = [st_img_prefix st_img_name_main '_maxbrsv' st_img_suffix];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end

    %% Graph 4, Optimal K' and B' Levels
    % compare results here to results from <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/solve/html/ff_ipwkbzr_evf.html
    % ff_ipwkbzr_evf>. Several key differences:
    %
    % # Each color line is thicker here, because there is in effect another
    % state that is relevant now in the 2nd stage, which is the
    % cash-on-hand percentage, which is implemented as a percentage of the
    % w = k' + b' choice that needs to go cover bridge loan. So different
    % percentages have the same color, hence thicker lines fore each color
    % # Jump between saving and borrowing, here, the borrowing and savings
    % interest rates differ
    % # Finally, the discontinuities in choices, they occur here because of
    % the formal menu of choices, the little squiggly up and downs are due
    % to households using informal choices to complement formal choices.
    %

    [~, ar_w_mz] = meshgrid(ar_z_wage_mesh_r_infbr_w1r2, ar_w_level_full);
    for sub_j=1:1:4

        if (bl_graph_onebyones)
            figure('PaperPosition', [0 0 7 4]);
        end

        if (sub_j==1)
            if(~bl_graph_onebyones)
                figure('PaperPosition', [0 0 14 4]);
                subplot(1,2,sub_j);
            end
            mt_y = mt_ev_condi_z_max_bp;
        end
        if (sub_j==2)
            if(~bl_graph_onebyones)
                subplot(1,2,sub_j);
            end

            mt_y = mt_ev_condi_z_max_kp;
        end
        if (sub_j==3)
            if(~bl_graph_onebyones)
                figure('PaperPosition', [0 0 14 4]);
                subplot(1,2,sub_j-2);
            end
            mt_y = zeros(size(mt_ev_condi_z_max_bp));
            mt_it_borr_idx = (mt_ev_condi_z_max_bp < 0);
            mt_y(mt_it_borr_idx) = -mt_ev_condi_z_max_bp(mt_it_borr_idx)/fl_b_bd;
            mt_y(~mt_it_borr_idx) = mt_ev_condi_z_max_bp(~mt_it_borr_idx)./ar_w_mz(~mt_it_borr_idx);
        end
        if (sub_j==4)
            if(~bl_graph_onebyones)
                subplot(1,2,sub_j-2);
            end
            mt_y = mt_ev_condi_z_max_kp./(ar_w_level_full'-fl_b_bd);
        end

        hold on;
        clr = jet(it_z_n);
        for m = 1:it_z_n
            chart(m) = scatter(ar_w_level_full, mt_y(:, m), 3, ...
                'Marker', 'O', ...
                'MarkerEdgeColor', clr(m,:), 'MarkerFaceAlpha', 0.75, ...
                'MarkerFaceColor', clr(m,:), 'MarkerEdgeAlpha', 0.75);
        end

        legend2plot = fliplr(ar_it_z_graph);
        legendCell = cl_st_legendCell;

        xline0 = xline(0);
        xline0.HandleVisibility = 'off';
        yline0 = yline(0);
        yline0.HandleVisibility = 'off';
        grid on;
        if (sub_j<=2)
            hline = refline([1 0]);
            hline.Color = 'k';
            hline.LineStyle = ':';
            hline.HandleVisibility = 'off';
        end

        if (sub_j==1)
            title('B Choices of W');
            ylabel('B Choices');
            xlabel({'Total Savings w=k+b'});
            legend(chart(legend2plot), legendCell(legend2plot), 'Location','northwest');
        end
        if (sub_j==2)
            title('K Choices of W');
            ylabel('K Choices');
            xlabel({'Total Savings w=k+b'});
            legend(chart(legend2plot), legendCell(legend2plot), 'Location','northwest');
        end

        if (sub_j==3)
            title('B Fraction of Borrow Max and Save');
            ylabel('B/bar(B) if br or B/W if sv');
            xlabel({'Total Savings w=k+b'});
            %             set(gca, 'YScale', 'log');
            ylim([-1.1 1.1]);
            legend(chart(legend2plot), legendCell(legend2plot), 'Location','northwest');
        end
        if (sub_j==4)
            title('K Fraction Choices of Total K Possible');
            ylabel('K/(W-bar(b)) ');
            xlabel({'Total Savings w=k+b'});
            %             set(gca, 'YScale', 'log');
            ylim([0 1.1]);
            legend(chart(legend2plot), legendCell(legend2plot), 'Location','northeast');
        end

    end

    % save file
    if (bl_img_save)
        mkdir(support_map('st_img_path'));
        st_file_name = [st_img_prefix st_img_name_main '_wkbopti' st_img_suffix];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end

end

%% Display Various Containers

if (bl_display_evf)

    %% Display 1 support_map
    fft_container_map_display(support_map, it_display_summmat_rowmax, it_display_summmat_colmax);
        
    %% Display 2 armt_map
    fft_container_map_display(armt_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 3 param_map
    fft_container_map_display(param_map, it_display_summmat_rowmax, it_display_summmat_colmax);
    
    %% Display 4 func_map
    fft_container_map_display(func_map, it_display_summmat_rowmax, it_display_summmat_colmax);
end

end
