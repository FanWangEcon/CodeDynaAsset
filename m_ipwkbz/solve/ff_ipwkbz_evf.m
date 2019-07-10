%% 2nd Stage Optimization for Risky + Safe Asset (Save + Borr) Interpolated-Percentage
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [mt_ev_condi_z_max, mt_ev_condi_z_max_idx, mt_ev_condi_z_max_kp, mt_ev_condi_z_max_bp] = ff_ipwkbz_evf(varargin)
%% FF_IPWKBZ_EVF solves the k' vs b' problem given aggregate savings
% This function follows the structure set up here:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_wkz_evf.html
% ff_wkz_evf> but now we solve the second stage with percentage choice grid
%
% We solve along a vector of w_n vector, that is an interpolation vector,
% not a vector of actual w choices picked in the first stage. k' choices
% are in terms of percentages. Compared to ff_wkz_evf where we only had an
% upper triangle of choices, now we have a full matrix of percentage
% choices.
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
% @return mt_ev_condi_z_max matrix choice_w_n by shock_n
% max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w)) conditional on z and w, at the
% optimal k' choice (w=k'+b') what is the expected utility? This is the
% value result from the 2nd stage problem. Note the result integrates over
% z'.
%
% @return mt_ev_condi_z_max_idx matrix choice_w_n by shock_n this is the
% argmax from max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w)). Given the vector of k'
% choices, which index maximized conditional on z and w integrating over
% z'/
%
% @return mt_ev_condi_z_max_kp matrix choice_w_level_n by shock_n the k'
% choice at max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w))
%
% @return mt_ev_condi_z_max_bp matrix choice_w_n by shock_n the b'=w-k'
% choice at max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w))
%
% @example
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkbz/paramfunc/ffs_ipwkbz_set_default_param.m ffs_ipwkbz_set_default_param>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkbz/paramfunc/ffs_ipwkbz_get_funcgrid.m ffs_ipwkbz_get_funcgrid>
%

%% Default

params_len = length(varargin);
bl_input_override = 0;
if (params_len == 5)
    bl_input_override = varargin{5};
end
if (bl_input_override)
    % override when called from outside
    [mt_val, param_map, support_map, armt_map, ~] = varargin{:};
else
    clear all;
    close all;
    
    % Not default parameters, but parameters that generate defaults
    it_param_set = 4;
    bl_input_override = true;
    [param_map, support_map] = ffs_ipwkbz_set_default_param(it_param_set);
    
    support_map('bl_graph_evf') = true;
    bl_display_evf = true;
    support_map('bl_display_evf') = bl_display_evf;        
    
    st_param_which = 'default';

    if (ismember(st_param_which, ['default']))

        param_map('it_ak_perc_n') = 250;
        
%         param_map('fl_z_r_borr_min') = 0.035;
%         param_map('fl_z_r_borr_max') = 0.095;
%         param_map('fl_z_r_borr_n') = 1;        
%         param_map('fl_z_r_borr_n') = 3;
%         param_map('it_z_wage_n') = 3;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
        
    elseif (strcmp(st_param_which, 'small'))

        param_map('fl_z_r_borr_n') = 2;
        param_map('it_z_wage_n') = 3;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

        param_map('fl_b_bd') = -20; % borrow bound, = 0 if save only
        param_map('fl_default_aprime') = 0;
        param_map('bl_default') = 0; % if borrowing is default allowed

        param_map('fl_w_min') = param_map('fl_b_bd');
        param_map('it_w_perc_n') = 10;
        param_map('it_ak_perc_n') = 10;

        param_map('fl_w_interp_grid_gap') = 2;
        param_map('fl_coh_interp_grid_gap') = 2;

        param_map('fl_z_r_borr_min') = 0.025;
        param_map('fl_z_r_borr_max') = 0.95;
        param_map('fl_z_r_borr_n') = 3;
        
    end
    
    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');    
    
    param_map('fl_w_interp_grid_gap') = (param_map('fl_w_max')-param_map('fl_b_bd'))/param_map('it_ak_perc_n');    

    [armt_map, func_map] = ffs_ipwkbz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
    
    % Generating Defaults
    params_group = values(param_map, {'it_z_n'});
    [it_z_n] = params_group{:};
    
    % Generating Defaults
    params_group = values(armt_map, {'mt_coh_wkb', 'ar_z_r_borr'});
    [mt_coh_wkb, ar_z_r_borr] = params_group{:};
    params_group = values(func_map, {'f_util_standin_coh'});
    [f_util_standin_coh] = params_group{:};
    
    % Note that for the testing function below, ar_z_r_borr does not need
    % to matter for testing, meaning V(coh, zw, zr_j) = V(coh, zw, zr_i).
    % With integration it matters. This is an important point, for just
    % last period debt, if no new borrowing choices are made, it does not
    % matter what new zr shocks are, just what last period rates are. But
    % once the problem is dynamic. But the object of interest here is:
    % EV(k', b', zw, zr), conditionally on the same k'/b', will zr have an
    % impact? yes it will, even just through interest rate on b'.
    mt_val = f_util_standin_coh(mt_coh_wkb(:), ar_z_r_borr);
    mt_val = reshape(mt_val, [size(mt_coh_wkb, 1), it_z_n]);
        
    % Display Parameters
    if (bl_display_evf)
        fft_container_map_display(param_map);
        fft_container_map_display(support_map);
    end
    
end

%% Parse Parameters
params_group = values(armt_map, {'ar_z_r_borr_mesh_wage_w1r2', 'ar_z_wage_mesh_r_borr_w1r2'});
[ar_z_r_borr_mesh_wage_w1r2, ar_z_wage_mesh_r_borr_w1r2] = params_group{:};

params_group = values(armt_map, {'mt_z_trans', 'ar_ak_perc', 'ar_w_level', 'ar_k_mesha', 'ar_a_meshk', 'ar_aplusk_mesh', 'mt_k'});
[mt_z_trans, ar_ak_perc, ar_w_level, ar_k_mesha, ar_a_meshk, ar_aplusk_mesh, mt_k] = params_group{:};
params_group = values(param_map, {'it_z_n', 'fl_z_r_borr_n', 'it_z_wage_n'});
[it_z_n, fl_z_r_borr_n, it_z_wage_n] = params_group{:};
params_group = values(param_map, {'fl_nan_replace', 'fl_b_bd'});
[fl_nan_replace, fl_b_bd] = params_group{:};

params_group = values(support_map, {'bl_graph_onebyones','bl_display_evf', 'bl_graph_evf'});
[bl_graph_onebyones, bl_display_evf, bl_graph_evf] = params_group{:};
params_group = values(support_map, {'bl_img_save', 'st_img_path', 'st_img_prefix', 'st_img_name_main', 'st_img_suffix'});
[bl_img_save, st_img_path, st_img_prefix, st_img_name_main, st_img_suffix] = params_group{:};

% append function name
st_func_name = 'ff_ipwkbz_evf';
st_img_name_main = [st_func_name st_img_name_main];

%% Integrate *E(V(coh(k',b',zr'),zw',zr')|zw,zr)*
% Each column for a different state z, to integrate:
% *E(V(coh(k',b',zr'),zw',zr')|zw,zr)*. Each column is a different shock,
% from the combinations of zw and zr shocks. Each row is a different unique
% level of reacheable cash-on-hand level, which is determined by the choice
% grid for w = k' + b', k' and b', as well as the borrowing shock vector
% zr. 
%
% The issue here is, unlike the productivity shock, where only the z'
% matters tomorrow, and z matters via conditional probability of p(z'|z),
% for the interest rate shock, both r and r' matter. For the z case, z'
% impacts the cash-on-hand, and z''. For r case, r impacts cash-on-hand
% tomorrow, since interest is known at the time when the loan is taken out,
% and r' also matters because it is the rate that decision maker next
% period faces when making b'' borrowing choices. 
%
% With the structure below, the interest rate r draw that households face
% today will impact the cash-on-hand tomorrow; the r' draw will impact
% tomorrow's value function through its effect on b'' choice; r impacts r'
% through conditional probability.
%
% Note that: mt_ev_condi_z = mt_val*mt_z_trans' work if we did not have the
% interest rate shock. With the interest rate shock, we have to proceed
% differently. 
%

% 1. Number of W/B/K Choice Combinations
it_ak_perc_n = length(ar_ak_perc);
it_w_interp_n = length(ar_w_level);
it_wak_n = it_w_interp_n*it_ak_perc_n;

% 2. Initialize mt_ev_condi_z = E(V(coh(k',b',zr'),zw',zr')|zw,zr)
% rows = it_wak_n
% cols = it_z_n
mt_ev_condi_z = zeros([it_wak_n, it_z_n]);

for it_z_r_borr_ctr = 1:1:fl_z_r_borr_n
    
    % Transition Row Subset: ((M^z) by (M^z x M^r))' for one m^r
    it_mt_z_trans_row_start = it_z_wage_n*(it_z_r_borr_ctr-1) + 1;
    it_mt_z_trans_row_end = it_mt_z_trans_row_start + it_z_wage_n - 1;    
    mt_z_trans_cur_z_r_borr = mt_z_trans(it_mt_z_trans_row_start:it_mt_z_trans_row_end, :);

    % Val Segment : ((M^z) by (M^z x M^r))' for one m^r
    it_mt_val_row_start = it_wak_n*(it_z_r_borr_ctr-1) + 1;
    it_mt_val_row_end = it_mt_val_row_start + it_wak_n - 1;
    mt_val_cur_z_r_borr = mt_val(it_mt_val_row_start:it_mt_val_row_end, :);
    
    % E(V(coh(k',b',zr'),zw',zr')|zw,zr) for one zr and all zw
    mt_ev_condi_z(:, it_mt_z_trans_row_start:it_mt_z_trans_row_end) = ...
        mt_val_cur_z_r_borr*mt_z_trans_cur_z_r_borr';
    
end

if(bl_display_evf)
    
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Expected Value: mt_ev_condi_z');
    disp("EV(k', b', zw, zr) = (V(coh(k',b',zr'),zw',zr')|zw,zr)");       
    disp("rows = k'/b' combos, cols = zw/zr combos");
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    it_col_n_keep = it_z_wage_n*2;
    it_row_n_keep = it_ak_perc_n*3;
    [it_row_n, it_col_n] = size(mt_ev_condi_z);
    [ar_it_cols, ar_it_rows] = fft_row_col_subset(it_col_n, it_col_n_keep, it_row_n, it_row_n_keep);    
    cl_st_full_cols = cellstr([num2str(ar_z_r_borr_mesh_wage_w1r2', 'r%3.2f;'), ...
                               num2str(ar_z_wage_mesh_r_borr_w1r2', 'w%3.2f')]);
    cl_st_full_rows = cellstr([num2str(ar_aplusk_mesh, 'w%3.2f'), ...
                               num2str(ar_k_mesha, 'k%3.2f'),...
                               num2str(ar_a_meshk, 'a%3.2f')]);
    tb_mt_exp_val = array2table(round(mt_ev_condi_z(ar_it_rows, ar_it_cols),6));
    cl_col_names = strcat('i', num2str(ar_it_cols'), ':', cl_st_full_cols(ar_it_cols));
    cl_row_names = strcat('i', num2str(ar_it_rows'), ':', cl_st_full_rows(ar_it_rows));
    tb_mt_exp_val.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_mt_exp_val.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);
       
    disp(size(mt_ev_condi_z));
    disp(tb_mt_exp_val(1:round(it_row_n_keep/2), :));
    disp(tb_mt_exp_val((round(it_row_n_keep/2)+1):it_row_n_keep, :));

end

%% Reshape *E(V(coh,z'|z,w))* to allow for maxing
% dim(mt_ev_condi_z): *IxJ by M*

[it_mt_bp_rown, it_mt_bp_coln] = size(mt_k);
mt_ev_condi_z_full = reshape(mt_ev_condi_z, [it_mt_bp_rown, it_mt_bp_coln*it_z_n]);

%% Maximize *max_{k'}(E(V(coh(k',b'=w-k'),z'|z,w))* optimal value and index
% Maximization, find optimal k'/b' combination given z and w=k'+b'

[ar_ev_condi_z_max, ar_ev_condi_z_max_idx] = max(mt_ev_condi_z_full);
mt_ev_condi_z_max = reshape(ar_ev_condi_z_max, [it_mt_bp_coln, it_z_n]);
mt_ev_condi_z_max_idx = reshape(ar_ev_condi_z_max_idx, [it_mt_bp_coln, it_z_n]);

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
    it_col_n_keep = it_z_wage_n*2;    
    it_row_n_keep = round(it_w_interp_n);
    [it_row_n, it_col_n] = size(mt_ev_condi_z_max);
    [ar_it_cols, ar_it_rows] = fft_row_col_subset(it_col_n, it_col_n_keep, it_row_n, it_row_n_keep);    
    cl_st_full_cols = cellstr([num2str(ar_z_r_borr_mesh_wage_w1r2', 'r%3.2f;'), ...
                               num2str(ar_z_wage_mesh_r_borr_w1r2', 'w%3.2f')]);
    cl_st_full_rows = cellstr([num2str(ar_w_level', 'w%3.2f')]);
    tb_mt_ev_condi_z_max = array2table(round(mt_ev_condi_z_max(ar_it_rows, ar_it_cols), 6));
    cl_col_names = strcat('i', num2str(ar_it_cols'), ':', cl_st_full_cols(ar_it_cols));
    cl_row_names = strcat('i', num2str(ar_it_rows'), ':', cl_st_full_rows(ar_it_rows));
    tb_mt_ev_condi_z_max.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_mt_ev_condi_z_max.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);    
    disp(size(mt_ev_condi_z_max));
    disp(tb_mt_ev_condi_z_max(1:round(it_row_n_keep/2), :));
    disp(tb_mt_ev_condi_z_max((round(it_row_n_keep/2)+1):it_row_n_keep, :));

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_max_idx: I by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    tb_mt_ev_condi_z_max_idx = array2table(mt_ev_condi_z_max_idx(ar_it_rows, ar_it_cols));
    tb_mt_ev_condi_z_max_idx.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_mt_ev_condi_z_max_idx.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);    
    disp(size(mt_ev_condi_z_max_idx));
    disp(tb_mt_ev_condi_z_max_idx(1:round(it_row_n_keep/2), :));
    disp(tb_mt_ev_condi_z_max_idx((round(it_row_n_keep/2)+1):it_row_n_keep, :));
    
end

%% Reindex K' and B' Choices for each State at the Optimal *w'=k'+b'* choice
% The K' and B' Optimal Choices Associated with EV opti
% dim(mt_ev_condi_z_max_kp): *I by M*
ar_add_grid = linspace(0, it_mt_bp_rown*(it_mt_bp_coln-1), it_mt_bp_coln);
mt_ev_condi_z_max_idx = mt_ev_condi_z_max_idx + ar_add_grid';

mt_ev_condi_z_max_kp = reshape(ar_k_mesha(mt_ev_condi_z_max_idx), [it_mt_bp_coln, it_z_n]);
mt_ev_condi_z_max_bp = reshape(ar_a_meshk(mt_ev_condi_z_max_idx), [it_mt_bp_coln, it_z_n]);

if(bl_display_evf)
    
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_max_kp: I by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    tb_ev_condi_z_max_kp = array2table(mt_ev_condi_z_max_kp(ar_it_rows, ar_it_cols));
    tb_ev_condi_z_max_kp.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_ev_condi_z_max_kp.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);    
    disp(size(mt_ev_condi_z_max_kp));
    disp(tb_ev_condi_z_max_kp(1:round(it_row_n_keep/2), :));
    disp(tb_ev_condi_z_max_kp((round(it_row_n_keep/2)+1):it_row_n_keep, :));
    
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z_max_bp: I by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    tb_ev_condi_z_max_bp = array2table(mt_ev_condi_z_max_bp(ar_it_rows, ar_it_cols));
    tb_ev_condi_z_max_bp.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_ev_condi_z_max_bp.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);    
    disp(size(mt_ev_condi_z_max_kp));
    disp(tb_ev_condi_z_max_bp(1:round(it_row_n_keep/2), :));
    disp(tb_ev_condi_z_max_bp((round(it_row_n_keep/2)+1):it_row_n_keep, :));

    
end

%% Graph

if (bl_graph_evf)
    
    %% Generate Limited Legends
    % 8 graph points, 2 levels of borrow rates, and 4 levels of rbr rates
    ar_it_z_r_borr = ([1 round((fl_z_r_borr_n)/2) (fl_z_r_borr_n)]);
    ar_it_z_wage = ([1 round((it_z_wage_n)/2) (it_z_wage_n)]);

    % combine by index
    mt_it_z_graph = ar_it_z_wage' + it_z_wage_n*(ar_it_z_r_borr-1);
    ar_it_z_graph = mt_it_z_graph(:)';

    % legends index final
    cl_st_legendCell = cellstr([num2str(ar_z_r_borr_mesh_wage_w1r2', 'zr=%3.2f;'), ...
                                num2str(ar_z_wage_mesh_r_borr_w1r2', 'zw=%3.2f')]);

    
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
            ar_x = ar_w_level;
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
    [ar_z_mw, ar_w_mz] = meshgrid(ar_z_wage_mesh_r_borr_w1r2, ar_w_level);
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
    
    [~, ar_w_mz] = meshgrid(ar_z_wage_mesh_r_borr_w1r2, ar_w_level);
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
            mt_y = mt_ev_condi_z_max_kp./(ar_w_level'-fl_b_bd);
        end
        
        hold on;
        chart = plot(ar_w_level, mt_y);
        clr = jet(numel(chart));
        
        if (length(ar_w_level) <= 100)
            scatter(ar_w_mz(:), mt_y(:), 3, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        end
        
        for m = 1:numel(chart)
            set(chart(m),'Color',clr(m,:))
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

end
