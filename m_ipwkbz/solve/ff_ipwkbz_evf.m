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
    support_map('bl_display_evf') = true;
        
    param_map('it_ak_perc_n') = 250;
    param_map('fl_w_interp_grid_gap') = (param_map('fl_w_max')-param_map('fl_b_bd'))/param_map('it_ak_perc_n');    

    [armt_map, func_map] = ffs_ipwkbz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
    
    % Generating Defaults
    params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'ar_z'});
    [ar_a_meshk, ar_k_mesha, ar_z] = params_group{:};
    params_group = values(func_map, {'f_util_standin', 'f_coh'});
    [f_util_standin, f_coh] = params_group{:};
    mt_val = f_util_standin(ar_z, ar_a_meshk, ar_k_mesha);
    mt_coh = f_coh(ar_z, ar_a_meshk, ar_k_mesha);
       
end

%% Parse Parameters
params_group = values(armt_map, {'mt_z_trans', 'ar_z',...
    'ar_w_level', 'ar_k_mesha', 'ar_a_meshk', 'mt_k'});
[mt_z_trans, ar_z, ar_w_level, ...
    ar_k_mesha, ar_a_meshk, mt_k] = params_group{:};
params_group = values(param_map, {'it_z_n', 'fl_nan_replace', 'fl_b_bd'});
[it_z_n, fl_nan_replace, fl_b_bd] = params_group{:};
params_group = values(support_map, {'bl_graph_onebyones','bl_display_evf', 'bl_graph_evf'});
[bl_graph_onebyones, bl_display_evf, bl_graph_evf] = params_group{:};
params_group = values(support_map, {'bl_img_save', 'st_img_path', 'st_img_prefix', 'st_img_name_main', 'st_img_suffix'});
[bl_img_save, st_img_path, st_img_prefix, st_img_name_main, st_img_suffix] = params_group{:};

% append function name
st_func_name = 'ff_ipwkbz_evf';
st_img_name_main = [st_func_name st_img_name_main];

%% Integrate *E(V(coh(k',b'), z')|z, w)*
% Each column for a different state z, each value *E(V(coh,z')|z)* integrated already
% Here, each column is a current z, more to right higher EV
% dim(mt_ev_condi_z): *Q by M*
% Note that: mt_ev_condi_z = mt_val*mt_z_trans' is a mistake, that would be
% what we do in the
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbz/paramfunc/html/ffs_ipwkbz_set_functions.html
% ffs_ipwkbz_set_functions> code where we loop over current z, and for each
% current z, grab out a particular row from the mt_z_trans that corresponds
% to a current shock's transition into all future states.
%
% here, each column of mt_val corresponds to a state z, think of that as
% future state z. The input mt_val is *V(coh, z)*, we need to integrate to
% get *E(V(coh,z')|z)*.
%

mt_ev_condi_z = mt_val*mt_z_trans';
if(bl_display_evf)
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('mt_ev_condi_z: Q by M');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(size(mt_ev_condi_z));
    disp(head(array2table(mt_ev_condi_z), 20));
    disp(tail(array2table(mt_ev_condi_z), 20));
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
ar_add_grid = linspace(0, it_mt_bp_rown*(it_mt_bp_coln-1), it_mt_bp_coln);
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

mt_ev_condi_z_max_kp = reshape(ar_k_mesha(mt_ev_condi_z_max_idx), [it_mt_bp_coln, it_z_n]);
mt_ev_condi_z_max_bp = reshape(ar_a_meshk(mt_ev_condi_z_max_idx), [it_mt_bp_coln, it_z_n]);

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
        
        legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/3)  numel(chart)]);
        legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));
        legend(chart(legend2plot), legendCell(legend2plot), 'Location','southeast');
        
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
        
        ar_it_z_graph = ([1 round((it_z_n)/4) round(2*((it_z_n)/4)) round(3*((it_z_n)/4)) (it_z_n)]);
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
        
        legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));
        legendCell{length(legendCell) + 1} = 'max-agg-save';
        legend(legendCell([ar_it_z_graph length(legendCell)]), 'Location','southeast');
        
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
    [ar_z_mw, ar_w_mz] = meshgrid(ar_z, ar_w_level);
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
    
    [~, ar_w_mz] = meshgrid(ar_z, ar_w_level);
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
        legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/3)  numel(chart)]);
        legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));
        
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
