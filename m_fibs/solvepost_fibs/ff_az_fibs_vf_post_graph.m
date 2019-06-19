%% Graph for For Inf Borr Save Outcomes for One Endo Asset Model
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

function ff_az_fibs_vf_post_graph(varargin)
%% FF_AZ_FIBS_VF_POST_GRAPH genereate 3 graphs
% Generates these graphs:
%
% # Stacked borrowing choices: separate subplot for different shocks; share
% y-axis scales; each chart, x-axis, current cash-on-hand levels, y-axis
% plot out results of several arrays. 1st array is formal choice, scattered
% height with unique color. 2nd array is the sum of formal + informal
% bridge. 3rd sum informal from formal+informal. 4th, formal savings. each
% use a different shape with increasing size so that later choices can show
% through. overall sum at the end should equal to aprime. 
%
% # Percentage credit market choices: similar file as above, except that
% now, rather than showing levels. Show percentages. divide each by
% aggregate asset level. Same Graph strategy. if 4 subplots, plotting
% with four levels of shocks. 
%
% # Participation regions, come up with all the possible formal and
% informal participation categories. What does the sorting map look like?
%
% Run this function directly with randomly generates matrixes for graphs
% and tables.
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param armt_map container container with states, choices and shocks
% grids that are inputs for grid based solution algorithm
%
% @param result_map container contains policy function matrix, value
% function matrix, iteration results; als coh consumption and other matrixes
%
% @example
%

%% Default

params_len = length(varargin);
bl_input_override = 0;
if (params_len == 5)
    bl_input_override = varargin{5};
end

if (bl_input_override)
    % if invoked from outside overrid fully
    [param_map, support_map, armt_map, result_map, ~] = varargin{:};
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
    params_group = values(armt_map, {'ar_a', 'ar_z'});
    [ar_a, ar_z] = params_group{:};
    params_group = values(func_map, {'f_util_standin', 'f_cons', 'f_coh'});
    [f_util_standin, f_cons, f_coh] = params_group{:};



    % Set Defaults
    mt_val = f_util_standin(ar_z, ar_a');
    mt_pol_a = zeros(size(mt_val)) + ar_a'*(cumsum(sort(ar_z))/sum(ar_z)*0.4 + 0.4);
    mt_cons = f_cons(ar_z, ar_a', mt_pol_a);
    mt_coh = f_coh(ar_z, ar_a');

    % Set Results Map
    result_map = containers.Map('KeyType','char', 'ValueType','any');
    result_map('mt_val') = mt_val;
    result_map('mt_pol_a') = mt_pol_a;
    result_map('mt_cons') = mt_cons;
    result_map('mt_coh') = mt_coh;
end

%% Parse Parameters

% param_map
params_group = values(param_map, {'it_z_n'});
[it_z_n] = params_group{:};

% support_map
params_group = values(support_map, {'bl_graph_onebyones', 'bl_graph_discrete', 'bl_graph_pol_lvl', 'bl_graph_pol_pct'});
[bl_graph_onebyones, bl_graph_discrete, bl_graph_pol_lvl, bl_graph_pol_pct] = params_group{:};
params_group = values(support_map, {'bl_img_save', 'st_img_path', 'st_img_prefix', 'st_img_name_main', 'st_img_suffix'});
[bl_img_save, st_img_path, st_img_prefix, st_img_name_main, st_img_suffix] = params_group{:};
params_group = values(support_map, {'st_title_prefix'});
[st_title_prefix] = params_group{:};

% armt_map
params_group = values(armt_map, {'ar_a', 'ar_z'});
[ar_a, ar_z] = params_group{:};

% result_map
params_group = values(result_map, {'mt_cons', 'mt_coh', 'mt_val', 'mt_pol_a',...
    'mt_pol_b_bridge', 'mt_pol_inf_borr_nobridge', 'mt_pol_for_borr', 'mt_pol_for_save'});
[mt_cons, mt_coh, mt_val, mt_pol_a, ...
    mt_pol_b_bridge, mt_pol_inf_borr_nobridge, mt_pol_for_borr, mt_pol_for_save] = params_group{:};

% How many zs to Graph
ar_it_z_graph = ([1 round((it_z_n)/4) round(2*((it_z_n)/4)) round(3*((it_z_n)/4)) (it_z_n)]);

%% Graph Optimal Discrete
% States: cash-on-hand, shock
% Choices: discrete optimal choices

if (bl_graph_discrete)

    if (~bl_graph_onebyones)
        figure('PaperPosition', [0 0 14 4]);
    end

    % start figure
    figure('PaperPosition', [0 0 7 4]);
    hold on;
    
    %% Graph Optimal Discrete Gen Discrete Outcomes    
    % * formal only (no bridge)
    % * Informal borrow only (no bridge)
    % * formal + informal borrow (no bridge)
    % * formal save + borrow (no bridge)
    % * bridge loan areas
    
    % Generate Binary Outcomes
    mt_it_for_borr_idx = (mt_pol_for_borr ~= 0);
    mt_it_for_save_idx = (mt_pol_for_save ~= 0);
    mt_it_inf_borr_nobridge_idx = (mt_pol_inf_borr_nobridge ~= 0);
    mt_it_b_bridge_idx = (mt_pol_b_bridge ~= 0);    
    
    % Generate Multinomial Outcomes
    mt_it_for_only_nbdg = ( mt_it_for_borr_idx & ~mt_it_for_save_idx & ~mt_it_inf_borr_nobridge_idx);
    mt_it_inf_only_nbdg = (~mt_it_for_borr_idx & ~mt_it_for_save_idx &  mt_it_inf_borr_nobridge_idx);
    mt_it_frin_brr_nbdg = ( mt_it_for_borr_idx & ~mt_it_for_save_idx &  mt_it_inf_borr_nobridge_idx);
    mt_it_fr_brrsv_nbdg = ( mt_it_for_borr_idx &  mt_it_for_save_idx & ~mt_it_inf_borr_nobridge_idx);
    
    % Generate x and y arrays
    [ar_z_mw, ar_w_mz] = meshgrid(ar_z, ar_a);
        
    %% Graph Optimal Discrete Graph Discrete Outcomes
    
    % colors etc
    cl_colors = {'blue', 'red', 'black', 'yellow', 'blue'};
    cl_legend = {'For Borr', 'Inf Borr', 'For+Inf Br', 'For+Br+Save', 'Bridge Loan'};
    cl_shapes = {'s','s','s','s','s'};
    cl_csizes = {10,10,10,10,10};
    
    % graphs to run       
    ar_it_graphs_run = 1:5;
    it_graph_counter = 0;
    for it_fig = ar_it_graphs_run
        
        % count cur graph
        it_graph_counter = it_graph_counter + 1;
        
        if (it_fig == 1)
            mt_cur_use_idx = mt_it_for_only_nbdg;
        end
        if (it_fig == 2)
            mt_cur_use_idx = mt_it_inf_only_nbdg;
        end
        if (it_fig == 3)
            mt_cur_use_idx = mt_it_frin_brr_nbdg;
        end
        if (it_fig == 4)
            mt_cur_use_idx = mt_it_fr_brrsv_nbdg;
        end
        if (it_fig == 5)
            mt_cur_use_idx = mt_it_b_bridge_idx;
        end
                
        % Generate x and y
        ar_x = mt_x(mt_cur_use_idx);
        ar_y = mt_y(mt_cur_use_idx);
        
        % Color and Size etc
        ar_color = cl_colors{it_fig};
        st_shape = cl_shapes{it_fig};
        it_csize = cl_csizes{it_fig};
        
        % Figure Collect
        ls_chart(it_graph_counter) = scatter(ar_x, ar_y, it_csize, ar_color, st_shape);
        cl_legend{it_graph_counter} = cl_legend{it_fig};
    end
    
    % legend
    legend(ls_chart, cl_legend, 'Location','northeast');
    
    % labeling
    title('Borrow and Save Regions')
    ylabel('Shocks')
    xlabel({'Total Savings w=k+b'})
    grid on;

    % save file
    if (bl_img_save)
        mkdir(support_map('st_img_path'));
        st_file_name = [st_img_prefix st_img_name_main '_fibs5' st_img_suffix];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end

end

%% Graph Choice Levels

if (bl_graph_pol_lvl)

    if (~bl_graph_onebyones)
        figure('PaperPosition', [0 0 14 8]);
        ar_sub_j = 1:1:4;
    else
        ar_sub_j = [1 3 2 4];
    end

    for sub_j = ar_sub_j

        if (sub_j==1 || sub_j == 3)
            
            % asset choice
            mt_outcome = mt_pol_a;
            
        end
        if (sub_j==2 || sub_j == 4)
            
            % consumption choice
            mt_outcome = mt_cons;
            
            % for borrowing models consumption could be at cmin, and next
            % period a' choice given default is a'=0, using the consumption
            % equation, this leads to not cmin but a negative consumption
            % level. so here adjust negative consumption to 0            
            mt_outcome(mt_cons <0) = 0;
            
        end

        if (~bl_graph_onebyones)
            subplot(2,2,sub_j)
        else
            figure('PaperPosition', [0 0 7 4]);
        end
        hold on;

        clr = jet(length(ar_it_z_graph));
        i_ctr = 0;
        for i = ar_it_z_graph
            i_ctr = i_ctr + 1;
            ar_opti_curz = mt_outcome(:, i);

            if (sub_j==1 || sub_j == 2)
                % levels
                ar_a_curz_use = ar_a';
                ar_opti_curz_use = ar_opti_curz';
            elseif (sub_j==3 || sub_j == 4)
                % logs
                ar_a_curz_use = log(ar_a' - min(ar_a) + 1);
                if (sub_j==3)
                    % borrow save
                    ar_opti_curz_use = log(ar_opti_curz' - min(ar_a) + 1);
                end
                if (sub_j == 4)
                    % consumption
                    ar_opti_curz_use = log(ar_opti_curz' + 1);
                end
            end

            scatter(ar_a_curz_use, ar_opti_curz_use, 5, ...
                'MarkerEdgeColor', clr(i_ctr,:), ...
                'MarkerFaceColor', clr(i_ctr,:));
        end

        if (sub_j==1)
            st_y_label = 'Savings/Borrowing';
            st_x_label = {'Asset (a) State'...
                          'if a''(a,z)>=a for all z, dist. shifts up'...
                          'if a''(a,z)<=a for all z, dist. shifts down'};
        end
        if (sub_j==2)
            st_y_label = 'Consumption (br cmin set to 0)';
            st_x_label = 'Asset (a) State';
        end
        if (sub_j==3)
            st_y_label = 'log(SaveBorr - min(a) + 1)';
            st_x_label = {'log(Asset State - min(a) + 1)'...
                          'if a''(a,z)>=a for all z, dist. shifts up'...
                          'if a''(a,z)<=a for all z, dist. shifts down'};
        end
        if (sub_j==4)
            st_y_label = 'log(Consumption + 1) (br cmin set to 0)';
            st_x_label = 'log(Asset State - min(a) + 1)';
        end


        grid on;
        legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));
        title([st_title_prefix st_y_label]);
        ylabel(st_y_label);
        xlabel(st_x_label);

        legend(legendCell(ar_it_z_graph), 'Location','northwest');

        hline = refline([1 0]);
        hline.Color = 'k';
        hline.LineStyle = ':';
        hline.HandleVisibility = 'off';

        if (sub_j==3 || sub_j == 4)
            xline0 = xline(log(0-min(ar_a)+1));
            xline0.HandleVisibility = 'off';
            yline0 = yline(log(0+1));
            yline0.HandleVisibility = 'off';
        else
            xline0 = xline(0);
            xline0.HandleVisibility = 'off';
            yline0 = yline(0);
            yline0.HandleVisibility = 'off';
        end
        grid on;
        grid minor;

    end

    % save file
    if (bl_img_save)
        mkdir(support_map('st_img_path'));
        st_file_name = [st_img_prefix st_img_name_main '_pol_lvl' st_img_suffix];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end

end

%% Graph Choice Percentages

if (bl_graph_pol_pct)

    if(~bl_graph_onebyones)
        figure('PaperPosition', [0 0 14 8]);
        ar_sub_j = 1:1:4;
    else
        ar_sub_j = [1 3 2 4];
    end

    for sub_j = ar_sub_j

        mt_outcome = zeros(size(mt_pol_a));
        mt_it_borr_idx = (mt_pol_a < 0);

        if (ismember(sub_j, [1,2]))
            bl_log_coh = 0;
            st_sv_suffix = '_xcoh';
            st_title_suffix = ' (x=coh)';
        else
            bl_log_coh = 1;
            st_sv_suffix = '_logxcoh';
            st_title_suffix = ' (x=log(coh))';
        end

        if (sub_j==1 || sub_j==3)
            mt_outcome(mt_it_borr_idx) = -mt_pol_a(mt_it_borr_idx)./min(ar_a);
            mt_outcome(~mt_it_borr_idx) = mt_pol_a(~mt_it_borr_idx)./mt_coh(~mt_it_borr_idx);
            st_y_label = 'aprime/min(ar_a) if br; aprime/cashonhand if sv';
            st_legend_loc = 'southeast';
            st_title = 'Save/Borrow % of Borrow Limit or COH';
        end
        if (sub_j==2 || sub_j==4)
            
            mt_cons_use = mt_cons;
            mt_cons_use(mt_cons <0) = 0;
            
            mt_outcome(mt_it_borr_idx) = mt_cons_use(mt_it_borr_idx)./(mt_coh(mt_it_borr_idx) - mt_pol_a(mt_it_borr_idx));
            mt_outcome(~mt_it_borr_idx) = mt_cons_use(~mt_it_borr_idx)./mt_coh(~mt_it_borr_idx);
            st_y_label = 'c/(coh-aprime) if br; c/cashonhand if sv';
            st_legend_loc = 'northeast';
            st_title = 'Consumption Choice As Fraction';
        end

        if (~bl_graph_onebyones)
            subplot(2,2,sub_j)
        else
            figure('PaperPosition', [0 0 7 4]);
        end
        hold on;

        clr = jet(length(ar_it_z_graph));
        i_ctr = 0;
        for i = ar_it_z_graph

            i_ctr = i_ctr + 1;
            ar_opti_curz = mt_outcome(:, i);

            if (bl_log_coh == 0)
                ar_x = ar_a;
            else
                ar_x = log(ar_a - min(ar_a) + 1);
            end

            scatter(ar_x, ar_opti_curz, 5, ...
                'MarkerEdgeColor', clr(i_ctr,:), ...
                'MarkerFaceColor', clr(i_ctr,:));
        end

        grid on;

        title([st_title_prefix st_title st_title_suffix])
        ylabel(st_y_label)

        legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));

        xlabel({'Asset State'})
        legend(legendCell(ar_it_z_graph), 'Location', st_legend_loc);
        %         xlim([min(ar_coh_curz)+1.5 15]);

        if (bl_log_coh == 0)
            xline0 = xline(0);
            xline0.HandleVisibility = 'off';
        else
            xline0 = xline(log(0-min(ar_a)+1));
            xline0.HandleVisibility = 'off';
        end

        yline0 = yline(0);
        yline0.HandleVisibility = 'off';
        yline0 = yline(1);
        yline0.HandleVisibility = 'off';

        % for save/borrow 100 and -100 percent
        if (sub_j==1)
            yline0 = yline(-1);
            yline0.HandleVisibility = 'off';
        end

        grid on;
        grid minor;
    end

    % save file
    if (bl_img_save)
        mkdir(support_map('st_img_path'));
        st_file_name = [st_img_prefix st_img_name_main '_pol_pct' st_img_suffix];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end

end

end
