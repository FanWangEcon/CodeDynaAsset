%% 
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository> 
% Table of Content.*

function ff_akz_vf_post_graph(varargin)
%% FF_AKZ_VF_POST_GRAPH genereate 3 graphs
% Generates three graphs:
%
% # Value Function Graph
% # Policy Function Consumption and Asset Choices, Level and log
% # Consumption and Asset as Percentages
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
%    ff_akz_vf_post_graph(param_map, support_map, armt_map, result_map);
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
    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override

    % Generate Default val and policy matrixes
    params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'ar_z', 'mt_coh_wkb'});
    [ar_a_meshk, ar_k_mesha, ar_z, mt_coh_wkb] = params_group{:};
    params_group = values(func_map, {'f_util_standin', 'f_cons', 'f_coh'});
    [f_util_standin, f_cons, f_coh] = params_group{:};    
        
    % Set Defaults
    mt_val = f_util_standin(ar_z, ar_a_meshk, ar_k_mesha);
    mt_pol_aksum = mt_coh_wkb.*(cumsum(sort(ar_z))/sum(ar_z)*0.4 + 0.4);
    mt_pol_a = mt_pol_aksum.*(0.7 - cumsum(sort(ar_z))/sum(ar_z)*0.3);
    mt_pol_k = mt_pol_aksum - mt_pol_a;
    mt_cons = f_cons(mt_coh_wkb, mt_pol_a, mt_pol_k);
    
    % Set Results Map
    result_map = containers.Map('KeyType','char', 'ValueType','any');
    result_map('mt_val') = mt_val;
    result_map('mt_pol_a') = mt_pol_a;
    result_map('mt_pol_k') = mt_pol_k;
    result_map('mt_cons') = mt_cons;
end

%% Parse Parameters

% param_map
params_group = values(param_map, {'fl_b_bd', 'it_z_n', 'fl_w_max'});
[fl_b_bd, it_z_n, fl_w_max] = params_group{:};

% support_map
params_group = values(support_map, {'bl_graph_onebyones', 'bl_graph_val', 'bl_graph_pol_lvl', 'bl_graph_pol_pct'});
[bl_graph_onebyones, bl_graph_val, bl_graph_pol_lvl, bl_graph_pol_pct] = params_group{:};
params_group = values(support_map, {'bl_img_save', 'st_img_path', 'st_img_prefix', 'st_img_name_main', 'st_img_suffix'});
[bl_img_save, st_img_path, st_img_prefix, st_img_name_main, st_img_suffix] = params_group{:};
params_group = values(support_map, {'st_title_prefix'});
[st_title_prefix] = params_group{:};

% armt_map
params_group = values(armt_map, {'ar_z', 'mt_coh_wkb'});
[ar_z, mt_coh_wkb] = params_group{:};

% result_map
params_group = values(result_map, {'mt_cons', 'mt_val', 'mt_pol_a', 'mt_pol_k'});
[mt_cons, mt_val, mt_pol_a, mt_pol_k] = params_group{:};

% How many zs to Graph
ar_it_z_graph = ([1 round((it_z_n)/4) round(2*((it_z_n)/4)) round(3*((it_z_n)/4)) (it_z_n)]);

%% Graphing Values
% when testing with random data using f_util_standin, shocks will not have
% impacts on z, we will see that lower shocks tend to have slightly lower
% coh values, but k,b,z effects on f_util_standin fully captured by coh. 

if (bl_graph_val)
    
    if(~bl_graph_onebyones)
        figure('PaperPosition', [0 0 7 4]);
    end    
    
    for sub_j=1:1:1
        
        if(sub_j==1)
            mt_outcome = mt_val;
            st_y_label = 'V(coh(a, k, z), z)';
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
            ar_x = mt_coh_wkb(:, i);
            ar_y = mt_outcome(:, i);
            scatter(ar_x, ar_y, 5, ...
                'MarkerEdgeColor', clr(i_ctr,:), ...
                'MarkerFaceColor', clr(i_ctr,:));
        end
        
        grid on;
        grid minor;        
        title([st_title_prefix 'COH, Shocks, and Value/Utility'])
        ylabel(st_y_label)
        xlabel({'Cash-on-Hand'})
        
        legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));
        xlinemax = xline(fl_w_max);
        xlinemax.Color = 'b';
        xlinemax.LineWidth = 1.5;
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
        st_file_name = [st_img_prefix st_img_name_main '_val' st_img_suffix];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end
    
end

%% Graphing Choice Levels

if (bl_graph_pol_lvl)
    
    if(~bl_graph_onebyones)
        figure('PaperPosition', [0 0 21 8]);
    end
    
    for sub_j=1:1:6
        
        if(sub_j==1 || sub_j == 4)
            mt_outcome = mt_pol_a;
        end
        if(sub_j==2 || sub_j == 5)
            mt_outcome = mt_pol_k;
        end
        if(sub_j==3 || sub_j == 6)
            mt_outcome = mt_cons;            
        end
                
        if(~bl_graph_onebyones)
            subplot(2,3,sub_j)
        else
            figure('PaperPosition', [0 0 7 4]);
        end               
        hold on;
        
        clr = jet(length(ar_it_z_graph));
        i_ctr = 0;
        for i = ar_it_z_graph
            i_ctr = i_ctr + 1;
            ar_opti_curz = mt_outcome(:, i);
                        
            if(sub_j==1 || sub_j == 2 || sub_j == 3)
                
                ar_a_curz_use = mt_coh_wkb(:,i)';
                ar_opti_curz_use = ar_opti_curz';                
                fl_w_max_line = fl_w_max;
                
            elseif(sub_j == 4 || sub_j == 5 || sub_j == 6)
                
                ar_a_curz_use = log(mt_coh_wkb(:,i)' - fl_b_bd + 1);
                fl_w_max_line = log(fl_w_max  - fl_b_bd + 1);
                
                if(sub_j == 4)
                    % borrow save
                    ar_opti_curz_use = log(ar_opti_curz' - fl_b_bd + 1);
                end
                
                if(sub_j == 5 || sub_j == 6)
                    % risky capital choice and consumption, both are >= 0
                    ar_opti_curz_use = log(ar_opti_curz' + 1);
                end
                
            end
            
            ar_x = ar_a_curz_use;
            ar_y = ar_opti_curz_use;
            
            scatter(ar_x, ar_y, 5, ...
                'MarkerEdgeColor', clr(i_ctr,:), ...
                'MarkerFaceColor', clr(i_ctr,:));
        end
        
        if(sub_j==1)
            st_y_label = 'Safe Savings/Borrowing';
            st_x_label = 'Cash-on-Hand';
        end
        if(sub_j==2)
            st_y_label = 'Riksy K investment';
            st_x_label = 'Cash-on-Hand';
        end
        if(sub_j==3)
            st_y_label = 'Consumption';
            st_x_label = 'Cash-on-Hand';
        end
        if(sub_j==4)
            st_y_label = 'log(SaveBorr - borrbound + 1)';
            st_x_label = 'log(COH - borrbound + 1)';
        end
        if(sub_j==5)
            st_y_label = 'log(Risky K + 1)';
            st_x_label = 'log(COH - borrbound + 1)';
        end
        if(sub_j==6)
            st_y_label = 'log(Consumption + 1)';
            st_x_label = 'log(COH - borrbound + 1)';
        end
        
        
        grid on;
        
        title([st_title_prefix 'COH, Shocks and Choices (Levels)']);
        ylabel(st_y_label);
        xlabel(st_x_label);
        
        legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));
        xlinemax = xline(fl_w_max_line);
        xlinemax.Color = 'b';
        xlinemax.LineWidth = 1.5;
        legendCell{length(legendCell) + 1} = 'max-agg-save';
        legend(legendCell([ar_it_z_graph length(legendCell)]), 'Location','northwest');
        
        hline = refline([1 0]);
        hline.Color = 'k';
        hline.LineStyle = ':';
        hline.HandleVisibility = 'off';
        if(sub_j==3 || sub_j == 4)
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

%% Graphing Choice Percentages

if (bl_graph_pol_pct)
    
    if(~bl_graph_onebyones)
        figure('PaperPosition', [0 0 14 8]);
    end    
    
    for sub_j=1:1:4
        
        mt_outcome = zeros(size(mt_pol_a));
        mt_it_borr_idx = (mt_pol_a < 0);
        
        if(sub_j==1)
            mt_outcome(mt_it_borr_idx) = -mt_pol_a(mt_it_borr_idx)./fl_b_bd;
            mt_outcome(~mt_it_borr_idx) = mt_pol_a(~mt_it_borr_idx)./mt_coh_wkb(~mt_it_borr_idx);
            st_y_label = 'aprime/borrbound if br; aprime/cashonhand if sv';
            st_legend_loc = 'southeast';
            st_title = 'Safe Savings As Fraction';
        end
        if(sub_j==2)
            mt_outcome(mt_it_borr_idx) = mt_pol_k(mt_it_borr_idx)./(mt_coh_wkb(mt_it_borr_idx) + mt_pol_a(mt_it_borr_idx));
            mt_outcome(~mt_it_borr_idx) = mt_pol_k(~mt_it_borr_idx)./mt_coh_wkb(~mt_it_borr_idx);
            st_y_label = 'kprime/(coh-aprime) if br; k/cashonhand if sv';            
            st_legend_loc = 'southeast';
            st_title = 'Risky Investment as Fraction';
        end        
        if(sub_j==3)
            %             If borrowing, how much is what is borrowing going to K?
            %             If saving, how much is what is total savings in K?
            mt_outcome(mt_it_borr_idx) = mt_pol_a(mt_it_borr_idx)./mt_pol_k(mt_it_borr_idx);
            mt_outcome(~mt_it_borr_idx) = mt_pol_a(~mt_it_borr_idx)./mt_pol_k(~mt_it_borr_idx);
            st_y_label = 'aprime/kprime';
            st_legend_loc = 'northwest';
        end
        if(sub_j==4)
            mt_outcome(mt_it_borr_idx) = mt_cons(mt_it_borr_idx)./(mt_coh_wkb(mt_it_borr_idx) + mt_pol_a(mt_it_borr_idx));
            mt_outcome(~mt_it_borr_idx) = mt_cons(~mt_it_borr_idx)./mt_coh_wkb(~mt_it_borr_idx);
            st_y_label = 'c/(coh-aprime) if br; c/cashonhand if sv';
            st_legend_loc = 'northeast';
            st_title = 'Consumption Choice As Fraction';
        end
        
        if(~bl_graph_onebyones)
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
            
            ar_x = mt_coh_wkb(:,i);
            ar_y = ar_opti_curz;
            
            scatter(ar_x, ar_y, 5, ...
                    'MarkerEdgeColor', clr(i_ctr,:), ...
                    'MarkerFaceColor', clr(i_ctr,:));
        end
        grid on;
        
        
        title([st_title_prefix st_title])
        ylabel(st_y_label)
        xlabel({'Cash-on-Hand'})
        
        legendCell = cellstr(num2str(ar_z', 'shock=%3.2f'));
        xlinemax = xline(fl_w_max);
        xlinemax.Color = 'b';
        xlinemax.LineWidth = 1.5;
        legendCell{length(legendCell) + 1} = 'max-agg-save';
        legend(legendCell([ar_it_z_graph length(legendCell)]), 'Location', st_legend_loc);
        
        %         xlim([min(ar_coh_curz)+1.5 15]);
        
        xline0 = xline(0);
        xline0.HandleVisibility = 'off';
        yline0 = yline(0);
        yline0.HandleVisibility = 'off';
        yline0 = yline(1);
        yline0.HandleVisibility = 'off';
        if(sub_j==1)
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