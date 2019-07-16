it_size_type = 102;
st_mat_prefix = '';
st_mat_name_main = 'res';
st_mat_suffix = ['_s' num2str(it_size_type)];
st_file_name = [st_mat_prefix st_mat_name_main st_mat_suffix];
st_mat_path = 'C:\Users\fan\CodeDynaAsset\m_fibs\m_ipwkbzr_test\mat\';
load(strcat(st_mat_path, st_file_name));

%% Controls
bl_graph_onebyones = true;

% cl_mt_outcomes_meansdperc_wthinfo{1}
% cl_mt_outcomes_meansdperc_wthinfo{2}
% cl_mt_outcomes_meansdperc_wthinfo{3}

%% Figure 1 Stats Formal and Informal
close all;
% for sub_j = ar_sub_j
%     
ar_sub_j = 1:1:3;

for sub_j = ar_sub_j
    
    if (bl_graph_onebyones)
    else
        figure('PaperPosition', [0 0 12 9]);
    end    
    
    mt_cur_data = cl_mt_outcomes_meansdperc_wthinfo{sub_j};
    it_y_sub1 = [6,7,8,9];
    it_y_sub2 = [11,12,13,14];

    if (sub_j==1)
        it_x_var_col = 3;

        st_x_label = 'mean informal r (poisson)';
        st_title = 'mean informal r';

    end
    if (sub_j==2)
        it_x_var_col = 4;
        
        st_x_label = 'formal borrow r';
        st_title = 'change formal borrow r';
        
    end
    if (sub_j==3)
        it_x_var_col = 5;
        
        st_x_label = 'formal menu gap';        
        st_title = 'change formal menu gap';
        
    end
    
    it_grp_gap = 14;
    ar_one = 6:it_grp_gap:size(mt_cur_data,1);
    ar_two = 7:it_grp_gap:size(mt_cur_data,1);
    ar_thr = 8:it_grp_gap:size(mt_cur_data,1);
    ar_fou = 9:it_grp_gap:size(mt_cur_data,1);
    
    ar_x = mt_cur_data(ar_one, it_x_var_col);
    
    it_mean_col = 6;
    mt_y_b_bridge = mt_cur_data(ar_one, it_mean_col);
    mt_y_inf_borr_nobridge = mt_cur_data(ar_two, it_mean_col);
    mt_y_for_borr = mt_cur_data(ar_thr, it_mean_col);
    mt_y_for_save = mt_cur_data(ar_fou, it_mean_col);
    
    %% Subplot 1
    if (bl_graph_onebyones)
        figure('PaperPosition', [0 0 7 4]);
    else
        subplot(1,2,1);
    end    
    hold on;
    
    st_legend_loc = 'southeast';
    
    blue = [57 106 177]./255;
    red = [204 37 41]./255;
    green = [62 150 81]./255;
    brown = [146 36 40]./255;
    cl_colors = {blue, red, green, brown};
    cl_scatter_shapes = {'s','x','o','d'};
    cl_linestyle = {'-','-','-','-'};
    it_sca_bs = 3;
    cl_scatter_csizes = {10*it_sca_bs, 10*it_sca_bs, 10*it_sca_bs, 10*it_sca_bs};
    it_line_bs = 2;
    cl_line_csizes = {1*it_line_bs, 1*it_line_bs, 1*it_line_bs, 1*it_line_bs};
    cl_legend = {'sum bridge','sum inf (not bridge)','sum for borr','sum for save'};

    it_graph_counter = 0;
    ls_chart = [];
    for it_fig = 1:4

        % Counter
        it_graph_counter = it_graph_counter + 1;

        % Color and Size etc
        it_csize = cl_scatter_csizes{it_fig};
        ar_color = cl_colors{it_fig};
        st_shape = cl_scatter_shapes{it_fig};
        st_lnsty = cl_linestyle{it_fig};
        st_lnwth = cl_line_csizes{it_fig};

        % plot scatter and include in legend
        if (it_fig == 1)
            ar_y = mt_y_b_bridge;
        elseif (it_fig == 2)
            ar_y = mt_y_inf_borr_nobridge;
        elseif (it_fig == 3)
            ar_y = mt_y_for_borr;
        elseif (it_fig == 4)
            ar_y = mt_y_for_save;
        end
        ls_chart(it_graph_counter) = scatter(ar_x', ar_y', it_csize, ar_color, st_shape);

        % plot line do not include in legend
        line = plot(ar_x, ar_y);
        line.HandleVisibility = 'off';
        line.Color = ar_color;
        line.LineStyle = st_lnsty;
        line.HandleVisibility = 'off';
        line.LineWidth = st_lnwth;

        cl_legend{it_graph_counter} = cl_legend{it_fig};
        
    end    
        
    st_y_label = 'Borrowing and Savings Volumns';    
    legend(ls_chart, cl_legend, 'Location', st_legend_loc);

    % 9. Titling etc
    grid on;
    title([st_title]);
    ylabel(st_y_label);
    xlabel(st_x_label)
    
    %% Subplot 2
    if (bl_graph_onebyones)
        figure('PaperPosition', [0 0 7 4]);
    else
        subplot(1,2,2);
    end
    hold on;
    
    st_legend_loc = 'northeast';
    
    it_grp_gap = 14;
    ar_one = 11:it_grp_gap:size(mt_cur_data,1);
    ar_two = 12:it_grp_gap:size(mt_cur_data,1);
    ar_thr = 13:it_grp_gap:size(mt_cur_data,1);
    ar_fou = 14:it_grp_gap:size(mt_cur_data,1);
    
    it_mean_col = 6;
    mt_y_inf_only_nbdg = mt_cur_data(ar_one, it_mean_col);
    mt_y_inf_frin_brr_nbdg = mt_cur_data(ar_two, it_mean_col);
    mt_y_for_fr_brrsv_nbdg = mt_cur_data(ar_thr, it_mean_col);
    mt_y_for_frmsavng_only = mt_cur_data(ar_fou, it_mean_col);
        
    blue = [57 106 177]./255;
    red = [204 37 41]./255;
    green = [62 150 81]./255;
    brown = [146 36 40]./255;
    cl_colors = {blue, red, green, brown};
    cl_scatter_shapes = {'s','x','o','d'};
    cl_linestyle = {'-','-','-','-'};
    it_sca_bs = 3;
    cl_scatter_csizes = {10*it_sca_bs, 10*it_sca_bs, 10*it_sca_bs, 10*it_sca_bs};
    it_line_bs = 2;
    cl_line_csizes = {1*it_line_bs, 1*it_line_bs, 1*it_line_bs, 1*it_line_bs};
    cl_legend = {'inf only','for + inf borr','for borr + save','save only'};

    it_graph_counter = 0;
    ls_chart = [];
    for it_fig = 1:4

        % Counter
        it_graph_counter = it_graph_counter + 1;

        % Color and Size etc
        it_csize = cl_scatter_csizes{it_fig};
        ar_color = cl_colors{it_fig};
        st_shape = cl_scatter_shapes{it_fig};
        st_lnsty = cl_linestyle{it_fig};
        st_lnwth = cl_line_csizes{it_fig};

        % plot scatter and include in legend
        if (it_fig == 1)
            ar_y = mt_y_inf_only_nbdg;
        elseif (it_fig == 2)
            ar_y = mt_y_inf_frin_brr_nbdg;
        elseif (it_fig == 3)
            ar_y = mt_y_for_fr_brrsv_nbdg;
        elseif (it_fig == 4)
            ar_y = mt_y_for_frmsavng_only;
        end
        ls_chart(it_graph_counter) = scatter(ar_x', ar_y', it_csize, ar_color, st_shape);

        % plot line do not include in legend
        line = plot(ar_x, ar_y);
        line.HandleVisibility = 'off';
        line.Color = ar_color;
        line.LineStyle = st_lnsty;
        line.HandleVisibility = 'off';
        line.LineWidth = st_lnwth;

        cl_legend{it_graph_counter} = cl_legend{it_fig};
        
    end    
    
    st_y_label = 'Borrowing and Savings Participation Shares';    
    legend(ls_chart, cl_legend, 'Location', st_legend_loc);

    % 9. Titling etc
    grid on;
    title([st_title]);
    ylabel(st_y_label);
    xlabel(st_x_label)    
    
    % snap figure in place    
    if (bl_graph_onebyones)
        snapnow;
    else
        st_img_path = 'C:\Users\fan\CodeDynaAsset\m_fibs\m_ipwkbzr_test\mat\';
        st_file_name = ['img1_s' num2str(it_size_type) '_m' num2str(sub_j)];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end
    
end