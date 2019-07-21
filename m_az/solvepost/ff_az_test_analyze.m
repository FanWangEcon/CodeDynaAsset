%% Analyze Distributional Results
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [tb_outcomes, support_map] = ff_az_test_analyze(varargin)
%% FF_AZ_TEST_ANALYZE post solution simulation
% Simulate the model along various dimensions, and produce graphical
% outputs to show aggregate/distributional statistics.
%
% @ar_it_plot_sets array integer array of which statistics graphs to
% generate, see line 128. Could add additional statistics to that
% conditional list.
%
% @param bl_simu_cross boolean cross vs gridd simulation. cross: with
% (x,y), vary X fixing y, vary Y fixing x. grid: all (x,y) \in (X,Y)
%
% @param it_size_type integer: 
%  param_map support_map param_tstar_map param_desc_map
% # it_size_type = 1 is quick grid
% # it_size_type = 2 is standard grid
% # it_size_type = 3 is denser grid
%
% @param cl_st_param_keys cell string cell array container parameters that
% we are simulating over
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param param_tstar_map container map of arrays with keys that are
% parameter keys. This could be specified outside with array values to
% override defaults here. 
%
% @param param_desc_map container map of strings for each parameter key.
%
% @param tb_outcomes table table of simulation outcomes for various
% statistics for each set of parameters. Columns are statistics, rows are
% outcome variables, groups of rows are different simulations.
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_test_gen.html ff_az_test_gen>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html ff_az_vf_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html ff_az_ds_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html ffs_az_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_get_funcgrid.html ffs_az_get_funcgrid>
%
% @seealso
%
% * _SPEED_ savings only overall benchmark speed testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_speed/html/fsi_az_ds_vecsv_speed.html fsi_az_ds_vecsv_speed>
% * _PREFERENCE_ savings only preference testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_pref/html/fsi_az_ds_vecsv_pref.html fsi_az_ds_vecsv_pref>
% * _PREFERENCE_ savings only preference testing cross:
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_pref/html/fsi_az_ds_vecsv_pref_cross.html fsi_az_ds_vecsv_pref_cross>
% * _SHOCK_ savings only shock testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_shock/html/fsi_az_ds_vecsv_shock.html fsi_az_ds_vecsv_shock>
% * _SHOCK_ savings only shock testing cross:
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_shock/html/fsi_az_ds_vecsv_shock_cross.html fsi_az_ds_vecsv_shock_cross>
% * _PRICE_ savings only wage and interest rate testing cross: adjust wage and savings rate
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_price/html/fsi_az_ds_vecsv_price_cross.html fsi_az_ds_vecsv_price_cross>
%

%% Default Parameter
ar_it_plot_sets = [1,2,4,5];
it_param_set = 9;
[param_map, support_map] = ffs_az_set_default_param(it_param_set);
support_map('bl_replacefile') = false;
support_map('bl_graph_onebyones') = true;
support_map('bl_display_graph_stats') = true;

%% Array Parameters
% initiate empty map

bl_simu_cross = true; % if false, simu full grid
it_size_type = 2;
cl_st_param_keys = {'fl_crra', 'fl_beta'};

param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
it_simu_vec_len = 5;
param_tstar_map('fl_crra') = linspace(1, 5, it_simu_vec_len);
param_tstar_map('fl_beta') = linspace(0.87, 0.97, it_simu_vec_len);

param_desc_map = containers.Map('KeyType','char', 'ValueType','any');
param_desc_map('fl_crra') = {'CRRA'};
param_desc_map('fl_beta') = {'Discount'};

%% Default
default_params = {ar_it_plot_sets ...
                    bl_simu_cross it_size_type cl_st_param_keys ...
                    param_map support_map param_tstar_map param_desc_map};

%% Parse Parameters 1
% override default set above if any parameters are updated

params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};

ar_it_plot_sets = default_params{1};
bl_simu_cross = default_params{2};
it_size_type = default_params{3};
cl_st_param_keys = default_params{4};

param_map = [param_map; default_params{5}];
support_map = [support_map; default_params{6}];
param_tstar_map = [param_tstar_map; default_params{7}];
param_desc_map = [param_desc_map; default_params{8}];

%% Parse Parameters 2

params_group = values(support_map, ...
    {'bl_graph_onebyones', 'bl_display_graph_stats'});
[bl_graph_onebyones, bl_display_graph_stats] = params_group{:};

%% Cross-Simulate Model Along Parameters

[tb_outcomes, support_map, param_desc_map] = ff_az_test_gen( ...
    bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map, param_desc_map);

%% Specify Outcome Variables to Plot
it_plot_n = length(ar_it_plot_sets);
[it_plot_rows, it_plot_cols] = deal(round(it_plot_n/3), 3);

cl_ar_st_variablenames = cell([it_plot_n,1]);
cl_ar_st_legend = cell([it_plot_n,1]);
cl_ar_st_colnames = cell([it_plot_n,1]);
cl_st_title = cell([it_plot_n,1]);
cl_st_ytitle = cell([it_plot_n,1]);

it_plot_ctr = 0;
for it_plot = ar_it_plot_sets
    
    it_plot_ctr = it_plot_ctr + 1;
    
    if (it_plot == 1)
        ar_st_colnames_plot =  {'p1', 'p25', 'p50', 'mean', 'p75', 'p99'};
        ar_st_variablenames_plot =  repmat({'cl_mt_pol_c'}, [1, length(ar_st_colnames_plot)]);
        ar_st_legend_plot =  ar_st_colnames_plot;
        st_title = 'Consumption Percentiles';
        st_ytitle = 'C Distribution';
    elseif (it_plot == 2)
        ar_st_colnames_plot =  {'mean', 'sd'};
        ar_st_variablenames_plot =  repmat({'cl_mt_pol_c'}, [1, length(ar_st_colnames_plot)]);
        ar_st_legend_plot =  ar_st_colnames_plot;
        st_title = 'Consumption Mean and SD';
        st_ytitle = 'C Mean and SD';
    elseif (it_plot == 3)
        ar_st_colnames_plot =  {'p1', 'p25', 'p50', 'mean', 'p75', 'p99'};
        ar_st_variablenames_plot =  repmat({'cl_mt_pol_a'}, [1, length(ar_st_colnames_plot)]);
        ar_st_legend_plot =  ar_st_colnames_plot;
        st_title = 'Savings Percentiles';
        st_ytitle = 'A Distribution';
    elseif (it_plot == 4)
        ar_st_colnames_plot =  {'mean', 'sd'};
        ar_st_variablenames_plot =  repmat({'cl_mt_pol_a'}, [1, length(ar_st_colnames_plot)]);
        ar_st_legend_plot =  ar_st_colnames_plot;
        st_title = 'Savings Mean and SD';
        st_ytitle = 'A Mean and SD';
    elseif (it_plot == 5)
        ar_st_variablenames_plot =  {'cl_mt_coh', 'cl_mt_pol_a', 'cl_mt_pol_c'};
        ar_st_legend_plot =  {'coh=wealth', 'savings', 'consumption'};
        ar_st_colnames_plot =  repmat({'mean'}, [1, length(ar_st_variablenames_plot)]);
        st_title = 'Aggregate Outcomes (wealth, savings, consumption)';
        st_ytitle = 'Aggregate Levels';
    end
    
    cl_ar_st_variablenames{it_plot_ctr} = ar_st_variablenames_plot;
    cl_ar_st_legend{it_plot_ctr} = ar_st_legend_plot;
    cl_ar_st_colnames{it_plot_ctr} = ar_st_colnames_plot;
    cl_st_title{it_plot_ctr} = st_title;
    cl_st_ytitle{it_plot_ctr} = st_ytitle;
    
end

%% Graph Outcomes
close all;

for it_pcombi_ctr = 1:length(cl_st_param_keys)
    
    st_param_key = cl_st_param_keys{it_pcombi_ctr};
    st_param_desc = param_desc_map(st_param_key);
    
    if (~bl_graph_onebyones)
        figure('PaperPosition', [0 0 it_plot_cols*7 it_plot_rows*4]);
    else
    end
    
    % get data
    mt_cur_data = tb_outcomes(strcmp(tb_outcomes.var_param_key, st_param_key), :);
    st_x_label = st_param_desc;
    
    %% Loop over Subplots (different sets of Outcomes)
    for it_plot=1:1:it_plot_n
        
        % Get x variable and label
        cl_legend = cl_ar_st_legend{it_plot};
        ar_st_variablenames_plot = cl_ar_st_variablenames{it_plot};
        ar_st_colnames_plot = cl_ar_st_colnames{it_plot};
        
        mt_graph_data = mt_cur_data(:, [{st_param_key}, ar_st_colnames_plot]);
        
        if (bl_display_graph_stats)
            disp(mt_graph_data)
        end
        
        %% Generate Graphs
        if (bl_graph_onebyones)
            figure('PaperPosition', [0 0 7 4]);
        else
            subplot(it_plot_rows,it_plot_cols,it_plot);
        end
        
        hold on;
        
        st_legend_loc = 'southeast';
        
        blue = [57 106 177]./255;
        red = [204 37 41]./255;
        black = [83 81 84]./255;
        green = [62 150 81]./255;
        brown = [146 36 40]./255;
        purple = [107 76 154]./255;
        cl_colors = {blue, red, black, ...
            green, brown, purple};
        cl_scatter_shapes = {'s','x','o','d','p','*'};
        cl_linestyle = {'-','-','-','-','-','-'};
        it_sca_bs = 3;
        cl_scatter_csizes = {10*it_sca_bs, 10*it_sca_bs, 10*it_sca_bs, 10*it_sca_bs, 10*it_sca_bs, 10*it_sca_bs};
        it_line_bs = 2;
        cl_line_csizes = {1*it_line_bs, 1*it_line_bs, 1*it_line_bs, 1*it_line_bs, 1*it_line_bs, 1*it_line_bs};
        
        it_graph_counter = 0;
        ls_chart = [];
        for it_fig = 1:length(cl_ar_st_variablenames{it_plot})
            
            % Counter
            it_graph_counter = it_graph_counter + 1;
            
            % Color and Size etc
            it_csize = cl_scatter_csizes{it_fig};
            ar_color = cl_colors{it_fig};
            st_shape = cl_scatter_shapes{it_fig};
            st_lnsty = cl_linestyle{it_fig};
            st_lnwth = cl_line_csizes{it_fig};
            
            % Access Y Outcomes
            ar_cur_rows = strcmp(mt_cur_data.variablenames, ar_st_variablenames_plot(it_fig));
            st_cur_stat_col = ar_st_colnames_plot{it_fig};
            
            % Access X and Y Values
            mt_graph_data = mt_cur_data(ar_cur_rows, {st_param_key, st_cur_stat_col});
            % Access Y Values
            ar_y = mt_graph_data{:, st_cur_stat_col};
            % Access X Values
            ar_x = mt_graph_data{:, st_param_key};
            
            % Plot Scatter
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
        
        legend(ls_chart, cl_legend, 'Location', st_legend_loc, 'color', 'none');
        
        % 9. Titling etc
        grid on;
        title(strrep(cl_st_title{it_plot}, '_', '\_'));
        ylabel(strrep(cl_st_ytitle{it_plot}, '_', '\_'));
        xlabel(strrep(st_x_label, '_', '\_'));
        
    end
    
    % snap figure in place
    if (bl_graph_onebyones)
        snapnow;
    else
        st_img_path = 'C:\Users\fan\CodeDynaAsset\m_fibs\m_ipwkbzr_test\mat\';
        st_file_name = ['img1_s' num2str(it_size_type) '_m' num2str(it_pcombi_ctr)];
        saveas(gcf, strcat(st_img_path, st_file_name));
    end
    
end

end
