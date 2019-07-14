%% Generate States, Choices and Shocks Grids and Get Functions (ABZR FIBS)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [armt_map, func_map] = ffs_abzr_fibs_get_funcgrid(varargin)
%% FFS_ABZ_FIBS_GET_FUNCGRID get funcs, params, states choices shocks grids
% centralized gateway for retrieving parameters, and solution grids and
% functions
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
%    [param_map, support_map] = ffs_abzr_fibs_set_default_param(it_param_set);
%    [armt_map, func_map] = ffs_abzr_fibs_get_funcgrid(param_map, support_map, bl_input_override);
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_fibs/paramfunc/ffs_abzr_fibs_set_functions.m ffs_abzr_fibs_set_functions>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_fibs/paramfunc_fibs/ffs_for_br_block.m ffs_for_br_block>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/tools/ffto_gen_tauchen_jhl.m ffto_gen_tauchen_jhl>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/tools/fft_gen_grid_loglin.m fft_gen_grid_loglin>
%
% @seealso
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_abzr/paramfunc/html/ffs_abzr_get_funcgrid.html ffs_abzr_get_funcgrid>
%

%% Default

if (~isempty(varargin))
    
    % override when called from outside
    [param_map, support_map] = varargin{:};
    
else
    close all
    % default internal run
    it_param_set = 4;
    [param_map, support_map] = ffs_abzr_fibs_set_default_param(it_param_set);
    support_map('bl_graph_funcgrids') = true;
    support_map('bl_display_funcgrids') = true;
    default_maps = {param_map, support_map};

    % numvarargs is the number of varagin inputted
    [default_maps{1:length(varargin)}] = varargin{:};
    param_map = [param_map; default_maps{1}];
    support_map = [support_map; default_maps{2}];
end

%% Parse Parameters 1

% param_map asset grid
params_group = values(param_map, {'fl_b_bd', 'bl_default', 'fl_a_min', 'fl_a_max', 'bl_loglin', 'fl_loglin_threshold', 'it_a_n'});
[fl_b_bd, bl_default, fl_a_min, fl_a_max, bl_loglin, fl_loglin_threshold, it_a_n] = params_group{:};

% param_map preference
params_group = values(param_map, {'fl_crra', 'fl_c_min'});
[fl_crra, fl_c_min] = params_group{:};

% param_map borrowing and price
params_group = values(param_map, {'bl_b_is_principle', 'fl_r_fbr', 'fl_r_fsv', 'fl_w'});
[bl_b_is_principle, fl_r_fbr, fl_r_fsv, fl_w] = params_group{:};

%% Parse Parameters 2

% param_map shock income
params_group = values(param_map, {'it_z_wage_n', 'fl_z_wage_mu', 'fl_z_wage_rho', 'fl_z_wage_sig'});
[it_z_wage_n, fl_z_wage_mu, fl_z_wage_rho, fl_z_wage_sig] = params_group{:};

% param_map shock borrowing interest
params_group = values(param_map, {'st_z_r_infbr_drv_ele_type', 'st_z_r_infbr_drv_prb_type', 'fl_z_r_infbr_poiss_mean', ...
    'fl_z_r_infbr_max', 'fl_z_r_infbr_min', 'fl_z_r_infbr_n'});
[st_z_r_infbr_drv_ele_type, st_z_r_infbr_drv_prb_type, fl_z_r_infbr_poiss_mean, ...
    fl_z_r_infbr_max, fl_z_r_infbr_min, fl_z_r_infbr_n] = params_group{:};

% param_map formal menu
params_group = values(param_map, {'st_forbrblk_type', 'fl_forbrblk_brmost', 'fl_forbrblk_brleast', 'fl_forbrblk_gap'});
[st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap] = params_group{:};

%% Parse Parameters 3

params_group = values(support_map, {'bl_graph_funcgrids', 'bl_display_funcgrids'});
[bl_graph_funcgrids, bl_display_funcgrids] = params_group{:};

%% Get Shock: Income Shock (ar1)

[~, mt_z_wage_trans, ar_wage_stationary, ar_z_wage] = ffto_gen_tauchen_jhl(fl_z_wage_mu,fl_z_wage_rho,fl_z_wage_sig,it_z_wage_n);

%% Get Shock: Interest Rate Shock (iid)

% get borrowing grid and probabilities
param_dsv_map = containers.Map('KeyType','char', 'ValueType','any');
param_dsv_map('st_drv_ele_type') = st_z_r_infbr_drv_ele_type;
param_dsv_map('st_drv_prb_type') = st_z_r_infbr_drv_prb_type;
param_dsv_map('fl_poiss_mean') = fl_z_r_infbr_poiss_mean;
param_dsv_map('fl_max') = fl_z_r_infbr_max;
param_dsv_map('fl_min') = fl_z_r_infbr_min;
param_dsv_map('fl_n') = fl_z_r_infbr_n;
[ar_z_r_infbr, ar_z_r_inf_prob] = fft_gen_discrete_var(param_dsv_map, true);

% iid transition matrix
mt_z_r_infbr_prob_trans = repmat(ar_z_r_inf_prob, [length(ar_z_r_inf_prob), 1]);

%% Get Shock: Mesh Shocks Together

% Kronecker product to get full transition matrix for the two shocks
mt_z_trans = kron(mt_z_r_infbr_prob_trans, mt_z_wage_trans);

% mesh the shock vectors
[mt_z_wage_mesh_r_infbr, mt_z_r_infbr_mesh_wage] = ndgrid(ar_z_wage, ar_z_r_infbr);
ar_z_r_infbr_mesh_wage = mt_z_r_infbr_mesh_wage(:)';
ar_z_wage_mesh_r_infbr = mt_z_wage_mesh_r_infbr(:)';

if (bl_display_funcgrids)
    
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Borrow R Shock: ar_z_r_inf_mesh_wage');
    disp('Prod/Wage Shock: mt_z_wage_mesh_r_infbr');
    disp('show which shock is inner and which is outter');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    
    tb_two_shocks = array2table([ar_z_r_infbr_mesh_wage;...
        ar_z_wage_mesh_r_infbr]');
    cl_col_names = ["Borrow R Shock (Meshed)", "Wage R Shock (Meshed)"];
    cl_row_names = strcat('zi=', string((1:length(ar_z_r_infbr_mesh_wage))));
    tb_two_shocks.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_two_shocks.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);
    
    it_row_display = it_z_wage_n*2;
    
    disp(size(tb_two_shocks));
    disp(head(tb_two_shocks, it_row_display));
    disp(tail(tb_two_shocks, it_row_display));
    
    
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Borrow Rate Transition Matrix: mt_z_r_infbr_prob_trans');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    it_col_n_keep = 15;
    it_row_n_keep = 15;    
    [it_row_n, it_col_n] = size(mt_z_r_infbr_prob_trans);
    [ar_it_cols, ar_it_rows] = fft_row_col_subset(it_col_n, it_col_n_keep, it_row_n, it_row_n_keep);    
    cl_st_full_rowscols = cellstr([num2str(ar_z_r_infbr', 'r%3.2f')]);
    tb_z_r_infbr_prob_trans = array2table(round(mt_z_r_infbr_prob_trans(ar_it_rows, ar_it_cols), 6));
    cl_col_names = strcat('zi=', num2str(ar_it_cols'), ':', cl_st_full_rowscols(ar_it_cols));
    cl_row_names = strcat('zi=', num2str(ar_it_rows'), ':', cl_st_full_rowscols(ar_it_cols));
    tb_z_r_infbr_prob_trans.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_z_r_infbr_prob_trans.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);
       
    disp(size(tb_z_r_infbr_prob_trans));
    disp(tb_z_r_infbr_prob_trans);
    
    
    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Wage Prod Shock Transition Matrix: mt_z_r_infbr_prob_trans');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    it_col_n_keep = 15;
    it_row_n_keep = 15;    
    [it_row_n, it_col_n] = size(mt_z_wage_trans);
    [ar_it_cols, ar_it_rows] = fft_row_col_subset(it_col_n, it_col_n_keep, it_row_n, it_row_n_keep);    
    cl_st_full_rowscols = cellstr([num2str(ar_z_wage', 'w%3.2f')]);
    tb_z_wage_trans = array2table(round(mt_z_wage_trans(ar_it_rows, ar_it_cols),6));    
    cl_col_names = strcat('zi=', num2str(ar_it_cols'), ':', cl_st_full_rowscols(ar_it_cols));
    cl_row_names = strcat('zi=', num2str(ar_it_rows'), ':', cl_st_full_rowscols(ar_it_cols));
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
    cl_st_full_rowscols = cellstr([num2str(ar_z_r_infbr_mesh_wage', 'r%3.2f;'), ...
                                   num2str(ar_z_wage_mesh_r_infbr', 'w%3.2f')]);
    tb_mt_z_trans = array2table(round(mt_z_trans(ar_it_rows, ar_it_cols),6));
    cl_col_names = strcat('i', num2str(ar_it_cols'), ':', cl_st_full_rowscols(ar_it_cols));
    cl_row_names = strcat('i', num2str(ar_it_rows'), ':', cl_st_full_rowscols(ar_it_cols));
    tb_mt_z_trans.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_mt_z_trans.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);
       
    disp(size(tb_mt_z_trans));
    disp(tb_mt_z_trans);
    
end

%% Get Equations

[f_util_log, f_util_crra, f_util_standin, f_inc, f_coh, f_cons_coh_fbis, f_cons_coh_save, f_bprime] = ...
    ffs_abzr_fibs_set_functions(fl_crra, fl_c_min, fl_r_fbr, fl_r_fsv, fl_w);

%% Get Formal Borrowing Blocks

[ar_forbrblk, ar_forbrblk_r] = ...
        ffs_for_br_block_gen(fl_r_fbr, st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap);

%% Get Asset and Choice Grid
% note this requires ar_z

if (bl_loglin)
    % C:\Users\fan\M4Econ\asset\grid\ff_grid_loglin.m
    ar_a = fft_gen_grid_loglin(it_a_n, fl_a_max, fl_a_min, fl_loglin_threshold);
else
    fl_r_inf_max = max(ar_z_r_infbr);
    
    [ar_a_inf, fl_borr_yminbd_inf, fl_borr_ymaxbd_inf] = ffs_abz_gen_borrsave_grid(...
        fl_b_bd, bl_default, ar_z_wage, fl_w, ...
        bl_b_is_principle, fl_r_inf_max, fl_a_min, fl_a_max, it_a_n);

    [ar_a_for, fl_borr_yminbd_for, fl_borr_ymaxbd_for] = ffs_abz_gen_borrsave_grid(...
        fl_b_bd, bl_default, ar_z_wage, fl_w, ...
        bl_b_is_principle, fl_r_fbr, fl_a_min, fl_a_max, it_a_n);

    if (min(ar_a_for) <= min(ar_a_inf))
        ar_a = ar_a_for;
        fl_borr_yminbd = fl_borr_yminbd_for;
        fl_borr_ymaxbd = fl_borr_ymaxbd_for;
    else
        ar_a = ar_a_inf;
        fl_borr_yminbd = fl_borr_yminbd_inf;
        fl_borr_ymaxbd = fl_borr_ymaxbd_inf;
    end

end


%% Store

armt_map = containers.Map('KeyType','char', 'ValueType','any');
armt_map('ar_a') = ar_a;
armt_map('mt_z_trans') = mt_z_trans;
armt_map('ar_z_r_infbr_mesh_wage') = ar_z_r_infbr_mesh_wage;
armt_map('ar_z_wage_mesh_r_infbr') = ar_z_wage_mesh_r_infbr;
armt_map('ar_forbrblk') = ar_forbrblk;
armt_map('ar_forbrblk_r') = ar_forbrblk_r;

func_map = containers.Map('KeyType','char', 'ValueType','any');
func_map('f_util_log') = f_util_log;
func_map('f_util_crra') = f_util_crra;
func_map('f_util_standin') = f_util_standin;
func_map('f_inc') = f_inc;
func_map('f_coh') = f_coh;
func_map('f_cons_coh_fbis') = f_cons_coh_fbis;
func_map('f_cons_coh_save') = f_cons_coh_save;
func_map('f_bprime') = f_bprime;

%% Graph: A, Shocks, COH, and Defaults
% # y-axis : coh(a,z)
% # x-axis : a
% # color: z
% # overlay: coh points points where there is default.

if (bl_graph_funcgrids)

    % mesh a and and z
    [mt_a_mesh_z, mt_z_mesh_a] = ndgrid(ar_a, ar_z_wage);

    % cash-on-hand given a and z
    mt_coh = f_coh(mt_z_mesh_a, mt_a_mesh_z);

    % loop over level vs log graphs
    for sub_j=1:1:1

        figure('PaperPosition', [0 0 7 4]);

        if (sub_j == 1)
            x_mat = mt_a_mesh_z;
            y_mat = mt_coh;
            st_title = 'coh(a,z)';
            st_ylabel = 'Cash-on-hand(a, z)';
            st_xlabel = 'Asset States (Choices)';

            fl_b_bd_graph = fl_b_bd;
            fl_borr_yminbd_graph = fl_borr_yminbd;
            fl_borr_ymaxbd_graph = fl_borr_ymaxbd;
        else
            x_mat = log(mt_a_mesh_z - min(min(mt_a_mesh_z)) + 1);
            y_mat = log(mt_coh - min(min(mt_coh)) + 1);
            st_title = 'coh(a,z) log scale';
            st_ylabel = 'log(Cash-on-hand(a, z) - min(coh) + 1)';
            st_xlabel = 'log(a - min(a) + 1)';

            fl_b_bd_graph = log(fl_b_bd  - min(min(mt_a_mesh_z)) + 1);
            fl_borr_yminbd_graph = log(fl_borr_yminbd - min(min(mt_a_mesh_z)) + 1);
            fl_borr_ymaxbd_graph = log(fl_borr_ymaxbd - min(min(mt_a_mesh_z)) + 1);
        end

        % plot main x and y
        chart = plot(x_mat, y_mat, 'blue');

        % add color based on z
        clr = jet(numel(chart));
        for m = 1:numel(chart)
            set(chart(m), 'Color', clr(m,:))
        end

        %     if (length(ar_w_level_full) <= 100)
        %         scatter(ar_a_meshk, ar_k_mesha, 3, 'filled', ...
        %             'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        %     end
        %     if (length(ar_w_level_full) <= 100)
        %         gf_invalid_scatter = scatter(ar_a_meshk_full(ar_bl_wkb_invalid),...
        %                                      ar_k_mesha_full(ar_bl_wkb_invalid),...
        %                 20, 'O', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black');
        %     end

        % add various borrowing bound lines

        % add 0 lines
        xline(0);
        yline(0);

        % add 45 degrees line
        hline = refline([1 0]);
        hline.Color = 'k';
        hline.LineStyle = ':';
        hline.HandleVisibility = 'off';
        hline.LineWidth = 2.5;

        title(st_title)
        ylabel(st_ylabel)

        grid on;
        grid minor;

        legend2plot = fliplr([1 round(numel(chart)/3) round((2*numel(chart))/4)  numel(chart)]);
        legendCell = cellstr(num2str(ar_z_wage', 'z=%3.2f'));
        chart(length(chart)+1) = hline;
        legendCell{length(legendCell) + 1} = 'if coh(a,z) >= a';
        legend2plot = [legend2plot length(legendCell)];

        % if borrow plot additional borrowing bound lines
        if (fl_b_bd >= 0 )
            ar_legend_ele = [legend2plot];
            xlabel({st_xlabel})
        else
            % add fl_b_bd exo borrow line
            if (fl_b_bd >= min(ar_a))
                xline_borrbound = xline(fl_b_bd_graph);
                xline_borrbound.HandleVisibility = 'on';
                xline_borrbound.LineStyle = '-';
                xline_borrbound.Color = 'black';
                xline_borrbound.LineWidth = 2.5;

                yline_borrbound = yline(fl_b_bd_graph);
                yline_borrbound.HandleVisibility = 'off';
                yline_borrbound.LineStyle = '-';
                yline_borrbound.Color = 'black';
                yline_borrbound.LineWidth = 1;
            end

            xline_yminbd = xline(fl_borr_yminbd_graph);
            xline_yminbd.HandleVisibility = 'on';
            xline_yminbd.LineStyle = '--';
            xline_yminbd.Color = 'red';
            xline_yminbd.LineWidth = 2.5;

            yline_yminbd = yline(fl_borr_yminbd_graph);
            yline_yminbd.HandleVisibility = 'off';
            yline_yminbd.LineStyle = '--';
            yline_yminbd.Color = 'red';
            yline_yminbd.LineWidth = 1;

            if (bl_default)
                xline_ymaxbd = xline(fl_borr_ymaxbd_graph);
                xline_ymaxbd.HandleVisibility = 'on';
                xline_ymaxbd.LineStyle = '--';
                xline_ymaxbd.Color = 'blue';
                xline_ymaxbd.LineWidth = 2.5;

                yline_ymaxbd = yline(fl_borr_ymaxbd_graph);
                yline_ymaxbd.HandleVisibility = 'on';
                yline_ymaxbd.LineStyle = ':';
                yline_ymaxbd.Color = 'blue';
                yline_ymaxbd.LineWidth = 2.5;
            end


            % add bound line legends
            it_addlines_cn = 0;
            if (fl_b_bd >= min(ar_a))
                it_addlines_cn = it_addlines_cn + 1;
                chart(length(chart)+1) = xline_borrbound;
                legendCell{length(legendCell) + 1} = 'exo. borrow bound fbbd';
            end
            it_addlines_cn = it_addlines_cn + 1;
            chart(length(chart)+1) = xline_yminbd;
            legendCell{length(legendCell) + 1} = 'neg min inc: -zmin*w/r (no default)';
            if (bl_default)
                it_addlines_cn = it_addlines_cn + 1;
                chart(length(chart)+1) = xline_ymaxbd;
                legendCell{length(legendCell) + 1} = 'neg max inc: -zmax*w/r (default)';
                it_addlines_cn = it_addlines_cn + 1;
                chart(length(chart)+1) = yline_ymaxbd;
                legendCell{length(legendCell) + 1} = 'must default if coh(a,z)<dot-line';
            end

            % draw legend
            ar_legend_ele = [legend2plot length(legendCell)-it_addlines_cn:1:length(legendCell)];
            xlabel({st_xlabel 'if coh(a,z) < a, then a''(a,z)<a'})
        end

        % draw legends
        legend(chart(unique(ar_legend_ele)), legendCell(unique(ar_legend_ele)), 'Location', 'northwest');

    end

end

%% Display

if (bl_display_funcgrids)

    disp('ar_z_wage');
    disp(size(ar_z_wage));
    disp(ar_z_wage);
    
    disp('mt_z_trans');
    disp(size(mt_z_trans));
    disp(mt_z_trans);

    disp('ar_forbrblk, ar_forbrblk_r');
    disp(size(ar_forbrblk));
    disp([ar_forbrblk;ar_forbrblk_r]');

    param_map_keys = keys(func_map);
    param_map_vals = values(func_map);
    for i = 1:length(func_map)
        st_display = strjoin(['pos =' num2str(i) '; key =' string(param_map_keys{i}) '; val =' func2str(param_map_vals{i})]);
        disp(st_display);
    end
end

end
