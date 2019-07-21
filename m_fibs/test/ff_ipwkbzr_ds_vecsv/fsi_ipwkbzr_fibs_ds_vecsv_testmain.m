%% Test Risky + Safe Asset (Save + Borr + FIBS + RShock) Interp-Perc Main Function
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function result_map = ff_ipwkbzr_testmain(varargin)
%% F main formal informal borrow save tester
% Controls four key groups of borrowing and savings parameters, allows for
% loops over various parameters, simulate model over multiple loops.
% Functiond designed to be invoked to compare when bl_bridge is true vs
% false.
%
% Note for the default parameters here, the formal and informal borrowing
% interest rates are the same. Additionally there is no bridge loan
% (bl_bridge = false). Given these parameters, this simulation should give
% the same results as the _ipwkbzr_ simulation with single borrowing source.
% (see @seealso below).
%
% @param it_simu_vec_len scalar number of elements each array to simulate
%
% @param ar_it_test_grp array which groups to simulate over, see below
%
% @seealso
%
% * *abz* test interest rate no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_ipwkbzr_ds_vecsv_nbc.html fsi_ipwkbzr_ds_vecsv_nbc>
% * *abz* test interest rate default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_ipwkbzr_ds_vecsv_default.html fsi_ipwkbzr_ds_vecsv_default>
%

%% Default Parameter Loops

% Loops to Vary for Simulations Below
it_size_type = 4;
ar_it_test_grp = [3, 8, 9];
it_simu_vec_len= 3;

%% Default Parameters

% Base parameters call
it_param_set = 9;
[param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);

% Preference Parameters
param_map('fl_beta') = 0.94;
param_map('fl_crra') = 1.5;

% Shock Parameters
param_map('fl_z_rho') = 0.8;
param_map('fl_z_sig') = 0.2;

% Key Borrowing Controls
param_map('bl_default') = true;
param_map('bl_bridge') = true;
param_map('bl_rollover') = true;

% Key borrowing float controls
param_map('fl_default_wprime') = 0;
param_map('fl_b_bd') = -20;
param_map('fl_c_min') = 0.02;

% Interest Rate Simulate Controls
param_map('fl_r_fsv') = 0.025;
param_map('fl_r_fbr') = 0.065;

% Informal Rate
param_map('st_z_r_infbr_drv_ele_type') = 'unif';
param_map('st_z_r_infbr_drv_prb_type') = 'poiss';
param_map('fl_z_r_infbr_poiss_mean') = 20;
param_map('fl_z_r_infbr_max') = 0.095;
param_map('fl_z_r_infbr_min') = 0.025;
param_map('fl_z_r_infbr_n') = 5;

% Borrowing Grid Simulate Controls
param_map('st_forbrblk_type') = 'unif';
param_map('fl_forbrblk_brmost') = -19;
param_map('fl_forbrblk_brleast') = -1;
param_map('fl_forbrblk_gap') = -1.5;

% Display Parameters
support_map('bl_display') = false;
support_map('bl_display_final') = false;
support_map('bl_time') = true;
support_map('bl_profile') = false;
support_map('bl_graph_funcgrids') = false;

support_map('bl_mat_test') = true;
st_matimg_path_root = support_map('st_matimg_path_root');
% test_borinf for default: ar_it_test_grp = [3, 8, 9]
support_map('st_mat_test_path') = [st_matimg_path_root '/test/ff_ipwkbzr_ds_vecsv/test_forinf/mat/'];

%% Array Parameters

% Default
default_params = {param_map support_map ...
                    it_size_type ar_it_test_grp it_simu_vec_len};

%% Parse Inputs
% override default set above if any parameters are updated

params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};
param_map = [param_map; default_params{1}];
support_map = [support_map; default_params{2}];
it_size_type = default_params{3};
ar_it_test_grp = default_params{4};
it_simu_vec_len = default_params{5};

% support_map
support_map('st_mat_test_prefix') = [''];
support_map('st_mat_test_name_main') = ['res'];
support_map('st_mat_test_suffix') = ['g' strrep(num2str(ar_it_test_grp), '  ', '') ...
    '_t' num2str(it_size_type) 'l' num2str(it_simu_vec_len)];

%% Set Solve Sizes

if (it_size_type == 1)
    
    % Basic Test Run    
    param_map('it_w_perc_n') = 10;
    param_map('it_ak_perc_n') = param_map('it_w_perc_n');
    param_map('it_coh_bridge_perc_n') = param_map('it_w_perc_n');

    param_map('fl_coh_interp_grid_gap') = 2;
    param_map('it_c_interp_grid_gap') = 10^-1;
    param_map('fl_w_interp_grid_gap') = 2;

    param_map('it_z_wage_n') = 2;
    param_map('fl_z_r_infbr_n') = 2;
    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');

    param_map('it_maxiter_val') = 5;
    
elseif (it_size_type == 2)
    
    % Usable Run
    param_map('it_maxiter_val') = 50;
    
elseif (it_size_type == 3)
    
    % Full Run
    
elseif (it_size_type == 4)
    
    % Denser Run
    param_map('it_w_perc_n') = 100;
    param_map('it_ak_perc_n') = param_map('it_w_perc_n');
    % ususally divide by 5
    param_map('it_coh_bridge_perc_n') = round(param_map('it_w_perc_n')/3);

end

%% Set up Arrays
% param_test_array_map

param_test_array_map = containers.Map('KeyType','char', 'ValueType','any');
param_test_array_map('ar_it_z_wage_n') = [5];
param_test_array_map('ar_fl_b_bd') = linspace(-20, -5, it_simu_vec_len);
param_test_array_map('ar_fl_c_min') = linspace(0.1, 0.001, it_simu_vec_len);
param_test_array_map('ar_fl_beta') = linspace(0.94, 0.98, it_simu_vec_len);
param_test_array_map('ar_fl_crra') = linspace(1, 2, it_simu_vec_len);
param_test_array_map('ar_fl_z_rho') = linspace(0.65, 0.95, it_simu_vec_len);
param_test_array_map('ar_fl_z_sig') = linspace(0.05, 0.35, it_simu_vec_len);

param_test_array_map('ar_fl_z_r_infbr_poiss_mean') = linspace(1, 20, it_simu_vec_len);
param_test_array_map('ar_fl_r_fbr') = linspace(0.030, 0.080, it_simu_vec_len);
param_test_array_map('ar_fl_forbrblk_gap') = linspace(-0.1, -2.5, it_simu_vec_len);

%% Parase Preference and Shock Parameters

% preference Parameters
params_group = values(param_map, {'fl_beta', 'fl_crra'});
[fl_beta, fl_crra] = params_group{:};

% shock Parameters
params_group = values(param_map, {'fl_z_rho', 'fl_z_sig'});
[fl_z_rho, fl_z_sig] = params_group{:};

%% Parse Borrowing and Savings Parameters
% param_map, borrow and save controls:
%
% # boolean controls
% # float controls
% # interest rates
% # formal menu/grid
%

% boolean controls
params_group = values(param_map, {'bl_default', 'bl_bridge', 'bl_rollover'});
[bl_default, bl_bridge, bl_rollover] = params_group{:};

% float controls
params_group = values(param_map, {'fl_default_wprime', 'fl_b_bd', 'fl_c_min'});
[fl_default_wprime, fl_b_bd, fl_c_min] = params_group{:};

% interest rates
params_group = values(param_map, {'fl_r_fsv', 'fl_z_r_infbr_poiss_mean', 'fl_r_fbr'});
[fl_r_fsv, fl_z_r_infbr_poiss_mean, fl_r_fbr] = params_group{:};

% formal menu controls
params_group = values(param_map, {'st_forbrblk_type', 'fl_forbrblk_brmost', 'fl_forbrblk_brleast', 'fl_forbrblk_gap'});
[st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap] = params_group{:};

% support map
params_group = values(support_map, {'bl_mat_test', 'st_mat_test_path', 'st_mat_test_prefix', 'st_mat_test_name_main', 'st_mat_test_suffix'});
[bl_mat_test, st_mat_test_path, st_mat_test_prefix, st_mat_test_name_main, st_mat_test_suffix] = params_group{:};

%% Parse Parameter Arrays

params_group = values(param_test_array_map, {...
    'ar_it_z_wage_n', ...
    'ar_fl_b_bd', 'ar_fl_c_min', 'ar_fl_z_r_infbr_poiss_mean', ...
    'ar_fl_beta', 'ar_fl_crra', 'ar_fl_z_rho', 'ar_fl_z_sig', ...
    'ar_fl_r_fbr', 'ar_fl_forbrblk_gap'});
[ar_it_z_wage_n, ...
    ar_fl_b_bd, ar_fl_c_min, ar_fl_z_r_infbr_poiss_mean, ...
    ar_fl_beta, ar_fl_crra, ar_fl_z_rho, ar_fl_z_sig, ...
    ar_fl_r_fbr, ar_fl_forbrblk_gap] = params_group{:};

%% Test Various Arrays of Parameters

close all;

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxx Preference xxxxxx');
disp(['fl_beta = ' num2str(fl_beta)]);
disp(['fl_crra = ' num2str(fl_crra)]);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxx shocks xxxxxx');
disp(['fl_z_rho = ' num2str(fl_z_rho)]);
disp(['fl_z_sig = ' num2str(fl_z_sig)]);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxx Boolean Controls xxxxxx');
disp(['bl_default = ' num2str(bl_default)]);
disp(['bl_bridge = ' num2str(bl_bridge)]);
disp(['bl_rollover = ' num2str(bl_rollover)]);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxx Float Controls xxxxxxxx');
disp(['fl_default_wprime = ' num2str(fl_default_wprime)]);
disp(['fl_b_bd = ' num2str(fl_b_bd)]);
disp(['fl_c_min = ' num2str(fl_c_min)]);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxx Interest Rates xxxxxxxx');
disp(['fl_r_fsv = ' num2str(fl_r_fsv)]);
disp(['fl_z_r_infbr_poiss_mean = ' num2str(fl_z_r_infbr_poiss_mean)]);
disp(['fl_r_fbr = ' num2str(fl_r_fbr)]);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxx Formal Menu Params xxxx');
disp(['st_forbrblk_type = ' num2str(st_forbrblk_type)]);
disp(['fl_forbrblk_brmost = ' num2str(fl_forbrblk_brmost)]);
disp(['fl_forbrblk_brleast = ' num2str(fl_forbrblk_brleast)]);
disp(['fl_forbrblk_gap = ' num2str(fl_forbrblk_gap)]);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
if (ismember(1, ar_it_test_grp))
    disp(['ar_fl_b_bd = ' num2str(ar_fl_b_bd)]);
elseif (ismember(2, ar_it_test_grp))
    disp(['ar_fl_c_min = ' num2str(ar_fl_c_min)]);
elseif (ismember(3, ar_it_test_grp))
    disp(['ar_fl_z_r_infbr_poiss_mean = ' num2str(ar_fl_z_r_infbr_poiss_mean)]);
elseif (ismember(4, ar_it_test_grp))
    disp(['ar_fl_beta = ' num2str(ar_fl_beta)]);
elseif (ismember(5, ar_it_test_grp))
    disp(['ar_fl_crra = ' num2str(ar_fl_crra)]);
elseif (ismember(6, ar_it_test_grp))
    disp(['ar_fl_z_rho = ' num2str(ar_fl_z_rho)]);
elseif (ismember(7, ar_it_test_grp))
    disp(['ar_fl_z_sig = ' num2str(ar_fl_z_sig)]);
elseif (ismember(8, ar_it_test_grp))
    disp(['ar_fl_r_fbr = ' num2str(ar_fl_r_fbr)]);    
elseif (ismember(9, ar_it_test_grp))
    disp(['ar_fl_forbrblk_gap = ' num2str(ar_fl_forbrblk_gap)]);    
end

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');

cl_tb_outcomes_meansdperc_wthinfo = cell([length(ar_it_test_grp), 1]);

for ar_it_test_ctr = 1:length(ar_it_test_grp)
    
    it_test_grp = ar_it_test_grp(ar_it_test_ctr);
    
    disp('---------------------------');
    disp('---------------------------');
    disp('---------------------------');
    disp('---------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');   
    
    %% Display Current Parameter been Varied
    % Conditional Changes
    if (it_test_grp == 1)
        % Change the borrowing bound
        disp(['Vary Borrow Bound']);
        disp(['ar_fl_b_bd = ' num2str(ar_fl_b_bd)]);
    elseif (it_test_grp == 2)
        % Change minimum consumption
        disp(['Vary Minimum Consumption']);
        disp(['ar_fl_c_min = ' num2str(ar_fl_c_min)]);
    elseif (it_test_grp == 3)
        disp(['Vary Interest rate']);
        disp(['ar_fl_z_r_infbr_poiss_mean = ' num2str(ar_fl_z_r_infbr_poiss_mean)]);
    elseif (it_test_grp == 4)
        disp(['Vary Discount']);
        disp(['ar_fl_beta = ' num2str(ar_fl_beta)]);
    elseif (it_test_grp == 5)
        disp(['Vary Risk Aversion']);
        disp(['ar_fl_crra = ' num2str(ar_fl_crra)]);
    elseif (it_test_grp == 6)
        disp(['Vary Shock Persistence']);
        disp(['ar_fl_z_rho = ' num2str(ar_fl_z_rho)]);
    elseif (it_test_grp == 7)
        disp(['Vary Shock Variance']);
        disp(['ar_fl_z_sig = ' num2str(ar_fl_z_sig)]);
    elseif (it_test_grp == 8)
        disp(['Vary Formal Borrowing Rate']);
        disp(['ar_fl_r_fbr = ' num2str(ar_fl_r_fbr)]);        
    elseif (it_test_grp == 9)
        disp(['Vary Formal Borrowing Blocks']);
        disp(['ar_fl_forbrblk_gap = ' num2str(ar_fl_forbrblk_gap)]);        
    end
    
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    tb_outcomes_meansdperc_wthinfo_stack = [];
    for it_cur_param = 1:1:it_simu_vec_len

        %% Adjust Value for Current Parameter been Varied
        
        % Reset Base Parameters, parameters already grabbed out, updating
        % param_map does not impact fl_b_bd etc..
        param_map('fl_b_bd') = fl_b_bd;
        param_map('fl_z_r_infbr_poiss_mean') = fl_z_r_infbr_poiss_mean;
        param_map('fl_c_min') = fl_c_min;

        % Conditional Changes
        if (it_test_grp == 1)
            % Change the borrowing bound
            disp(['xxxxx fl_b_bd = ' num2str(ar_fl_b_bd(it_cur_param)) ' xxxxx']);
            % Update fl_b_bd
            param_map('fl_b_bd') = ar_fl_b_bd(it_cur_param);
        elseif (it_test_grp == 2)
            % Change minimum consumption
            disp(['xxxxx fl_c_min = ' num2str(ar_fl_c_min(it_cur_param)) ' xxxxx']);
            % Update fl_c_min
            param_map('fl_c_min') = ar_fl_c_min(it_cur_param);
        elseif (it_test_grp == 3)
            % Change the informal interest rate
            disp(['xxxxx fl_z_r_infbr_poiss_mean = ' num2str(ar_fl_z_r_infbr_poiss_mean(it_cur_param)) ' xxxxx']);
            % Update
            param_map('fl_z_r_infbr_poiss_mean') = ar_fl_z_r_infbr_poiss_mean(it_cur_param);
        elseif (it_test_grp == 4)
            % Change the discount
            disp(['xxxxx fl_beta = ' num2str(ar_fl_beta(it_cur_param)) ' xxxxx']);
            % Update
            param_map('fl_beta') = ar_fl_beta(it_cur_param);
        elseif (it_test_grp == 5)
            % Change the risk aversion
            disp(['xxxxx fl_crra = ' num2str(ar_fl_crra(it_cur_param)) ' xxxxx']);
            % Update
            param_map('fl_crra') = ar_fl_crra(it_cur_param);
        elseif (it_test_grp == 6)
            % Change the shock persistence
            disp(['xxxxx fl_z_rho = ' num2str(ar_fl_z_rho(it_cur_param)) ' xxxxx']);
            % Update
            param_map('fl_z_rho') = ar_fl_z_rho(it_cur_param);
        elseif (it_test_grp == 7)
            % Change the shock variance
            disp(['xxxxx fl_z_sig = ' num2str(ar_fl_z_sig(it_cur_param)) ' xxxxx']);
            % Update
            param_map('fl_z_sig') = ar_fl_z_sig(it_cur_param);
        elseif (it_test_grp == 8)
            % Change the shock variance
            disp(['xxxxx fl_r_fbr = ' num2str(ar_fl_r_fbr(it_cur_param)) ' xxxxx']);
            % Update
            param_map('fl_r_fbr') = ar_fl_r_fbr(it_cur_param);
        elseif (it_test_grp == 9)
            % Change the shock variance
            disp(['xxxxx fl_forbrblk_gap = ' num2str(ar_fl_forbrblk_gap(it_cur_param)) ' xxxxx']);
            % Update
            param_map('fl_forbrblk_gap') = ar_fl_forbrblk_gap(it_cur_param);
        end
        
        %% Simulate Model
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');
        it_w_perc_n = param_map('it_w_perc_n');
        it_z_n = param_map('it_z_n');
        disp(['xxxxx it_w_perc_n = ' num2str(it_w_perc_n) ', it_z_n = ' num2str(it_z_n) ' xxxxx']);
        % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_ipwkbzr_fibs_get_funcgrid.html ffs_ipwkbzr_fibs_get_funcgrid>
        [armt_map, func_map] = ffs_ipwkbzr_fibs_get_funcgrid(param_map, support_map);
        % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_ipwkbzr_fibs_vf_vecsv.html ff_ipwkbzr_fibs_vf_vecsv>
        result_map = ff_ipwkbzr_fibs_vf_vecsv(param_map, support_map, armt_map, func_map);
        % Call Distribution CProgram
        result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map);

        %% Store Results to Table
        % Append Simulation Loop Statistics 
        tb_outcomes_meansdperc = result_map('tb_outcomes_meansdperc');            
        ar_simu_info = [it_test_grp, it_cur_param, param_map('fl_z_r_infbr_poiss_mean'), param_map('fl_r_fbr'), param_map('fl_forbrblk_gap')];
        tb_simu_info = array2table(ar_simu_info + zeros([size(tb_outcomes_meansdperc,1), length(ar_simu_info)]));
        cl_col_names = {'it_test_grp', 'it_cur_param', ...
                        'fl_z_r_infbr_poiss_mean', 'fl_r_fbr', 'fl_forbrblk_gap'};
        tb_simu_info.Properties.VariableNames = cl_col_names;
        variablenames = tb_outcomes_meansdperc.Properties.RowNames;
        tb_simu_info.Properties.RowNames = variablenames;

        % Collect statistics
        tb_outcomes_meansdperc_wthinfo = [tb_simu_info tb_outcomes_meansdperc];
        tb_outcomes_meansdperc_wthinfo.Properties.RowNames = ...
            strcat(variablenames, '_p', num2str(it_cur_param));
        tb_outcomes_meansdperc_wthinfo = addvars(tb_outcomes_meansdperc_wthinfo, variablenames, 'Before', 1);
        
        %% Combine Results from Different Parameters in Shared Table
        if (it_cur_param == 1)
            tb_outcomes_meansdperc_wthinfo_stack = tb_outcomes_meansdperc_wthinfo;
        else 
            tb_outcomes_meansdperc_wthinfo_stack = [tb_outcomes_meansdperc_wthinfo_stack; tb_outcomes_meansdperc_wthinfo];
        end
                
    end
    
    cl_tb_outcomes_meansdperc_wthinfo{ar_it_test_ctr} = tb_outcomes_meansdperc_wthinfo_stack;

end

%% Save Mat

if (bl_mat_test)
    clear armt_map result_map
    if ~exist(st_mat_test_path,'dir'); mkdir(st_mat_test_path); end
    st_file_name = [st_mat_test_prefix st_mat_test_name_main st_mat_test_suffix];
    save(strcat(st_mat_test_path, st_file_name));
end

% close all
close all;

end
