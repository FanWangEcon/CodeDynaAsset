%% Test For+inf+Borr+Save Borrow Bound Interest Rate Etc
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function result_map = fsi_abz_fibs_ds_vecsv_testmain(varargin)
%% FSI_ABZ_FIBS_DS_VECSV_TESTMAIN main formal informal borrow save tester
% Controls four key groups of borrowing and savings parameters, allows for
% loops over various parameters, simulate model over multiple loops.
% Functiond designed to be invoked to compare when bl_bridge is true vs
% false.
%
% Note for the default parameters here, the formal and informal borrowing
% interest rates are the same. Additionally there is no bridge loan
% (bl_bridge = false). Given these parameters, this simulation should give
% the same results as the _abz_ simulation with single borrowing source.
% (see @seealso below).
%
% @param it_simu_vec_len scalar number of elements each array to simulate
%
% @param ar_it_test_grp array which groups to simulate over, see below
%
% @seealso
%
% * *abz* test interest rate no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc.html fsi_abz_ds_vecsv_nbc>
% * *abz* test interest rate default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default.html fsi_abz_ds_vecsv_default>
%

%% Default Parameters

% Base parameters call
it_param_set = 9;
[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

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
param_map('fl_default_aprime') = 0;
param_map('fl_b_bd') = -12.5;
param_map('fl_c_min') = 0.025;

% Interest Rate Simulate Controls
param_map('fl_r_fsv') = 0.025;
param_map('fl_r_inf') = 0.10;
param_map('fl_r_inf_bridge') = param_map('fl_r_inf');
param_map('fl_r_fbr') = 0.10;

% Borrowing Grid Simulate Controls
param_map('st_forbrblk_type') = 'seg3';
param_map('fl_forbrblk_brmost') = -19;
param_map('fl_forbrblk_brleast') = -1;
param_map('fl_forbrblk_gap') = -1.5;

% Display Parameters
support_map('bl_display') = false;
support_map('bl_display_final') = false;
support_map('bl_time') = true;
support_map('bl_profile') = false;
support_map('bl_graph_funcgrids') = false;

%% Default Parameter Loops

% Loops to Vary for Simulations Below
it_simu_vec_len = 4;
ar_it_test_grp = [1,2,3];

% param_test_array_map 
param_test_array_map = containers.Map('KeyType','char', 'ValueType','any');
param_test_array_map('ar_it_a_n') = [750];
param_test_array_map('ar_it_z_n') = [15];
param_test_array_map('ar_fl_b_bd') = linspace(-20, -5, it_simu_vec_len);
param_test_array_map('ar_fl_c_min') = linspace(0.1, 0.001, it_simu_vec_len);
param_test_array_map('ar_fl_r_inf') = linspace(0.20, 0.01, it_simu_vec_len);
param_test_array_map('ar_fl_beta') = linspace(0.94, 0.98, it_simu_vec_len);
param_test_array_map('ar_fl_crra') = linspace(1, 2, it_simu_vec_len);
param_test_array_map('ar_fl_z_rho') = linspace(0.65, 0.95, it_simu_vec_len);
param_test_array_map('ar_fl_z_sig') = linspace(0.05, 0.35, it_simu_vec_len);

% Default
default_params = {param_map support_map ...
                    ar_it_test_grp it_simu_vec_len ...
                    param_test_array_map};

%% Parse Inputs
% override default set above if any parameters are updated

params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};
param_map = [param_map; default_params{1}];
support_map = [support_map; default_params{2}];
ar_it_test_grp = default_params{3};
it_simu_vec_len = default_params{4};
param_test_array_map = [param_test_array_map; default_params{5}];

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
params_group = values(param_map, {'fl_default_aprime', 'fl_b_bd', 'fl_c_min'});
[fl_default_aprime, fl_b_bd, fl_c_min] = params_group{:};

% interest rates
params_group = values(param_map, {'fl_r_fsv', 'fl_r_inf', 'fl_r_inf_bridge', 'fl_r_fbr'});
[fl_r_fsv, fl_r_inf, fl_r_inf_bridge, fl_r_fbr] = params_group{:};

% formal menu controls
params_group = values(param_map, {'st_forbrblk_type', 'fl_forbrblk_brmost', 'fl_forbrblk_brleast', 'fl_forbrblk_gap'});
[st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap] = params_group{:};

%% Parse Parameter Arrays

params_group = values(param_test_array_map, {...
    'ar_it_a_n', 'ar_it_z_n', ...
    'ar_fl_b_bd', 'ar_fl_c_min', 'ar_fl_r_inf', ...
    'ar_fl_beta', 'ar_fl_crra', 'ar_fl_z_rho', 'ar_fl_z_sig'});
[ar_it_a_n, ar_it_z_n, ...
    ar_fl_b_bd, ar_fl_c_min, ar_fl_r_inf, ...
    ar_fl_beta, ar_fl_crra, ar_fl_z_rho, ar_fl_z_sig] = params_group{:};

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
disp(['fl_default_aprime = ' num2str(fl_default_aprime)]);
disp(['fl_b_bd = ' num2str(fl_b_bd)]);
disp(['fl_c_min = ' num2str(fl_c_min)]);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxx Interest Rates xxxxxxxx');
disp(['fl_r_fsv = ' num2str(fl_r_fsv)]);
disp(['fl_r_inf = ' num2str(fl_r_inf)]);
disp(['fl_r_inf_bridge = ' num2str(fl_r_inf_bridge)]);
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
    disp(['ar_fl_r_inf = ' num2str(ar_fl_r_inf)]);
elseif (ismember(4, ar_it_test_grp))
    disp(['ar_fl_beta = ' num2str(ar_fl_beta)]);
elseif (ismember(5, ar_it_test_grp))
    disp(['ar_fl_crra = ' num2str(ar_fl_crra)]);
elseif (ismember(6, ar_it_test_grp))
    disp(['ar_fl_z_rho = ' num2str(ar_fl_z_rho)]);
elseif (ismember(7, ar_it_test_grp))
    disp(['ar_fl_z_sig = ' num2str(ar_fl_z_sig)]);
end
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');


for it_test_grp = ar_it_test_grp

    disp('---------------------------');
    disp('---------------------------');
    disp('---------------------------');
    disp('---------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');

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
        disp(['ar_fl_r_inf = ' num2str(ar_fl_r_inf)]);
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
    end

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');


    for it_cur_param = 1:1:it_simu_vec_len

        % Reset Base Parameters, parameters already grabbed out, updating
        % param_map does not impact fl_b_bd etc..
        param_map('fl_b_bd') = fl_b_bd;
        param_map('fl_r_inf') = fl_r_inf;
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
            disp(['xxxxx fl_r_inf = ' num2str(ar_fl_r_inf(it_cur_param)) ' xxxxx']);
            % Update
            param_map('fl_r_inf') = ar_fl_r_inf(it_cur_param);
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
        end

        % Loop over potentially different accuracy levels and run file
        for it_accuracy = 1:length(ar_it_a_n)
            % Accuracy Regular
            param_map('it_a_n') = ar_it_a_n(it_accuracy);
            param_map('it_z_n') = ar_it_z_n(it_accuracy);
            disp(['xxxxx it_a_n = ' num2str(ar_it_a_n(it_accuracy)) ', it_z_n = ' num2str(ar_it_z_n(it_accuracy)) ' xxxxx']);
            % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_fibs_get_funcgrid.html ffs_abz_fibs_get_funcgrid>
            bl_input_override = true;
            [armt_map, func_map] = ffs_abz_fibs_get_funcgrid(param_map, support_map, bl_input_override);
            % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_fibs_vf_vecsv.html ff_abz_fibs_vf_vecsv>
            result_map = ff_abz_fibs_vf_vecsv(param_map, support_map, armt_map, func_map);
            % Call Distribution CProgram
            result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
        end

    end

end

% close all
close all;

end
