%% ff_vf_az solves the one asset one shock asset problem
function [mt_val, mt_pol, flag] = ff_az_vf(varargin)

%% Check Parameters
params_len = length(varargin);
if params_len > 4
    error('ff_az_vf:Can only have 4 container map parameters');
end

%% Parameters  Defaults
it_param_set = 1;
[param_map, support_map] = ffs_az_default_param(it_param_set);
[armt_map] = ffs_az_grids(param_map, support_map);
default_params = {param_map support_map armt_map};

%% Override
% if varargin only has param_map and support_map,
[default_maps{1:params_len}] = varargin{:};
param_map = [param_map; default_maps{1}];
support_map = [support_map; default_maps{2}];
if params_len >= 1 && params_len <= 2
    [armt_map, func_map] = get_grid_func(param_map, support_map);
else
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

%% Retrieve Parameters from Map
params_group = values(param_map, {'fl_crra', 'fl_beta', 'fl_rho', 'fl_sig'});
[fl_crra, fl_beta, fl_rho, fl_sig] = params_group{:};
params_group = values(param_map, {'it_z_n', 'fl_z_mu', 'fl_z_rho', 'fl_z_sig', 'ar_z', 'mt_z_trans'});
[it_z_n, ar_z, mt_z_trans] = params_group{:};
params_group = values(param_map, {'it_a_n', 'ar_a'});
[it_a_n, ar_a] = params_group{:};
params_group = values(param_map, {'fl_w', 'fl_r'});
[fl_w, fl_r] = params_group{:};
params_group = values(param_map, {'it_maxiter_val', 'fl_tol_val', 'fl_tol_pol'});
[it_maxiter_val, fl_tol_val, fl_tol_pol] = params_group{:};

% support
params_group = values(support_map, {'bl_display', 'bl_graph', 'bl_graph_onebyones'});
[bl_display, bl_graph, bl_graph_onebyones] = params_group{:};
params_group = values(support_map, {'bl_time', 'bl_profile', 'st_profile_path'});
[bl_time, bl_profile, st_profile_path] = params_group{:};


%% Profiling start
bl_profile = true;
if (bl_profile)
    close all;
    profile off;
    profile on;
end


%% Initialize value and policy function
mt_val0 = zeros(length(ar_a),length(ar_z));
mt_val1 = mt_val0;
mt_pol1 = zeros(length(ar_a),length(ar_z));

fl_diff = 1;
it_iter = 1;

%% Define Equations and Parameters
f_cons = @(z, a, aprime)(fl_w*z + (1+fl_r)*a - aprime);
f_incm = @(z, a)        (fl_w*z + fl_r*a);
f_util = @(c)           ((c)^(1-fl_crra)-1)/(1-fl_crra);
u_neg_c = -1000;

%% Timing Starts
if (bl_time_vf_okz); tic; end
%% Loop Solution
while fl_diff > fl_tol_val && it_iter <= it_maxiter_val
    for i = 1:length(ar_z) %outer loop- current z
        z = ar_z(i);
        for j = 1:length(ar_a) %first inner loop - current assets
            a = ar_a(j);
            vnew = zeros(size(ar_a));
            for x = 1:length(ar_a) %second inner loop - chosen a
                aprime = ar_a(x);
%                 c = (fl_w*z + (1+fl_r)*a - aprime);
%                 c = cons_az(fl_w, fl_r, a, z, aprime);
                  c = f_cons(z, a, aprime);
                if c > 0
					if (fl_crra == 1)
						vnew(x) = log(c);
					else
						vnew(x) = f_util(c);
					end
                    for y = 1:length(ar_z)
                        vnew(x) = vnew(x) + fl_beta*mt_z_trans(i,y)*mt_val0(x,y);
                    end
                else
                    vnew(x) = u_neg_c;
                end
            end
            max_a_position = find(vnew == max(vnew));
            mt_val1(j,i) = vnew(max_a_position(1));
            mt_pol1(j,i) = ar_a(max_a_position(1));
        end
    end
    it_iter = it_iter + 1;
    fl_diff = norm(mt_val1 - mt_val0);
    mt_val0 = mt_val1;
end

%% Outputs
[mt_val, mt_pol] = deal(mt_val1, mt_pol1);
if (support_map('bl_save_mat_vfi'))
    % use broadcasting
    ar_a_mby1 = reshape(ar_a, [length(ar_a), 1]);
    ar_z_1byn = reshape(ar_z, [1, length(ar_z)]);
    mt_cons = f_cons(ar_z_1byn, ar_a_mby1, mt_pol);
    mt_incm = f_incm(ar_z_1byn, ar_a_mby1);
%     mt_cons = cons_az(fl_w, fl_r, ar_a_mby1, ar_z_1byn, mt_pol);
%     mt_incm = incm_az(fl_w, fl_r, ar_a_mby1, ar_z_1byn);
end

%% Flag
if fl_diff <= fl_tol_val && it_iter > it_maxiter_val
%     tolerance did not reach, iteration bound reached
    flag = 1;
elseif fl_diff > fl_tol_val && it_iter <= it_maxiter_val
%     tolerance reached, iter bound did not reach
    flag = 2;
else
%     both reached
    flag = 0;
end


%% Timing Ends
if (bl_time_vf_okz); toc; end

%% Profiling
if (bl_profile)
    profile off
    profile viewer
    st_profile_file = ['C:/Users/fan/ThaiForInfLuuRobFan/matlab/inf_okz/solve/_profile/vf_okz_interp_p' num2str(it_cm)];
    profsave(profile('info'), st_profile_file);
end

%% Graphing Final
if (bl_graph_vf_okz)
    vf_okz_grh(param_map, support_map, armt_map, mt_val, mt_pol_wkb, mt_pol_kap)
end

%% Save Workspace in Mat File
if (support_map('bl_save_mat_vfi'))
    st_file_name = ['vf_' support_map('st_file_name')];
    save(ff_sup_save_prep(support_map('st_path_folder'), st_file_name));
end


end
