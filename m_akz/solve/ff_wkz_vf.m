%% 
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository> 
% Table of Content.*

function result_map = ff_wkz_vf(varargin)
%% FF_WKZ_VF solve infinite horizon exo shock + endo asset problem
% This program solves the infinite horizon dynamic savings and risky
% capital asset problem with some ar1 shock. This is the two step solution
% version of
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf.html
% ff_akz_vf>. See
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_wkz_evf.html
% ff_wkz_evf> for details about the second stage. 
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param armt_map container container with states, choices and shocks
% grids that are inputs for grid based solution algorithm
%
% @param func_map container container with function handles for
% consumption cash-on-hand etc.
%
% @return result_map container contains policy function matrix, value
% function matrix, iteration results, and policy function, value function
% and iteration results tables. 
%
% keys included in result_map:
%
% * mt_val matrix states_n by shock_n matrix of converged value function grid
% * mt_pol_a matrix states_n by shock_n matrix of converged policy function grid
% * ar_val_diff_norm array if bl_post = true it_iter_last by 1 val function
% difference between iteration
% * ar_pol_diff_norm array if bl_post = true it_iter_last by 1 policy
% function difference between iterations
% * mt_pol_perc_change matrix if bl_post = true it_iter_last by shock_n the
% proportion of grid points at which policy function changed between
% current and last iteration for each element of shock
%
% @example
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/paramfunc/ffs_akz_set_default_param.m ffs_akz_set_default_param>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/paramfunc/ffs_akz_get_funcgrid.m ffs_akz_get_funcgrid>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solvepost/ff_akz_vf_post.m ff_akz_vf_post>
%


%% Default
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button

it_param_set = 3;
bl_input_override = true;
[param_map, support_map] = ffs_akz_set_default_param(it_param_set);

% parameters can be set inside ffs_akz_set_default_param or updated here
param_map('it_w_n') = 50;
param_map('it_z_n') = 15;

% get armt and func map
[armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
default_params = {param_map support_map armt_map func_map};

%% Parse Parameters 1

% if varargin only has param_map and support_map,
params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};
param_map = [param_map; default_params{1}];
support_map = [support_map; default_params{2}];
if params_len >= 1 && params_len <= 2
    % If override param_map, re-generate armt and func if they are not
    % provided
    bl_input_override = true;
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% append function name
st_func_name = 'ff_wkz_vf';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'ar_w', 'ar_z'});
[ ar_w, ar_z] = params_group{:};
params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'mt_coh_wkb', 'it_ameshk_n'});
[ar_a_meshk, ar_k_mesha, mt_coh_wkb, it_ameshk_n] = params_group{:};
% func_map
params_group = values(func_map, {'f_util_log', 'f_util_crra', 'f_cons'});
[f_util_log, f_util_crra, f_cons] = params_group{:};
% param_map
params_group = values(param_map, {'fl_r_save', 'fl_r_borr', 'fl_w',...
    'it_z_n', 'fl_crra', 'fl_beta', 'fl_c_min'});
[fl_r_save, fl_r_borr, fl_wage, it_z_n, fl_crra, fl_beta, fl_c_min] = params_group{:};
params_group = values(param_map, {'it_maxiter_val', 'fl_tol_val', 'fl_tol_pol', 'it_tol_pol_nochange'});
[it_maxiter_val, fl_tol_val, fl_tol_pol, it_tol_pol_nochange] = params_group{:};
% support_map
params_group = values(support_map, {'bl_profile', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display', 'it_display_every', 'bl_post'});
[bl_profile, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display, it_display_every, bl_post] = params_group{:};

%% Initialize Output Matrixes

mt_val_cur = zeros(length(ar_a_meshk),length(ar_z));
mt_val = mt_val_cur - 1;
mt_pol_a = zeros(length(ar_a_meshk),length(ar_z));
mt_pol_a_cur = mt_pol_a - 1;
mt_pol_k = zeros(length(ar_a_meshk),length(ar_z));
mt_pol_k_cur = mt_pol_k - 1;

%% Initialize Convergence Conditions

bl_vfi_continue = true;
it_iter = 0;
ar_val_diff_norm = zeros([it_maxiter_val, 1]);
ar_pol_diff_norm = zeros([it_maxiter_val, 1]);
mt_pol_perc_change = zeros([it_maxiter_val, it_z_n]);

%% Iterate Value Function
% Loop solution with 4 nested loops
%
% # loop 1: over exogenous states
% # loop 2: over endogenous states
% # loop 3: over choices
% # loop 4: add future utility, integration--loop over future shocks
%

% Start Profile
if (bl_profile)
    close all;
    profile off;
    profile on;
end

% Start Timer
if (bl_time)
    tic;
end

% Value Function Iteration
while bl_vfi_continue
    it_iter = it_iter + 1;
    
    %% Solve Second Stage Problem k*(w,z)
    % This is the key difference between this function and
    % <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_functions.html
    % ffs_akz_set_functions> which solves the two stages jointly
    
    bl_input_override = true;
    [mt_ev_condi_z_max, ~, mt_ev_condi_z_max_kp, mt_ev_condi_z_max_bp] = ...
        ff_wkz_evf(mt_val_cur, param_map, support_map, armt_map, bl_input_override);
    
    %% Solve First Stage Problem w*(z) given k*(w,z)
    % loop 1: over exogenous states
    for it_z_i = 1:length(ar_z)
        
        % Get 2nd Stage Arrays
        ar_ev_condi_z_max_z = mt_ev_condi_z_max(:, it_z_i);        
        ar_w_kstar_z = mt_ev_condi_z_max_kp(:, it_z_i);
        ar_w_astar_z = mt_ev_condi_z_max_bp(:, it_z_i);        
        
        % loop 2: over endogenous states
        for it_coh_j = 1:length(ar_a_meshk)
            % Get cash-on-hand which include k,b,z
            fl_coh = mt_coh_wkb(it_coh_j, it_z_i);
            
            % loop 3: over choices, only w vector
            % we choose w(z), know from ff_wkz_evf k*(w,z), b*=w-k*
            ar_val_cur = zeros(size(ar_w));
            for it_cohp_k = 1:length(ar_w)                                
                fl_w_kstar_z = ar_w_kstar_z(it_cohp_k);
                fl_w_astar_z = ar_w_astar_z(it_cohp_k);
                
                % consumption
                fl_c = f_cons(fl_coh, fl_w_astar_z, fl_w_kstar_z);
                
                % current utility
                if (fl_crra == 1)
                    ar_val_cur(it_cohp_k) = f_util_log(fl_c);
                    fl_u_neg_c = f_util_log(fl_c_min);
                else
                    ar_val_cur(it_cohp_k) = f_util_crra(fl_c);
                    fl_u_neg_c = f_util_crra(fl_c_min);
                end
                
                % loop 4: add future utility, integration already done in
                % ff_wkz_evf
                fl_ev_condi_z_max_z = ar_ev_condi_z_max_z(it_cohp_k);
                ar_val_cur(it_cohp_k) = ar_val_cur(it_cohp_k) + fl_beta*fl_ev_condi_z_max_z;
                
                % Replace if negative consumption
                if fl_c <= 0
                    ar_val_cur(it_cohp_k) = fl_u_neg_c;
                end
                
            end
            
            % maximization over loop 3 choices for loop 1+2 states
            it_max_lin_idx = find(ar_val_cur == max(ar_val_cur));
            mt_val(it_coh_j,it_z_i) = ar_val_cur(it_max_lin_idx(1));
            mt_pol_a(it_coh_j,it_z_i) = ar_w_astar_z(it_max_lin_idx(1));
            mt_pol_k(it_coh_j,it_z_i) = ar_w_kstar_z(it_max_lin_idx(1));
            
        end
    end
    
    %% Check Tolerance and Continuation
    
    % Difference across iterations
    ar_val_diff_norm(it_iter) = norm(mt_val - mt_val_cur);
    ar_pol_diff_norm(it_iter) = norm(mt_pol_a - mt_pol_a_cur) + norm(mt_pol_k - mt_pol_k_cur);
    ar_pol_a_perc_change = sum((mt_pol_a ~= mt_pol_a_cur))/(it_ameshk_n);
    ar_pol_k_perc_change = sum((mt_pol_k ~= mt_pol_k_cur))/(it_ameshk_n);    
    mt_pol_perc_change(it_iter, :) = mean([ar_pol_a_perc_change;ar_pol_k_perc_change]);
    
    % Update
    mt_val_cur = mt_val;
    mt_pol_a_cur = mt_pol_a;
    mt_pol_k_cur = mt_pol_k;
    
    % Print Iteration Results
    if (bl_display && (rem(it_iter, it_display_every)==0))
        fprintf('VAL it_iter:%d, fl_diff:%d, fl_diff_pol:%d\n', ...
            it_iter, ar_val_diff_norm(it_iter), ar_pol_diff_norm(it_iter));
        tb_valpol_iter = array2table([mean(mt_val_cur,1);...
                                      mean(mt_pol_a_cur,1); ...
                                      mean(mt_pol_k_cur,1); ...
                                      mt_val_cur(it_ameshk_n,:); ...
                                      mt_pol_a_cur(it_ameshk_n,:); ...
                                      mt_pol_k_cur(it_ameshk_n,:)]);
        tb_valpol_iter.Properties.VariableNames = strcat('z', string((1:size(mt_val_cur,2))));
        tb_valpol_iter.Properties.RowNames = {'mval', 'map', 'mak', 'Hval', 'Hap', 'Hak'};
        disp('mval = mean(mt_val_cur,1), average value over a')
        disp('map  = mean(mt_pol_a_cur,1), average choice over a')
        disp('mkp  = mean(mt_pol_k_cur,1), average choice over k')
        disp('Hval = mt_val_cur(it_ameshk_n,:), highest a state val')
        disp('Hap = mt_pol_a_cur(it_ameshk_n,:), highest a state choice')
        disp('mak = mt_pol_k_cur(it_ameshk_n,:), highest k state choice')                
        disp(tb_valpol_iter);
    end
    
    % Continuation Conditions:
    % 1. if value function convergence criteria reached
    % 2. if policy function variation over iterations is less than
    % threshold
    if (it_iter == (it_maxiter_val + 1))
        bl_vfi_continue = false;
    elseif ((it_iter == it_maxiter_val) || ...
            (ar_val_diff_norm(it_iter) < fl_tol_val) || ...
            (sum(ar_pol_diff_norm(max(1, it_iter-it_tol_pol_nochange):it_iter)) < fl_tol_pol))
        % Fix to max, run again to save results if needed
        it_iter_last = it_iter;
        it_iter = it_maxiter_val;        
    end
    
end

% End Timer
if (bl_time)
    toc;
end

% End Profile
if (bl_profile)
    profile off
    profile viewer
    st_file_name = [st_profile_prefix st_profile_name_main st_profile_suffix];
    profsave(profile('info'), strcat(st_profile_path, st_file_name));
end

%% Process Optimal Choices

result_map = containers.Map('KeyType','char', 'ValueType','any');
result_map('mt_val') = mt_val;
result_map('mt_pol_a') = mt_pol_a;
result_map('mt_pol_k') = mt_pol_k;

if (bl_post)
    bl_input_override = true;
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);
    result_map = ff_akz_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
end

end
