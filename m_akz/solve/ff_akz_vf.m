%% Risky + Safe Asset Dyna Prog Concurrent Solution (Loop)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository> 
% Table of Content.*

%%
function result_map = ff_akz_vf(varargin)
%% FF_AKZ_VF solve infinite horizon exo shock + endo asset problem
% This program solves the infinite horizon dynamic savings and risky
% capital asset problem with some ar1 shock. The two assets could be a safe
% bond and a risky stock if the risky asset has constant return to scale.
% Alternatively, the risky asset could be a capital investment asset with
% decreasing return to scale. There is risky component to the risky capital
% investment, but one could also potentially resale a fraction of the risky
% capital after depreciation, giving the household/investor/entrepreur a
% safe minimum return to risky investment. The state variable is the
% cash-on-hand, which is determined by risky and safe asset as well as the
% shock jointly.
%
% This problem is more computationally intensive to solve
% compared to the single asset problem, shown here:
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf.html
% ff_az_vf>. Here I show the looped solution. As before, as we expand the 
% state and choices spaces to have a looped slow version of the code so
% that one could debug and make sure the model works. 
%
% Vectorized, optimized vectorized versions of the code are shown here:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vec.html
% ff_akz_vf_vec> and here:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vecsv.html
% ff_akz_vf_vecsv>. These improve upon the speed, yet the
% speed is still not as fast as we need for some applications. 
%
% The basic forms of the one and two asset problems have analytically
% tractable mathmatical solutions. The benefit of these grid based solution
% algorithm is that we are able to add very flexible constraints to the
% problem or incorporate additional discrete choices that make the
% underlying problem non-continous and non-differentiable and mathamtically
% untractable.
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
% * mt_pol_a matrix states_n by shock_n matrix of converged policy function
% grid safe asset
% * mt_pol_k matrix states_n by shock_n matrix of converged policy function
% grid risky asset
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
%   it_param_set = 2;
%   [param_map, support_map] = ffs_akz_set_default_param(it_param_set);   
%   % Simulation Accuracy
%   param_map('it_w_n') = 750;
%   param_map('it_ak_n') = param_map('it_w_n');
%   param_map('it_z_n') = 11;
%   % Display Parameters
%   support_map('bl_display') = false;
%   support_map('bl_display_final') = false;
%   support_map('bl_time') = true;
%   support_map('bl_profile') = false;
%   % Call Program with external parameters that override defaults
%   ff_akz_vf(param_map, support_map);
% 
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html ffs_akz_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html ffs_akz_get_funcgrid>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solvepost/html/ff_akz_vf_post.html ff_akz_vf_post>
%
% @seealso
%
% * concurrent (safe + risky) loop: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf.html ff_akz_vf>
% * concurrent (safe + risky) vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vec.html ff_akz_vf_vec>
% * concurrent (safe + risky) optimized-vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vecsv.html ff_akz_vf_vecsv>
% * two-stage (safe + risky) loop: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_wkz_vf.html ff_wkz_vf>
% * two-stage (safe + risky) vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_wkz_vf_vec.html ff_wkz_vf_vec>
% * two-stage (safe + risky) optimized-vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_wkz_vf_vecsv.html ff_wkz_vf_vecsv>
% * two-stage + interpolate (safe + risky) loop: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_vf.html ff_iwkz_vf>
% * two-stage + interpolate (safe + risky) vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_vf_vec.html ff_iwkz_vf_vec>
% * two-stage + interpolate (safe + risky) optimized-vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_vf_vecsv.html ff_iwkz_vf_vecsv>
%

%% Default
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button

it_param_set = 4;
[param_map, support_map] = ffs_akz_set_default_param(it_param_set);

% Note: param_map and support_map can be adjusted here or outside to override defaults
% param_map('it_w_n') = 50;
% param_map('it_z_n') = 15;

% get armt and func map
[armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map); % 1 for override
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
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% append function name
st_func_name = 'ff_akz_vf';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'mt_z_trans', 'ar_z'});
[ mt_z_trans, ar_z] = params_group{:};
params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'mt_coh_wkb', 'it_ameshk_n'});
[ar_a_meshk, ar_k_mesha, mt_coh_wkb, it_ameshk_n] = params_group{:};
% func_map
params_group = values(func_map, {'f_util_log', 'f_util_crra', 'f_cons', 'f_coh'});
[f_util_log, f_util_crra, f_cons, f_coh] = params_group{:};
% param_map
params_group = values(param_map, {'it_z_n', 'fl_crra', 'fl_beta', 'fl_c_min'});
[it_z_n, fl_crra, fl_beta, fl_c_min] = params_group{:};
params_group = values(param_map, {'it_maxiter_val', 'fl_tol_val', 'fl_tol_pol', 'it_tol_pol_nochange'});
[it_maxiter_val, fl_tol_val, fl_tol_pol, it_tol_pol_nochange] = params_group{:};
% support_map
params_group = values(support_map, {'bl_profile', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display_defparam', 'bl_display', 'it_display_every', 'bl_post'});
[bl_profile, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display_defparam, bl_display, it_display_every, bl_post] = params_group{:};
params_group = values(support_map, {'it_display_summmat_rowmax', 'it_display_summmat_colmax'});
[it_display_summmat_rowmax, it_display_summmat_colmax] = params_group{:};

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
    
    %% Solve Optimization Problem Current Iteration
    
    % loop 1: over exogenous states
    for it_z_i = 1:length(ar_z)
        
        % loop 2: over endogenous states
        for it_coh_j = 1:length(ar_a_meshk)
            % Get cash-on-hand which include k,b,z
            fl_coh = mt_coh_wkb(it_coh_j, it_z_i);
            
            % loop 3: over choices
            ar_val_cur = zeros(size(ar_a_meshk));
            for it_cohp_k = 1:length(ar_a_meshk)
                fl_ap = ar_a_meshk(it_cohp_k);
                fl_kp = ar_k_mesha(it_cohp_k);
                
                % consumption
                fl_c = f_cons(fl_coh, fl_ap, fl_kp);
                
                % current utility
                if (fl_crra == 1)
                    ar_val_cur(it_cohp_k) = f_util_log(fl_c);
                    fl_u_neg_c = f_util_log(fl_c_min);
                else
                    ar_val_cur(it_cohp_k) = f_util_crra(fl_c);
                    fl_u_neg_c = f_util_crra(fl_c_min);
                end
                
                % loop 4: add future utility, integration--loop over future shocks
                for it_zp_q = 1:length(ar_z)
                    ar_val_cur(it_cohp_k) = ar_val_cur(it_cohp_k) + fl_beta*mt_z_trans(it_z_i,it_zp_q)*mt_val_cur(it_cohp_k,it_zp_q);
                end
                
                % Replace if negative consumption
                if fl_c <= 0
                    ar_val_cur(it_cohp_k) = fl_u_neg_c;
                end
                
            end
            
            % maximization over loop 3 choices for loop 1+2 states
            it_max_lin_idx = find(ar_val_cur == max(ar_val_cur));
            mt_val(it_coh_j,it_z_i) = ar_val_cur(it_max_lin_idx(1));
            mt_pol_a(it_coh_j,it_z_i) = ar_a_meshk(it_max_lin_idx(1));
            mt_pol_k(it_coh_j,it_z_i) = ar_k_mesha(it_max_lin_idx(1));
            
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

mt_coh = f_coh(ar_z, ar_a_meshk, ar_k_mesha);
result_map('cl_mt_coh') = {mt_coh, zeros(1)};
result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
result_map('cl_mt_pol_k') = {mt_pol_k, zeros(1)};
result_map('cl_mt_pol_c') = {f_cons(mt_coh, mt_pol_a, mt_pol_k), zeros(1)};
result_map('ar_st_pol_names') = ["cl_mt_coh", "cl_mt_pol_a", "cl_mt_pol_k", "cl_mt_pol_c"];

if (bl_post)
    bl_input_override = true;
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);
    result_map = ff_akz_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
end

%% Display Various Containers

if (bl_display_defparam)

    %% Display 1 support_map
    fft_container_map_display(support_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 2 armt_map
    fft_container_map_display(armt_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 3 param_map
    fft_container_map_display(param_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 4 func_map
    fft_container_map_display(func_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 5 result_map
    fft_container_map_display(result_map, it_display_summmat_rowmax, it_display_summmat_colmax);

end

end
