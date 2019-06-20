%% Derive Asset and Choices/Outcomes Distribution (Vectorized)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_az_ds_vec(varargin)
%% FF_AZ_DS_VEC finds the stationary asset distributions Vectorized
% Building on the Asset Dynamic Programming Problem
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
% ff_az_vf_vecsv>, here we solve for the asset distribution using
% vectorized codes.
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html
% ff_az_ds> shows looped codes for finding asset distribution. The solution
% is the same. Both *ff_az_ds* and *ff_az_ds_vec* using
% optimized-vectorized dynamic programming code from ff_az_vf_vecsv. The
% idea here is that in addition to vectornizing the dynamic programming
% funcion, we can also vectorize the distribution code here. 
%
% Distributions of Interest:
%
% * $p(a,z)$
% * $p(Y=y, z) = \sum_{a} \left( 1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$
% * $p(Y=y, a) = \sum_{z} \left( 1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$
% * $p(Y=y) = \sum_{a,z} \left( 1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$
%
% Statistics include:
%
% * $\mu_y = \sum_{y} p(Y=y) \cdot y$
% * $\sigma_y = \sqrt{ \sum_{y} p(Y=y) \cdot \left( y - \mu_y \right)^2}$
% * $p(y=0)$
% * $p(y=\max(y))$
% * percentiles: $min_{y} \left\{ P(Y \le y) - percentile \mid P(Y \le y) \ge percentile \right\}$
% * fraction of outcome held by up to percentiles: $E(Y<y)/E(Y)$
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
% new keys included in result_map in addition to the output from
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
% ff_az_vf_vecsv> are various distribution statistics for each model
% outcome, keys include *cl_mt_pol_a*, *cl_mt_pol_c*, *cl_mt_pol_coh*, etc.
%
% @example
%
%    % Get Default Parameters
%    it_param_set = 6;
%    [param_map, support_map] = ffs_az_set_default_param(it_param_set);
%    % Change Keys in param_map
%    param_map('it_a_n') = 500;
%    param_map('it_z_n') = 11;
%    param_map('fl_a_max') = 100;
%    param_map('fl_w') = 1.3;
%    % Change Keys support_map
%    support_map('bl_display') = false;
%    support_map('bl_post') = true;
%    support_map('bl_display_final') = false;
%    % Call Program with external parameters that override defaults
%    ff_az_ds_vec(param_map, support_map);
%
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html ff_az_vf_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_ds_post_stats.html ff_az_ds_post_stats>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html fft_disc_rand_var_stats>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html fft_disc_rand_var_mass2outcomes>
%
% @seealso
%
% * derive distribution loop: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html ff_az_ds>
% * derive distribution vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vec.html ff_az_ds_vec>
%

%% Default
% Program can be externally invoked with _az_, _abz_ or various other
% programs. By default, program invokes using _az_ model programs:
%
% # it_subset = 5 is basic invoke quick test
% # it_subset = 6 is invoke full test
% # it_subset = 7 is profiling invoke
% # it_subset = 8 is matlab publish
% # it_subset = 9 is invoke operational (only final stats) and coh graph
%

params_len = length(varargin);
bl_input_override = 0;
if (params_len == 6)
    bl_input_override = varargin{6};
end

if (bl_input_override)
    % if invoked from outside override fully
    [param_map, support_map, armt_map, func_map, result_map, ~] = varargin{:};

else
    % default invoke
    close all;

    it_param_set = 8;
    bl_input_override = true;

    % 1. Generate Parameters
    [param_map, support_map] = ffs_az_set_default_param(it_param_set);

    % Note: param_map and support_map can be adjusted here or outside to override defaults
    % param_map('it_a_n') = 750;
    % param_map('it_z_n') = 15;

    % 2. Generate function and grids
    [armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override

    % 3. Solve value and policy function using az_vf_vecsv, if want to solve
    % other models, solve outside then provide result_map as input
    [result_map] = ff_az_vf_vecsv(param_map, support_map, armt_map, func_map);

end

%% Parse Parameters

% append function name
st_func_name = 'ff_az_ds_vec';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

% result_map
% ar_st_pol_names is from section _Process Optimal Choices_ in the value
% function code.
params_group = values(result_map, {'cl_mt_pol_a', 'mt_pol_idx', 'ar_st_pol_names'});
[cl_mt_pol_a, mt_pol_idx, ar_st_pol_names] = params_group{:};
mt_pol_a = deal(cl_mt_pol_a{1});

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans', 'ar_z'});
[ar_a, mt_z_trans, ar_z] = params_group{:};

% param_map
params_group = values(param_map, {'it_a_n', 'it_z_n'});
[it_a_n, it_z_n] = params_group{:};
params_group = values(param_map, {'it_maxiter_dist', 'fl_tol_dist'});
[it_maxiter_dist, fl_tol_dist] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile_dist', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display_dist', 'it_display_every', 'bl_display_final_dist', 'bl_post'});
[bl_profile_dist, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display_dist, it_display_every, bl_display_final_dist, bl_post] = params_group{:};

%% Start Profiler and Timer

% Start Profile
if (bl_profile_dist)
    close all;
    profile off;
    profile on;
end

% Start Timer
if (bl_time)
    tic;
end

%% *f(a,z)*: Initialize Output Matrixes
% Initialize the distribution to be uniform

mt_dist_az_init = ones(length(ar_a),length(ar_z))/length(ar_a)/length(ar_z);
mt_dist_az_cur = mt_dist_az_init;
mt_dist_az_zeros = zeros(length(ar_a),length(ar_z));

%% *f(a,z)*: Initialize Convergence Conditions

bl_histiter_continue = true;
it_iter = 0;
ar_dist_diff_norm = zeros([it_maxiter_dist, 1]);
mt_dist_perc_change = zeros([it_maxiter_dist, it_z_n]);

%% *f(a,z)*: Derive Stationary Distribution
% Iterate over the discrete joint random variable variables (a,z)
while (bl_histiter_continue)

    it_iter = it_iter + 1;

    %% *f(a,z)*: Vectorized Solution
    % this is the only part of the code that differs from
    % <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html
    % ff_az_ds.html> the looped code.
    
    % 1. initialize empty
    mt_dist_az = mt_dist_az_zeros;

    % 2. One loop remains
    for i = 1:it_z_n

        % 3. Get Unique Index (future states receive from multiple current states)
        [ar_idx_full, ~, ar_idx_of_unique] = unique(mt_pol_idx(:,i));
        mt_zi_prob = mt_dist_az_cur(:,i) * mt_z_trans(i,:);

        % 4. Cumulative probability received at state from zi
        [mt_idx_of_unique_mesh, mt_col_idx] = ndgrid(ar_idx_of_unique, 1:size(mt_zi_prob,2));
        mt_zi_cumu_prob = accumarray([mt_idx_of_unique_mesh(:) mt_col_idx(:)], mt_zi_prob(:));

        % 5. Adding up
        mt_dist_az(ar_idx_full, :) = mt_zi_cumu_prob +  mt_dist_az(ar_idx_full,:);
    end
    

    %% *f(a,z)*: Check Tolerance and Continuation

    % Difference across iterations
    ar_dist_diff_norm(it_iter) = norm(mt_dist_az - mt_dist_az_cur);
    mt_dist_perc_change(it_iter, :) = sum((mt_dist_az ~= mt_dist_az))/(it_a_n);

    % Update
    mt_dist_az_cur = mt_dist_az;

    % Print Iteration Results
    if (bl_display_dist && (rem(it_iter, it_display_every)==0))
        fprintf('Dist it_iter:%d, fl_dist_diff:%d\n', it_iter, ar_dist_diff_norm(it_iter));
        tb_hist_iter = array2table([sum(mt_dist_az_cur,1); std(mt_dist_az_cur,1); ...
                                    mt_dist_az_cur(1,:); mt_dist_az_cur(it_a_n,:)]);
        tb_hist_iter.Properties.VariableNames = strcat('z', string((1:size(mt_dist_az,2))));
        tb_hist_iter.Properties.RowNames = {'mdist','sddist', 'Ldist', 'Hdist'};
        disp('mdist = sum(mt_dist_az_cur,1) = sum_{a}(p(a)|z)')
        disp('sddist = std(mt_pol_a_cur,1) = std_{a}(p(a)|z)')
        disp('Ldist = mt_dist_az_cur(1,:) = p(min(a)|z)')
        disp('Hdist = mt_dist_az_cur(it_a_n,:) = p(max(a)|z)')
        disp(tb_hist_iter);
    end

    % Continuation Conditions:
    if (it_iter == (it_maxiter_dist + 1))
        bl_histiter_continue = false;
    elseif ((it_iter == it_maxiter_dist) || ...
            (ar_dist_diff_norm(it_iter) < fl_tol_dist))
        it_iter_last = it_iter;
        it_iter = it_maxiter_dist;
    end

end

%% End Time and Profiler

% End Timer
if (bl_time)
    toc;
end

% End Profile
if (bl_profile_dist)
    profile off
    profile viewer
    st_file_name = [st_profile_prefix st_profile_name_main st_profile_suffix];
    profsave(profile('info'), strcat(st_profile_path, st_file_name));
end

%% *f(y), f(c), f(a)*: Generate Key Distributional Statistics for Each outcome
% Having derived f(a,z) the probability mass function of the joint discrete
% random variables, we now obtain distributional statistics. Note that we
% know f(a,z), and we also know relevant policy functions a'(a,z), c(a,z),
% or other policy functions. We can simulate any choices that are a
% function of the random variables (a,z), using f(a,z)

bl_input_override = true;
result_map = ff_az_ds_post_stats(support_map, result_map, mt_dist_az, bl_input_override);

end
