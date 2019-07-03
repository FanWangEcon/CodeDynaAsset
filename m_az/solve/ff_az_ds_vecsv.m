%% Derive Asset and Choices/Outcomes Distribution (Analytical)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_az_ds_vecsv(varargin)
%% FF_AZ_DS_VECSV finds the stationary asset distributions analytically
% Here, we implement the iteration free semi-analytical method for finding
% asset distributions. The method analytically give the exact
% stationary distribution induced by the policy function from the dynamic
% programming problem, conditional on discretizations.
%
% See the appedix of
% <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3316939 Wang (2019)>
% which develops details on how this works. Suppose endo state size is N and
% shock is size M, then our transition matrix is (NxM) by (NxM). We know
% the coh(a,z) value associated with each element of the (NxM) by 1 array.
% We also know f(a'(a,z),z'|z) transition probability. We contruct a markov
% chain that has (NxM) states. Specifically:
%
% * We need to transform: mt_pol_idx. This matrix is indexing 1 through N,
% we need for it to index 1 through (NxM).
% * Then we need to duplicate the transition matrix fro shocks.
% * Transition Matrix is *sparse*
%
% Once we have the all states meshed markov transition matrix, then we can
% use standard methods to find the stationary distribution. Three options
% are offered here that provide identical solutions:
%
% # The Eigenvector Approach: very fast
% # The Projection Approach: medium
% # The Power Approach: very slow (especially with sparse matrix)
%
% The program here builds on the Asset Dynamic Programming Problem
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
%    ff_az_ds_vecsv(param_map, support_map);
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
% * derive distribution semi-analytical: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html ff_az_ds_vecsv>
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

    param_map('st_analytical_stationary_type') = 'eigenvector';
%     param_map('st_analytical_stationary_type') = 'projection';
%     param_map('st_analytical_stationary_type') = 'power';
    
    % 2. Generate function and grids
    [armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override

    % 3. Solve value and policy function using az_vf_vecsv, if want to solve
    % other models, solve outside then provide result_map as input
    [result_map] = ff_az_vf_vecsv(param_map, support_map, armt_map, func_map);

end

%% Parse Parameters

% append function name
st_func_name = 'ff_az_ds_vecsv';
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
params_group = values(armt_map, {'mt_z_trans'});
[mt_z_trans] = params_group{:};

% param_map
params_group = values(param_map, { 'it_trans_power_dist', 'st_analytical_stationary_type'});
[it_trans_power_dist, st_analytical_stationary_type] = params_group{:};


% support_map
params_group = values(support_map, {'bl_profile_dist', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display_final_dist', 'bl_post'});
[bl_profile_dist, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display_final_dist, bl_post] = params_group{:};

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

%% Get Size of Endogenous and Exogenous State
% The key idea is that all information for policy function is captured by
% _mt_pol_idx_ matrix, its rows are the number of endogenous states, and
% its columns are the exogenous shocks.

[it_endostates_rows_n, it_exostates_cols_n] = size(mt_pol_idx);

%% 1. Generate Max Index in (NxM) from (N) array.
% Suppose we have: 
%
%   mt_pol_idx =
% 
%      1     1     2
%      2     2     2
%      3     3     3
%      4     4     4
%      4     5     5
%
% These become:
%
%   mt_pol_idx_mesh_max = 
% 
%      1     6    11
%      2     7    12
%      3     8    13
%      4     9    14
%      5    10    15
%      1     6    11
%      2     7    12
%      3     8    13
%      4     9    14
%      5    10    15
%      1     6    11
%      2     7    12
%      3     8    13
%      4     9    14
%      5    10    15
%

% mt_pol_idx_mesh_max is (NxM) by M, mt_pol_idx is N by M
mt_pol_idx_mesh_max = mt_pol_idx(:) + (0:1:(it_exostates_cols_n-1))*it_endostates_rows_n;

%% 2. Transition Probabilities from (M by M) to (NxM) by M
%
%   mt_trans_prob =
% 
%     0.9332    0.0668    0.0000
%     0.9332    0.0668    0.0000
%     0.9332    0.0668    0.0000
%     0.9332    0.0668    0.0000
%     0.9332    0.0668    0.0000
%     0.0062    0.9876    0.0062
%     0.0062    0.9876    0.0062
%     0.0062    0.9876    0.0062
%     0.0062    0.9876    0.0062
%     0.0062    0.9876    0.0062
%     0.0000    0.0668    0.9332
%     0.0000    0.0668    0.9332
%     0.0000    0.0668    0.9332
%     0.0000    0.0668    0.9332
%     0.0000    0.0668    0.9332
%

mt_trans_prob = reshape(repmat(mt_z_trans(:)', [it_endostates_rows_n, 1]), [it_endostates_rows_n*it_exostates_cols_n, it_exostates_cols_n]);

%% 3. Fill mt_pol_idx_mesh_idx to mt_full_trans_mat
% Try to always use sparse matrix, unless grid sizes very small, keeping
% non-sparse code here for comparison. Sparse matrix is important for
% allowing the code to be fast and memory efficient. Otherwise this method
% is much slower than iterative method.

it_sparse_threshold = 100*7;

if (it_endostates_rows_n*it_exostates_cols_n > it_sparse_threshold)
    
    %% 3.1 Sparse Matrix Approach
    i = mt_pol_idx_mesh_max(:);
    j = repmat((1:1:it_endostates_rows_n*it_exostates_cols_n),[1,it_exostates_cols_n])';
    v = mt_trans_prob(:);
    m = it_endostates_rows_n*it_exostates_cols_n;
    n = it_endostates_rows_n*it_exostates_cols_n;
    mt_full_trans_mat = sparse(i, j, v, m, n);
    
else
    
    %% 3.2 Full Matrix Approach
    %
    %   ar_lin_idx_start_point =
    % 
    %      0    15    30    45    60    75    90   105   120   135   150   165   180   195   210
    %
    %   mt_pol_idx_mesh_idx_meshfull =
    % 
    %      1     6    11
    %     17    22    27
    %     33    38    43
    %     49    54    59
    %     65    70    75
    %     76    81    86
    %     92    97   102
    %    108   113   118
    %    124   129   134
    %    140   145   150
    %    151   156   161
    %    167   172   177
    %    183   188   193
    %    199   204   209
    %    215   220   225
    %

    % Each row's linear index starting point
    ar_lin_idx_start_point = ((it_endostates_rows_n*it_exostates_cols_n)*(0:1:(it_endostates_rows_n*it_exostates_cols_n-1)));

    % mt_pol_idx_mesh_idx_meshfull is (NxM) by M
    % Full index in (NxM) to (NxM) transition Matrix
    mt_pol_idx_mesh_idx_meshfull = mt_pol_idx_mesh_max + ar_lin_idx_start_point';

    % Fill mt_pol_idx_mesh_idx to mt_full_trans_mat
    mt_full_trans_mat = zeros([it_endostates_rows_n*it_exostates_cols_n, it_endostates_rows_n*it_exostates_cols_n]);    
    mt_full_trans_mat(mt_pol_idx_mesh_idx_meshfull(:)) = mt_trans_prob(:);
    
end

%% 4. Stationary Distribution *Method A*, Eigenvector Approach
% Given that markov chain we have constructured for all state-space
% elements, we can now find the stationary distribution using standard
% <https://en.wikipedia.org/wiki/Markov_chain#Stationary_distribution_relation_to_eigenvectors_and_simplices
% eigenvector> approach.

if (strcmp(st_analytical_stationary_type, 'eigenvector'))
    [V, ~] = eigs(mt_full_trans_mat,1,1);
    ar_stationary = V/sum(V);
end

%% 5. Stationary Distribution *Method B*, Projection
% This uses Projection. 

if (strcmp(st_analytical_stationary_type, 'projection'))
    
    % a. Transition - I
    % P = P*T
    % 0 = P*T - P
    % 0 = P*(T-1)                
    % Q = trans_prob - np.identity(state_count);

    mt_diag = eye(it_endostates_rows_n*it_exostates_cols_n);
    if (it_endostates_rows_n*it_exostates_cols_n > it_sparse_threshold)
        % if larger, use sparse matrix
        mt_diag = sparse(mt_diag);
    end
    mt_Q = mt_full_trans_mat' - mt_diag;

    % b. add all 1 as final column, (because P*1 = 1) 
    % one_col = np.ones((state_count,1))
    % Q = np.column_stack((Q, one_col))

    ar_one = ones([it_endostates_rows_n*it_exostates_cols_n,1]);
    mt_Q = [mt_Q, ar_one];

    % c. b is the LHS 
    % b = [0,0,0,...,1]
    % b = np.zeros((1, (state_count+1)))
    % b[0, state_count] = 1

    ar_b = zeros([1, it_endostates_rows_n*it_exostates_cols_n+1]);
    if (it_endostates_rows_n*it_exostates_cols_n > it_sparse_threshold)
        % if larger, use sparse matrix
        ar_b = sparse(ar_b);
    end    
    ar_b(it_endostates_rows_n*it_exostates_cols_n+1) = 1;

    % d. solve
    % b = P*Q
    % b*Q^{T} = P*Q*Q^{T}
    % P*Q*Q^{T} = b*Q^{T}
    % P = (b*Q^{T})[(Q*Q^{T})^{-1}]
    % Q_t = np.transpose(Q)
    % b_QT = np.dot(b, Q_t)
    % Q_QT = np.dot(Q, Q_t)

    % inv_mt_Q_QT = inv(mt_Q*mt_Q');
    ar_stationary = (ar_b*mt_Q')/(mt_Q*mt_Q');

end

%% 6. Stationary Distribution *Method C*, Power
% Takes markov chain to Nth power. This is the slowest.

if (strcmp(st_analytical_stationary_type, 'power'))
    
    mt_stationary_full = (mt_full_trans_mat)^it_trans_power_dist;
    ar_stationary = mt_stationary_full(:,1);
end

%% 7. Stationary Vector to Stationary Matrix in Original Dimensions

mt_dist_az = reshape(ar_stationary, size(mt_pol_idx));

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
% function of the random variables (a,z), using f(a,z). We call function
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_ds_post_stats.html
% ff_az_ds_post_stats> which uses
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html
% fft_disc_rand_var_stats> and
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html
% fft_disc_rand_var_mass2outcomes> to compute various statistics of
% interest.

bl_input_override = true;
result_map = ff_az_ds_post_stats(support_map, result_map, mt_dist_az, bl_input_override);

end
