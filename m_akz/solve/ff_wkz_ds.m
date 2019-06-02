function [ar_dist_a, mt_dist_az, mt_pol, flag] = ds_az(param_map, support_map)

%% Solve Value/Pol Function
[~, mt_pol, ~] = vf_az(param_map, support_map);

%% Process Paramaters
% Shock Parameters
it_z_n = param_map('it_z_n');
ar_z =  param_map('ar_z');
mt_z_trans = param_map('mt_z_trans');

% Choice vector
it_a_n = param_map('it_a_n');
ar_a = param_map('ar_a');

% History Accuracy
it_maxiter_dist = param_map('it_maxiter_dist');
fl_tol_dist = param_map('fl_tol_dist');

% Preferences
fl_crra = param_map('fl_crra');
fl_rho = param_map('fl_rho');
fl_sig = param_map('fl_sig');

%% Simulating History
% Initialize distribution
mt_D0 = ones(length(ar_a),length(ar_z))/length(ar_a)/length(ar_z);

it_iter = 1;
fl_diff = 1;

while it_iter <= it_maxiter_dist && fl_diff > fl_tol_dist
    mt_D1 = zeros(length(ar_a),length(ar_z));
    for i = 1:length(ar_z)
        for j = 1:length(ar_a)
            aprime = mt_pol(j,i);
            loc_a_prime = find(ar_a == aprime);
            pi = mt_D0(j,i)*mt_z_trans(i,:);
            mt_D1(loc_a_prime,:) = mt_D1(loc_a_prime,:) + pi;
        end
    end
    fl_diff = norm(mt_D1 - mt_D0);
    it_iter = it_iter + 1;
    mt_D0 = mt_D1;
end

%% Output
ar_dist_a = sum(mt_D1,2);
mt_dist_az = mt_D1;

%% Flag
if fl_diff <= fl_tol_dist && it_iter > it_maxiter_dist
%     tolerance did not reach, iteration bound reached
    flag = 1;
elseif fl_diff > fl_tol_dist && it_iter <= it_maxiter_dist
%     tolerance reached, iter bound did not reach
    flag = 2;
else
%     both reached
    flag = 0;
end

%% Save mat
if (support_map('bl_save_mat_dist'))
    st_file_name = ['ds_' support_map('st_file_name')];
    save(ff_sup_save_prep(support_map('st_path_folder'), st_file_name));
end

end
