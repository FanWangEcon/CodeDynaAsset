%% Test Asset Grid Count (Save + Borr Dynamic Programming Problem)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Testing the
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
% program for solving the savings and borrowing dynamic programming problem.
%
% Loop over the number of asset points. 100 points means 100 grid points
% for the choice grid as well as the asset grid. The choice state grid is
% the same as the asset state grid.
%
% Here we see how does changing the asset/choice grid count from 100 to 750
% and to 2250 points change the model solutions.
%
% @seealso
%
% * _PRECISION_: savings only quick vs benchmark testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_vf_vecsv/test_precision/html/fsi_az_vf_vecsv_main.html fsi_az_vf_vecsv_main>
% * _PRECISION_: savings only asset grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_vf_vecsv/test_precision/html/fsi_az_vf_vecsv_a_n.html fsi_az_vf_vecsv_a_n>
% * _PRECISION_: savings only shock grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_vf_vecsv/test_precision/html/fsi_az_vf_vecsv_z_n.html fsi_az_vf_vecsv_z_n>
% * _BORROW GRID_: borrowing choice grid with default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ffs_abz_get_funcgrid/test_borr/html/ffs_abz_get_funcgrid_defnodfalt.html ffs_abz_get_funcgrid_defnodfalt>
% * _BORROW_: borrow and default small grid testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_borr/html/ff_abz_vf_vecsv_default_small.html ff_abz_vf_vecsv_default_small>
% * _BORROW_: borrow and default large grid testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_borr/html/ff_abz_vf_vecsv_default_large.html ff_abz_vf_vecsv_default_large>
% * _PRECISION_: borr + save quick vs benchmark testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_precision/html/fsi_abz_vf_vecsv_main.html fsi_abz_vf_vecsv_main>
% * _PRECISION_: borr + save asset grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_precision/html/fsi_abz_vf_vecsv_a_n.html fsi_abz_vf_vecsv_a_n>
% * _PRECISION_: borr + save shock grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_precision/html/fsi_abz_vf_vecsv_z_n.html fsi_abz_vf_vecsv_z_n>
%

%% Simulate Model with Different Asset Grid Count
close all

ar_it_a_n = [100, 750, 2250];
% ar_it_w_n = [25, 50];

for it_a_n = ar_it_a_n

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['it_a_n = ' num2str(it_a_n)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    it_param_set = 4;
    [param_map, support_map] = ffs_abz_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_a_n') = it_a_n;
    param_map('it_z_n') = 15;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    % Call Program
    ff_abz_vf_vecsv(param_map, support_map);
    % Snap
    snapnow;

end

% close all
close all;
clear all;
