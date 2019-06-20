%% Test Shock Grid Cout (Savings Dynamic Programming Problem)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Testing the
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
% ff_az_vf_vecsv> program for solving the savings only dynamic
% programming problem.
%
% Loop over the number of shock points. 9, 15, and 27 shocks points share
% the same AR1 shock process parameters following the benchmark simulation.
% The benchmark parameters can be adjusted below or inside the default
% parameters function:
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html
% ffs_az_set_default_param>.
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

%% Simulate Model with Different Shock Grid Count
close all

ar_it_z_n = [9, 15, 27];

for it_z_n = ar_it_z_n

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['it_z_n = ' num2str(it_z_n)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    it_param_set = 4;
    [param_map, support_map] = ffs_az_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_a_n') = 750;
    param_map('it_z_n') = it_z_n;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    % Call Program
    ff_az_vf_vecsv(param_map, support_map);

    % Snap
    snapnow;

end

% close all
close all;
clear all;