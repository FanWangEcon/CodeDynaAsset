%% Set Model Functions (ABZ FIBS)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [f_util_log, f_util_crra, f_util_standin, f_inc, f_coh, f_cons_coh_fbis, f_cons_coh_save, f_bprime] = ...
    ffs_abz_fibs_set_functions(varargin)
%% FFS_ABZ_SET_FIBS_FUNCTIONS setting model functions
% Define functions for the abz problem with formal and informal borrowing
% and savings. See here:
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_set_functions.html
% ffs_abz_set_functions> for *abz* function settings without without formal
% and informal considerations.
%
% @param fl_crra float crra utility
%
% @param fl_c_min float minimum consumption
%
% @param fl_r_fbr float formal borrowing interest rate
%
% @param fl_r_fsv float savings interest rates
%
% @param fl_w float wage rate
%
% @return f_util_log handle log utility
%
% @return f_util_crra handle crra utility general
%
% @return f_util_standin handle log utility with coh for testing graphing codes
%
% @return f_inc income handle equation wage and interests
%
% @return f_coh_fbis handle cash on hand equation given current period
% shock and formal and informal borrowing and savings choices.
%
% @return f_cons_coh_fbis handle consumption given cash on hand and formal
% and informal borrowing and savings choices.
%
% @return f_bprime handle aggregate borrowing and savings principles
%
% @example
%
%   [f_util_log, f_util_crra, f_util_standin, f_inc, f_coh_fbis, f_cons_coh_fbis] = ...
%        ffs_abz_set_functions(fl_crra, fl_c_min, fl_r_fbr, fl_r_fsv, fl_w);
%

%% Default

[fl_crra, fl_c_min] = deal(1, 0.001);
[fl_r_fbr, fl_r_fsv, fl_w] = deal(0.035, 0.025, 1.28);
bl_display_setfunc = false;
bl_b_is_principle = false;
default_params = {fl_crra fl_c_min fl_r_fbr fl_r_fsv fl_w bl_display_setfunc};

%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[fl_crra, fl_c_min, fl_r_fbr, fl_r_fsv, fl_w, bl_display_setfunc] = default_params{:};

%% Define Utility Functions

f_util_log = @(c) log(c);
f_util_crra = @(c) (((c).^(1-fl_crra)-1)./(1-fl_crra));

%% Define Income Equation

if (bl_b_is_principle)
    f_inc = @(ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save) ...
        (ar_z*fl_w ...
        + (ar_for_borr.*(fl_r_fbr) ...
        + ar_inf_borr.*(fl_r_inf) ...
        + ar_for_save.*(fl_r_fsv)));
else
    f_inc = @(ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save) ...
        (ar_z*fl_w ...
        + ((ar_for_borr./(1+fl_r_fbr))*fl_r_fbr ...
        + (ar_inf_borr./(1+fl_r_inf))*fl_r_inf ...
        + (ar_for_save./(1+fl_r_fsv))*fl_r_fsv));
end

if (bl_display_setfunc)
    disp('f_inc')
    disp(f_inc)
    [ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save] = deal(1, 0.06, -1, -0.33, 0);
    fl_inc = f_inc(ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save);
    fprintf('fl_inc=%.3f\nwith fl_r_inf:%.3f, ar_for_borr:%.3f, ar_inf_borr:%.3f, ar_for_save:%.3f',...
        fl_inc, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save);
end

%% Define Cash on Hand Current Period

if (bl_b_is_principle)
    f_coh = @(ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save) ...
        (ar_z*fl_w ...
        + (ar_for_borr.*(1+fl_r_fbr) ...
        + ar_inf_borr.*(1+fl_r_inf) ...
        + ar_for_save.*(1+fl_r_fsv)));
else
    f_coh = @(ar_z, ar_b) (ar_z*fl_w + ar_b);
end

if (bl_display_setfunc)
    disp('f_coh_fbis')
    disp(f_coh)
    if (bl_b_is_principle)
        [ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save] = deal(1, 0.06, -1, -0.33, 0);
        fl_coh_fbis = f_coh(ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save);
        fprintf('fl_coh_fbis=%.3f\nwith fl_r_inf:%.3f, ar_for_borr:%.3f, ar_inf_borr:%.3f, ar_for_save:%.3f',...
            fl_coh_fbis, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save);
    end
end


%% Define Consumption

if (bl_b_is_principle)

    % ar_for_borr, ar_inf_borr, ar_for_save are choices
    f_cons_coh_fbis = @(coh, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save) ...
        (coh + (ar_for_borr + ar_inf_borr + ar_for_save));

    f_cons_coh_save = @(coh, ar_for_save) (coh - ar_for_save);

    f_bprime = @(fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save) ...
        (ar_for_borr + ar_inf_borr + ar_for_save);

else
    % ar_bprime_in_c: array of b prime choice in current consumption units
    f_cons_coh_fbis = @(coh, ar_bprime_in_c) (coh + ar_bprime_in_c);

    f_cons_coh_save = @(coh, ar_for_save) (coh - ar_for_save./(1+fl_r_fsv));

    f_bprime = @(fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save) ...
        (ar_for_borr./(1+fl_r_fbr) ...
        + ar_inf_borr./(1+fl_r_inf) ...
        + ar_for_save./(1+fl_r_fsv));

end

if (bl_display_setfunc)
    disp('f_cons_coh_fbis')
    disp(f_cons_coh_fbis)

    if (bl_b_is_principle)
        % current period states
        [ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save] = deal(1, 0.06, -1, -0.33, 0);
        fl_coh_fbis = f_coh(ar_z, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save);

        % current period choices
        [ar_for_borr, ar_inf_borr, ar_for_save] = deal(-1, -0.33, 0);
        fl_cons_coh_fbis = f_cons_coh_fbis(fl_coh_fbis, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save);

        % display
        fprintf('fl_cons_coh_fbis=%.3f\nwith fl_r_inf:%.3f, ar_for_borr:%.3f, ar_inf_borr:%.3f, ar_for_save:%.3f',...
            fl_cons_coh_fbis, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save);
    end

end

%% Defome Stand-in Utility for Testing
f_coh_simple = @(ar_z, ar_b) (ar_z*fl_w + ar_b);
f_util_standin = @(z, b) f_util_log(f_coh_simple(z, b).*(f_coh_simple(z, b) > 0) + fl_c_min.*(f_coh_simple(z, b) <= 0));

end
