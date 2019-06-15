    %% Compute Consumption given Borrowing and Savings
    if (fl_ap < 0)                    

        % Calculate Bridge Loan Borrowing
        if (fl_coh < 0)
            bl_input_override = true;
            [fl_aprime_nobridge, fl_c_bridge] = ffs_fibs_inf_bridge(...
                bl_b_is_principle, fl_r_inf, fl_ap, fl_coh, ...
                bl_display_infbridge, bl_input_override);

        else

            fl_aprime_nobridge = fl_ap;
            fl_c_bridge = 0;

        end

        % Find Optimal Formal Informal Borrow Save Combo
        % calculate consumption gain from formal + informal
        % borrowing and savings choices. 
        bl_input_override = true;
        fl_max_c_nobridge = ffs_fibs_min_c_cost(...
            bl_b_is_principle, fl_r_inf, fl_r_fsv, ...
            ar_forbrblk_r, ar_forbrblk, ...
            fl_aprime_nobridge, bl_display_minccost, bl_input_override);

        % Compute Consumption given Formal and Informal joint
        % consumption with formal borrow menu + bridge loans.
        fl_c = f_cons_coh_fbis(fl_coh, fl_max_c_nobridge + fl_c_bridge);

    else

        % consumption with savings                    
        fl_c = f_cons_coh_save(fl_coh, fl_ap);

    end
