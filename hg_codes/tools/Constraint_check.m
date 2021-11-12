function [Constraints, Constraint_uppers, Constraint_Max, Constraint_Min]  = Constraint_check(Internal_Stresses, Safty_Factor, Yield_strength)


% skin checks (yield)
Skin_yield_check=Safty_Factor*Internal_Stresses.skin_VonMise/Yield_strength;


% skin checks (buckling)
Skin_buckling_check1=Safty_Factor*Internal_Stresses.skin_principal_normal./Internal_Stresses.skin_critical_buck_normal;

Skin_buckling_check2=Safty_Factor*Internal_Stresses.skin_principal_shear./Internal_Stresses.skin_critical_buck_shear;


% spar cap checks (yield)
Spar_cap_yield_check=Safty_Factor*Internal_Stresses.sparcap_VonMise/Yield_strength;

% spar cap checks (column buckling)
Spar_cap_CB_check=Safty_Factor*Internal_Stresses.sparcap_VonMise./Internal_Stresses.sparcap_ColumnBuckling;

% spar cap checks (sheet buckling)
Spar_cap_SheetB_check=Safty_Factor*Internal_Stresses.sparcap_VonMise./Internal_Stresses.sparcap_SheetBuckling;


% spar web checks (yield)
Spar_web_yield_check=Safty_Factor*Internal_Stresses.sparweb_VonMise/Yield_strength;


% spar web checks (shear buckling)
Spar_web_buckling_check=Safty_Factor*Internal_Stresses.sparweb_Shear./Internal_Stresses.sparweb_ShearBuckling;


% stringer checks (crippling)
Strg_crippling_check=Safty_Factor*Internal_Stresses.Strg_compression./Internal_Stresses.Strg_crippling;


% stringer checks (column buckling)
Strg_CB_check=Safty_Factor*Internal_Stresses.Strg_compression./Internal_Stresses.Strg_column_buckle;


%% Result output

% all constraints
Constraints.Skin_yield_check=Skin_yield_check;

Constraints.Skin_buckling_check1=Skin_buckling_check1;

Constraints.Skin_buckling_check2=Skin_buckling_check2;


Constraints.Spar_cap_yield_check=Spar_cap_yield_check;

Constraints.Spar_cap_CB_check=Spar_cap_CB_check;

Constraints.Spar_cap_SheetBuck_check=Spar_cap_SheetB_check;


Constraints.Spar_web_yield_check=Spar_web_yield_check;

Constraints.Spar_web_buckling_check=Spar_web_buckling_check;


Constraints.Strg_crippling_check=Strg_crippling_check;

Constraints.Strg_CB_check=Strg_CB_check;


% upper values for thickness adjusting 

% if skin failed by normal stress buckling, add spar thickness 
% only added skin thickness if it is failed by shear buckling, in that
% case, spar cannot help anyway. 

Constraint_uppers.Skin_Constraints_Upper=Skin_buckling_check2;

Constraint_uppers.Spar_cap_Constraints_Upper=max([Skin_buckling_check1; Spar_cap_yield_check; Spar_cap_CB_check; Constraints.Spar_cap_SheetBuck_check]);

Constraint_uppers.Spar_web_Constraints_Upper=max([Spar_web_yield_check;Spar_web_buckling_check]);

Constraint_uppers.Strg_Constraints_Upper=max([Strg_crippling_check; Strg_CB_check]);


% end loop condition
Constraint_Max=max([Constraint_uppers.Skin_Constraints_Upper, Constraint_uppers.Spar_web_Constraints_Upper, Constraint_uppers.Strg_Constraints_Upper]);

Constraint_Min=min([Constraint_uppers.Skin_Constraints_Upper, Constraint_uppers.Spar_web_Constraints_Upper, Constraint_uppers.Strg_Constraints_Upper]);


end