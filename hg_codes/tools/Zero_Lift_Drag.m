
function CD0_All = Zero_Lift_Drag(Wing_Structural_Data,Wing_Area)

    
   Bwidth = Wing_Structural_Data.Wingbox_property.Bwidth;
   Bheight = Wing_Structural_Data.Wingbox_property.Bheight;


    Re_c=1e6; % Ref. Reynold's number
    % Re=2.56e7; % Reynold's number
    mu=1.47e-5; % air viscosity
    c=4.18; % mean aerodynamic chord
    rho=0.4; % air density
    U=230; % free stream velocity
    Ma=0.78; % mach number
    x_tc=0.4; % relative chord position for the max. thickness
    Q=1;    % Interface factor: 1 for wing
    tc=Bheight./(2*Bwidth); % thickness to chord ratio along the span
    AreaR=2*(Bheight+2*Bwidth)./(2*Bwidth); % wetted area to reference area (planform area)
    Roughness=1.33e-5;
    Dyn_Pressure=0.5*rho*U^2;
    
    % Cut off Re number
    Xc=mu*Re_c/(rho*c*U);
    Re_wingcut=38.21*(c/Roughness)^1.053;
    
    % Reynold's number
    Re_wing_=rho*U*c/mu;
    Re_wing=min(Re_wing_,Re_wingcut);
    
    % drag coefficient calculation
    C_lam=1.328/sqrt(Re_wing);
    C_turb=0.455/(log10(Re_wing)^2.58*(1+0.144*Ma^2)^0.65);
    Cf=Xc*C_lam+(1-Xc)*C_turb;
    
    % Form factor to include the effect of pressure drag
    FF=(1+(0.6/x_tc)*tc+100*tc.^4).*(1.34*Ma^0.18*(cos(27*pi/180)^0.28));
    
    Y_pt=Wing_Structural_Data.Wingbox_property.NodesY;
    
    Cd0=Cf.*FF.*Q.*AreaR;
    
    drag_span=2*Cd0.*(2*Bwidth)*Dyn_Pressure;

    Wing_drag_force=trapz(Y_pt,drag_span);
    
    wing_cd0=Wing_drag_force/(Dyn_Pressure*Wing_Area);
    
    % fuselage
    
    Fuselage_length=45;
    Fuselage_diameter=4;
    FF_fuselage=1+60/(Fuselage_length/Fuselage_diameter)^3+Fuselage_length/(400*Fuselage_diameter);
%     AreaR_fuselage=(2*pi*Fuselage_diameter*0.5*Fuselage_length)/(Fuselage_diameter*Fuselage_length);
    AreaR_fuselage=(2*pi*Fuselage_diameter*0.5*Fuselage_length)/Wing_Area;
    
    % Cut off Re number
    Xc_fuse=mu*Re_c/(rho*Fuselage_length*U);
    Re_fusecut=38.21*(Fuselage_length/Roughness)^1.053;
    
    % Reynold's number
    Re_fuselage_=rho*U*Fuselage_length/mu;
    Re_fuselage=min(Re_fuselage_,Re_fusecut);
    
    % drag coefficient calculation
    C_lam_fuselage=1.328/sqrt(Re_fuselage);
    C_turb_fuselage=0.455/(log10(Re_fuselage)^2.58*(1+0.144*Ma^2)^0.65);
    Cf_fuselage=Xc*C_lam_fuselage+(1-Xc_fuse)*C_turb_fuselage;
    
    Cd0f=Cf_fuselage.*FF_fuselage.*Q.*AreaR_fuselage;
    
    % Engine
    nacelles_length=3.5;
    nacelles_diameter=2;
    FF_nacelles=1+0.35/(nacelles_length/nacelles_diameter);
    Cf_nacelles=C_turb;
%     AreaR_nacelles=(2*pi*nacelles_diameter*0.5*nacelles_length)/(nacelles_length*nacelles_diameter);
    AreaR_nacelles=(2*pi*nacelles_diameter*0.5*nacelles_length)/(Wing_Area);
    
    Q_en=1.3;
    Cd0n=Cf_nacelles.*FF_nacelles.*Q_en.*AreaR_nacelles;
    
    % wave drag 
    CD_wave=0.001;
    
    % total zero-lift drag
%     CD0_all=(Wing_Area*wing_cd0+0.8*565*Cd0f+14*Cd0n)/Wing_Area + CD_wave;
    
    CD0_All=wing_cd0 + Cd0f + Cd0n + CD_wave;
    
end