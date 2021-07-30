function [Wingbox_mass,Secondary_mass,Total_mass]=Wing_Mass_Calc_v1(Wingbox_property, Param)


    %% Box-mass
    a1 = Wingbox_property.A(1:end-1);
    a2 = Wingbox_property.A(2:end);
    mean_area=(a1+a2)/2;
    
    l2=Wingbox_property.NodesY(2:end);
    l1=Wingbox_property.NodesY(1:end-1);
    seg_length=l2-l1;
    
    Vol_wing=sum(mean_area.*seg_length);
    
    rho_al=2800;
    
    Wingbox_mass=Vol_wing*rho_al;
    

    %% secondary masses

    MTOW=Param.Masses.MTOW; %kg
    Surfacearea=Param.Wing.HalfArea;
    b_ref=50;% m 
    b_st=16;% m spanwise distance 
    k_fle=1.3;
    dese_ref=56; % N/m^2
    W_ref=1e6; %kg
    c_st=4.2; %m take as the mean chord length
    f_spar=0.15;
    kl_sla=0.8;
    kd_sla=0.15;
    S_ref=10; %m^2
    S_st=57.8;
    LE_Sweep=Param.Wing.LE_Sweep;
    k_sup=1.2;
    k_slot=1;
    k_bal=1.54;
    
    V_ac=0.78*340;
    q=0.5*0.4135*(V_ac)^2;
    q_ref=30^3;
    

    % fixed leading edge
    dense_fle=3.15*k_fle*dese_ref*(q/q_ref)^0.25*((MTOW*b_st)/(W_ref*b_ref))^0.145;
    l_fle=b_st/cos(LE_Sweep*pi/180);
    d_fle=c_st*f_spar;
    W_fle=dense_fle*l_fle*d_fle;

    % movable leading edge
    l_sla=kl_sla*b_st/cos(LE_Sweep*3.14/180);
    d_sla=c_st*kd_sla;
    dense_mle=4.83*dese_ref*(l_sla*d_sla/S_ref)^0.183;
    W_mle=dense_mle*l_sla*d_sla;

    % fixed trailing edge
    dense_fte=2.6*dese_ref*(MTOW*b_st/(W_ref*b_ref))^0.0544;
    l_fte=(b_st-c_st*sin(LE_Sweep*pi/180))/cos(LE_Sweep*pi/180);
    d_fte=c_st*(1-0.65);
    W_fte=dense_fte*l_fte*d_fte;

    % movable trailing edge 
    S_ail=Surfacearea*0.03;
    S_fla=Surfacearea*0.06;
    S_spo=Surfacearea*0.16; 
    dense_fla=1.7*k_sup*k_slot*dese_ref*(1+(MTOW/W_ref)^0.35);
    dense_ail=3*dese_ref*k_bal*(S_ail/S_ref)^0.044;
    dense_spo=2.2*dese_ref*(S_spo/S_ref)^0.032;
    
    W_mte=dense_fla*S_fla+dense_ail*S_ail + dense_spo*S_spo;

    W_secondary=W_fle+W_mle+W_fte+W_mte; 
    
    Secondary_mass=W_secondary/9.81; % weight to mass convert


   %% 
   
   Total_mass=Wingbox_mass+Secondary_mass;


end