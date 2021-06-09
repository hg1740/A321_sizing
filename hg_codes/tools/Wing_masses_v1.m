function [Wingbox_mass, Secondary_mass, Total_mass]=Wing_masses_v1(x,Height,Width,etaRL,MTOW,Surfacearea,LE_Sweep,Semi_span)

    % Note: here the 'Wingbox_mass' is the structure mass of one wing box

    num=numel(x)/3;

    thickness1=x(1 : num);
    thickness2=x(num+1 : num*2);
    Astrg=x(num*2+1 : num*3);


    %% Primary mass: Wing mass
    strg_n=0.24;

    NumStrg=Width/strg_n;
%     Height=[0.711000000000000,0.660875855073693,0.610751710147386,0.560627565221078,0.510503420294771,0.460379275368464,0.410255130442157,0.360130985515850,0.345731821009572,0.331332656503294,0.316933491997016,0.302534327490738,0.288135162984460,0.273735998478182,0.259336833971904,0.244937669465626,0.230538504959348,0.216139340453070,0.201740175946792,0.187341011440514,0.172941846934236,0.158542682427958,0.144143517921680,0.129744353415403,0.115345188909125];
%     Width=[3,2.84277500415601,2.68555000831201,2.52832501246801,2.37110001662402,2.21387502078002,2.05665002493603,1.89942502909203,1.82673333439662,1.75404163970121,1.68134994500579,1.60865825031038,1.53596655561496,1.46327486091955,1.39058316622413,1.31789147152872,1.24519977683331,1.17250808213789,1.09981638744248,1.02712469274706,0.954432998051649,0.881741303356235,0.809049608660820,0.736357913965406,0.663666219269992];
%     etaRL=17.1113*[0,0.0373518769000545,0.0747037538001090,0.112055630700163,0.149407507600218,0.186759384500272,0.224111261400327,0.261463138300381,0.304906483106241,0.348349827912101,0.391793172717961,0.435236517523821,0.478679862329681,0.522123207135541,0.565566551941401,0.609009896747261,0.652453241553121,0.695896586358981,0.739339931164840,0.782783275970700,0.826226620776560,0.869669965582420,0.913113310388280,0.956556655194140,1];

    cs_=(thickness1.*Height+thickness2.*Width)*2+Astrg.*NumStrg*2;
    a1=cs_(1:end-1);
    a2=cs_(2:end);
    mean_area=(a1+a2)/2;

    l2=etaRL(2:end);
    l1=etaRL(1:end-1);
    seg_length=l2-l1;

    Vol_wing=sum(mean_area.*seg_length);
    rho_al=2700;

    Wingbox_mass=Vol_wing*rho_al;
    
    %% Secondary mass: LE flexible & fixed, TE felxible & fixed
    
    % Surface area here is the area of one wing without fuselage part at the
    % root
    
    % The secondary mass is calculated based the initial guess of the MTOW.
    
    b_ref=50;% m
    dese_ref=56; % N/m^2
    W_ref=1e6; %kg
    c_st=4.2; %m take as the mean chord length
    f_spar=0.15;
    S_ref=10; %m^2
    
    k_fle=1.3;
    k_sup=1.2;
    k_slot=1.5;
    k_bal=1.54;
    kl_sla=0.9;
    kd_sla=0.15;
    
    V_ac=0.78*340;
    q=0.5*0.4135*(V_ac)^2;
    q_ref=30^3;
    
    Span=Semi_span*2+4;
    
    % fixed leading edge
    dense_fle=3.15*k_fle*dese_ref*(q/q_ref)^0.25*((10*MTOW*Span)/(W_ref*b_ref))^0.145;
    l_fle=Semi_span/cos(LE_Sweep*pi/180);
    d_fle=c_st*f_spar;
    W_fle=dense_fle*l_fle*d_fle;
    
    % movable leading edge
    l_sla=kl_sla*Semi_span/cos(LE_Sweep*3.14/180);
    d_sla=c_st*kd_sla;
    dense_mle=4.83*dese_ref*(2*l_sla*d_sla/S_ref)^0.183;
    W_mle=dense_mle*l_sla*d_sla;
    
    % fixed trailing edge
    dense_fte=2.6*dese_ref*(10*MTOW*Span/(W_ref*b_ref))^0.0544;
    l_fte=0.5*(Span-c_st*sin(LE_Sweep*pi/180))/cos(LE_Sweep*pi/180);
    d_fte=c_st*(1-0.65);
    W_fte=dense_fte*l_fte*d_fte;
    
    % movable trailing edge
    S_ail=Surfacearea*0.05;
    S_fla=Surfacearea*0.072;
    S_spo=Surfacearea*0.16;
    dense_fla=1.7*k_sup*k_slot*dese_ref*(1+(10*MTOW/W_ref)^0.35);
    dense_ail=3*dese_ref*k_bal*(2*S_ail/S_ref)^0.044;
    dense_spo=2.2*dese_ref*(2*S_spo/S_ref)^0.032;
    
    W_mte=dense_fla*S_fla+dense_ail*S_ail + dense_spo*S_spo;
    
    % Miscellaneous
    t_paint=0.1;
    rho_paint=1;
    K_paint=2;
    
    W_tip=150*(10*MTOW/W_ref)^0.65;
    W_paint=t_paint*rho_paint*K_paint*Surfacearea;
    
    W_misc=W_tip+W_paint;
    
    % masses on the wing
    W_secondary=W_fle+W_mle+W_fte+W_mte+W_misc;
    Secondary_mass=W_secondary/10;
    
    %% Total mass
    
    Total_mass=Wingbox_mass+Secondary_mass;
    


end