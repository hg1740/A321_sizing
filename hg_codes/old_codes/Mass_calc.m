function [Wing_mass,Total_mass]=Mass_calc(x)

%     thickness1=x(1:25);
%     thickness2=x(26:50);
%     Astrg=x(51:75);
    
    thickness1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
        x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
        x(21),x(22),x(23),x(24),x(25)];
    
    thickness2=[x(26),x(27),x(28),x(29),x(30),x(31),x(32),x(33),x(34),x(35)...
        x(36),x(37),x(38),x(39),x(40),x(41),x(42),x(43),x(44),x(45)...
        x(46),x(47),x(48),x(49),x(50)];
    
    Astrg=[x(51),x(52),x(53),x(54),x(55),x(56),x(57),x(58),x(59),x(60)...
        x(61),x(62),x(63),x(64),x(65),x(66),x(67),x(68),x(69),x(70)...
        x(71),x(72),x(73),x(74),x(75)];
    
    %% secondary masses

%     MTOW=77000; %kg
%     b_ref=50;% m 
%     b_st=17.05;% m spanwise distance 
%     k_fle=1.3;
%     dese_ref=56; % N/m^2
%     W_ref=10e6; %kg
%     c_st=4; %m take as the mean chord length
%     f_spar=0.15;
%     kl_sla=0.8;
%     kd_sla=0.15;
%     S_ref=10; %m^2
%     S_st=57.8;
%     
%     V_ac=0.78*340;
%     q=0.5*0.7*(V_ac)^2;
%     q_ref=30^3;
%     
% 
%     % fixed leading edge
%     dense_fle=3.15*k_fle*dese_ref*(q/q_ref)^0.25*((MTOW*b_st)/(W_ref*b_ref))^0.145;
%     l_fle=b_st/cos(25*3.14/180);
%     d_fle=c_st*f_spar;
%     W_fle=dense_fle*l_fle*d_fle;
% 
%     % movable leading edge
%     l_sla=kl_sla*b_st/cos(25*3.14/180);
%     d_sla=c_st*kd_sla;
%     dense_mle=4.83*dese_ref*(l_sla*d_sla/S_ref)^0.183;
%     W_mle=dense_mle*l_sla*d_sla;
% 
%     % fixed trailing edge
%     dense_fte=2.6*dese_ref*(MTOW*b_st/(W_ref*b_ref))^0.0544;
%     l_fte=(b_st-c_st*sin(25*3.14/180))/cos(25*3.14/180);
%     d_fte=c_st*(1-0.65);
%     W_fte=dense_fte*l_fte*d_fte;
% 
%     % movable trailing edge 
%     S_ail=0.8*1.7;
%     S_fla=0.8*5;
%     dense_fla=1.7*1.2*1*dese_ref*(1+(MTOW/W_ref)^0.35);
%     dense_ail=3*dese_ref*1.54*(S_ail/S_ref)^0.044;
%     W_mte=dense_fla*S_fla+dense_ail*S_ail;

Wing_mass=8097/2;
W_fle=Wing_mass*0.05;
W_mle=Wing_mass*0.05;
W_fte=Wing_mass*0.08;
W_mte=Wing_mass*0.16;

W_secondary=W_fle+W_mle+W_fte+W_mte;


    %% Primary mass: Wing mass

    NumStrg=[12,11,11,10,9,8,8,7,7,7,6,6,6,6,6,5,5,5,5,5,4,4,4,4,4];
    Height=[0.616,0.596,0.577,0.56,0.54,0.52,0.50,0.48,0.46,0.44,0.42,0.398,0.38,0.36,0.34,0.318,0.298,0.278,0.258,0.238,0.218,0.198,0.179,0.1586,0.139];
    Width=[3,2.83, 2.66, 2.5, 2.32, 2.15, 1.98,1.8,1.76,1.7,1.66,1.6,1.56,1.51,1.46,1.42,1.37,1.32,1.274,1.22,1.171,1.122,1.07,1.024,0.975];
    etaRL=[0,0.75,1.5,2.25,3.0,3.75,4.5,5.25,6.0,6.76,7.52,8.27,9.0,9.78,10.54,11.29,12.0,12.8,13.56,14.32,15.0,15.83,16.58,17.34,18.09];

    cs_=(thickness1.*Height+thickness2.*Width)*2+Astrg.*NumStrg*2;
    a1=cs_(1:end-1);
    a2=cs_(2:end);
    mean_area=(a1+a2)/2;

    l2=etaRL(2:end);
    l1=etaRL(1:end-1);
    seg_length=l2-l1;

    Vol_wing=sum(mean_area.*seg_length);
    rho_al=2700;

    Wing_mass=Vol_wing*rho_al;
    
    %% total mass
    fuel_mass=4675;
    Total_mass=Wing_mass+W_secondary+fuel_mass;


end