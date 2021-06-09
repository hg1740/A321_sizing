function [Wing_mass,Total_mass]=Mass_calc_v1_FWT(x)

%     thickness1=x(1:25);
%     thickness2=x(26:50);
%     Astrg=x(51:75);
    
    thickness1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
        x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
        x(21)];

    thickness2=[x(22),x(23),x(24),x(25),x(26),x(27),x(28),x(29),x(30),x(31),x(32)...
        x(33),x(34),x(35),x(36),x(37),x(38),x(39),x(40),x(41),x(42)];
     

    Astrg=[x(43),x(44),x(45),x(46),x(47),x(48),x(49),x(50),x(51),x(52),x(53),x(54)...
        x(55),x(56),x(57),x(58),x(59),x(60),x(61),x(62),x(63)];
    
    
    
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

    NumStrg=[12,11,11,10,10,9,8,8,8,7,7,6,6,6,5,5,5,4,4,4,3];
    Height=[0.711000000000000,0.663706105896848,0.616412211793696,0.569118317690545,0.521824423587393,0.474530529484241,0.427236635381089,0.379942741277937,0.363036014302919,0.346129287327901,0.329222560352882,0.312315833377864,0.295409106402845,0.278502379427827,0.261595652452809,0.244688925477790,0.227782198502772,0.210875471527754,0.193968744552735,0.177062017577717,0.160155290602698];
    Width=[3,2.85770248740050,2.71540497480099,2.57310746220149,2.43080994960198,2.28851243700248,2.14621492440297,2.00391741180347,1.92065391017970,1.83739040855593,1.75412690693216,1.67086340530839,1.58759990368462,1.50433640206085,1.42107290043708,1.33780939881331,1.25454589718954,1.17128239556577,1.08801889394201,1.00475539231824,0.921491890694467];
    etaRL=[0,0.578571428571429,1.15714285714286,1.73571428571429,2.31428571428571,2.89285714285714,3.47142857142857,4.05000000000000,4.89230769230769,5.73461538461539,6.57692307692308,7.41923076923077,8.26153846153846,9.10384615384615,9.94615384615385,10.7884615384615,11.6307692307692,12.4730769230769,13.3153846153846,14.1576923076923,15];

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