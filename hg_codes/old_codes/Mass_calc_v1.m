function [Wing_mass,Total_mass]=Mass_calc_v1(x)

%     thickness1=x(1:25);
%     thickness2=x(26:50);
%     Astrg=x(51:75);
    
    thickness1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
        x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
        x(21),x(22)];

    thickness2=[x(23),x(24),x(25),x(26),x(27),x(28),x(29),x(30),x(31),x(32)...
        x(33),x(34),x(35),x(36),x(37),x(38),x(39),x(40),x(41),x(42)...
        x(43),x(44)];

    Astrg=[x(45),x(46),x(47),x(48),x(49),x(50),x(51),x(52),x(53),x(54)...
        x(55),x(56),x(57),x(58),x(59),x(60),x(61),x(62),x(63),x(64)...
        x(65),x(66)];
    
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

    NumStrg=[12,11,11,10,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,4,3,3];
    Height=[0.711000000000000,0.654695466486852,0.598390932973703,0.542086399460555,0.485781865947406,0.429477332434258,0.373172798921109,0.357676186633189,0.342179574345268,0.326682962057348,0.311186349769428,0.295689737481508,0.280193125193587,0.264696512905667,0.249199900617747,0.233703288329826,0.218206676041906,0.202710063753985,0.187213451466065,0.171716839178145,0.156220226890225,0.140723614602304];
    Width=[3,2.82803516079563,2.65607032159126,2.48410548238689,2.31214064318252,2.14017580397815,1.96821096477378,1.89097603915990,1.81374111354603,1.73650618793215,1.65927126231827,1.58203633670439,1.50480141109052,1.42756648547664,1.35033155986276,1.27309663424888,1.19586170863501,1.11862678302113,1.04139185740725,0.964156931793372,0.886922006179495,0.809687080565617];
    etaRL=[0,0.675000000000000,1.35000000000000,2.02500000000000,2.70000000000000,3.37500000000000,4.05000000000000,4.78000000000000,5.51000000000000,6.24000000000000,6.97000000000000,7.70000000000000,8.43000000000000,9.16000000000000,9.89000000000000,10.6200000000000,11.3500000000000,12.0800000000000,12.8100000000000,13.5400000000000,14.2700000000000,15];

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