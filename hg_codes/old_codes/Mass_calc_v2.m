function [Wing_mass,Total_mass]=Mass_calc_v2(x)

%     thickness1=x(1:25);
%     thickness2=x(26:50);
%     Astrg=x(51:75);
    
%     thickness1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
%         x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
%         x(21),x(22)];
% 
%     thickness2=[x(23),x(24),x(25),x(26),x(27),x(28),x(29),x(30),x(31),x(32)...
%         x(33),x(34),x(35),x(36),x(37),x(38),x(39),x(40),x(41),x(42)...
%         x(43),x(44)];
% 
%     Astrg=[x(45),x(46),x(47),x(48),x(49),x(50),x(51),x(52),x(53),x(54)...
%         x(55),x(56),x(57),x(58),x(59),x(60),x(61),x(62),x(63),x(64)...
%         x(65),x(66)];

    num=numel(x)/3;

    thickness1=x(1 : num);
    thickness2=x(num+1 : num*2);
    Astrg=x(num*2+1 : num*3);
    
%     Payload_max=25300;
%     Fuel_max=25000;
% 
%     MTOW=93500;
%     OWE=48500;
%     MWE=44057;
%     Wing_mass=8906;
% 
%     Engine_mass=7362/2;
%     Pylon=1239/2;
%     Fuselage_mass=8986;
    
    
    
    
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

    NumStrg=[12,11,11,10,9,9,8,7,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2];
    Height=[0.711000000000000,0.660875855073693,0.610751710147386,0.560627565221078,0.510503420294771,0.460379275368464,0.410255130442157,0.360130985515850,0.345731821009572,0.331332656503294,0.316933491997016,0.302534327490738,0.288135162984460,0.273735998478182,0.259336833971904,0.244937669465626,0.230538504959348,0.216139340453070,0.201740175946792,0.187341011440514,0.172941846934236,0.158542682427958,0.144143517921680,0.129744353415403,0.115345188909125];
    Width=[3,2.84277500415601,2.68555000831201,2.52832501246801,2.37110001662402,2.21387502078002,2.05665002493603,1.89942502909203,1.82673333439662,1.75404163970121,1.68134994500579,1.60865825031038,1.53596655561496,1.46327486091955,1.39058316622413,1.31789147152872,1.24519977683331,1.17250808213789,1.09981638744248,1.02712469274706,0.954432998051649,0.881741303356235,0.809049608660820,0.736357913965406,0.663666219269992];
    etaRL=17.1113*[0,0.0373518769000545,0.0747037538001090,0.112055630700163,0.149407507600218,0.186759384500272,0.224111261400327,0.261463138300381,0.304906483106241,0.348349827912101,0.391793172717961,0.435236517523821,0.478679862329681,0.522123207135541,0.565566551941401,0.609009896747261,0.652453241553121,0.695896586358981,0.739339931164840,0.782783275970700,0.826226620776560,0.869669965582420,0.913113310388280,0.956556655194140,1];

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
    fuel_mass=9500;
    Total_mass=Wing_mass+W_secondary+fuel_mass;


end