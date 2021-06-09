
%% optimisation 

% %Initial thicknesses1
% t11=0.001;t21=0.001;t31=0.001;t41=0.001;t51=0.001;
% t61=0.001;t71=0.001;t81=0.001;t91=0.001;t101=0.001;
% t111=0.001;t121=0.001;t131=0.001;t141=0.001;t151=0.001;
% t161=0.001;t171=0.001;t181=0.001;t191=0.001;t201=0.001;
% t211=0.001;t221=0.001;t231=0.001;t241=0.001;t251=0.001;
% 
% %Initial thicknesses2
% t12=0.001;t22=0.001;t32=0.001;t42=0.001;t52=0.001;
% t62=0.001;t72=0.001;t82=0.001;t92=0.001;t102=0.001;
% t112=0.001;t122=0.001;t132=0.001;t142=0.001;t152=0.001;
% t162=0.001;t172=0.001;t182=0.001;t192=0.001;t202=0.001;
% t212=0.001;t222=0.001;t232=0.001;t242=0.001;t252=0.001;
% 
% 
% %Initial Asg
% 
% Asg1=0.001,Asg2,Asg3,Asg4,Asg5,Asg6,Asg7,Asg8,Asg9,Asg10,Asg11,Asg12,Asg13,Asg14,Asg15,Asg16,Asg17,Asg18,Asg19,Asg20,Asg21,Asg22,Asg23,Asg24,Asg25
% 
% x0=[t11,t21,t31,t41,t51,t61,t71,t81,t91,t101,t111,t121,t131,t141,t151,t161,t171,t181,t191,t201,t211,t221,t231,t241,t251...
%     t12,t22,t32,t42,t52,t62,t72,t82,t92,t102,t112,t122,t132,t142,t152,t162,t172,t182,t192,t202,t212,t222,t232,t242,t252...
%     Asg1,Asg2,Asg3,Asg4,Asg5,Asg6,Asg7,Asg8,Asg9,Asg10,Asg11,Asg12,Asg13,Asg14,Asg15,Asg16,Asg17,Asg18,Asg19,Asg20,Asg21,Asg22,Asg23,Asg24,Asg25];

%% code testing

% x=[ones(1,25)*0.003,ones(1,25)*0.003,ones(1,25)*0.00001];
% [wing_mass,total_mass]=Mass_calc(x0);
% [RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_strg, Rsigma_crip, Rsigma_col]=BeamStress_calc(x0);
% F061=NastranMethods1.extractNastranResults('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\ignore\A320_half_model_SOL144.f06','ReadF06',true,'ReadHDF5',false);
% 
% % extract forces
% M_P1=[F061.f06data.Bendingmoment.UMPLN1(69:92),F061.f06data.Bendingmoment.LMPLN1(92)]; % in-plane moment
% M_P2=[F061.f06data.Bendingmoment.UMPLN2(69:92),F061.f06data.Bendingmoment.LMPLN2(92)]; % out of plane moment
% T=[F061.f06data.Bendingmoment.UTORQUE1(69:92),F061.f06data.Bendingmoment.LTORQUE1(92)];% torque
% 
% S_P1=[F061.f06data.Bendingmoment.USPLN1(69:92),F061.f06data.Bendingmoment.LSPLN1(92)]; % in plane shear
% S_P2=[F061.f06data.Bendingmoment.USPLN2(69:92),F061.f06data.Bendingmoment.LSPLN2(92)]; % out of plane shear
% 
% plot(M_P2)
% plot(T)
% plot(S_P2)
% plot(RVon_skn)
% plot(RVon_skn)


%% -----------------------------------------------------------

%  x0=[ones(1,25)*0.005,ones(1,25)*0.005,ones(1,25)*0.0005];


lb_t1=ones(1,25)*0.0001;
lb_t2=ones(1,25)*0.0001;
lb_A=ones(1,25)*1e-7;

lb=[lb_t1,lb_t2,lb_A];

ub_t1=ones(1,25)*0.1;
ub_t2=ones(1,25)*0.1;
ub_A=ones(1,25)*0.1;

ub=[ub_t1,ub_t2,ub_A];

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[xopt, fva,ef,output,lamda]=fmincon(@objective,x0,[],[],[],[],lb,ub,@constraint,options);



% Define objective function 

function obj = objective(x)
 [wing_mass,total_mass]=Mass_calc(x);
 obj=wing_mass;
 
 disp(obj)
end

% Constraint
function [c, ceq]=constraint(x)

[RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_strg, Rsigma_crip, Rsigma_col]=BeamStress_calc(x);

%tensile
cond1=max(RVon_skn/5e8); 
cond2=max(RVon_spr/5e8);

%compressive
cond3=max(RVon_skn./Rsigmab_skn);
cond4=max(Rsigma_strg./Rsigma_crip);

c= max([cond1,cond2,cond3,cond4]) - 1;
ceq=[];

end





 

















