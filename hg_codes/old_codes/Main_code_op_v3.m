


%% -----------------------------------------------------------

x0=[-0.002,-0.001,-0.001,0.005,...
    -0.002,-0.001,-0.001,0.005,...
    -0.00002,-0.00001,-0.00001,0.00005];

lb=ones(1,12)*1e-7;

ub=ones(1,12)*1;

% x=0:0.01:1
% 
% a=-0.6;b=0.6;d=0.1;c=d-a-b;
% 
% f=@(x)a*x.^3+b*x.^2+c*x+d;
% 
% y=f(x)
% 
% plot(x,y)


% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% [xopt, fva,ef,output,lamda]=fmincon(@objective,x0,[],[],[],[],lb,[],@constraint,[]);

[xopt, fva,ef,output,lamda]=fmincon(@objective,x0,[],[],[],[],lb,ub,@constraint,[]);

% Define objective function 
function obj = objective(x)
 [wing_mass,total_mass]=Mass_calc_dev(x);
 obj=wing_mass;
 disp(obj)
 disp(x)
end

% Constraint
function [c, ceq]=constraint(x)

[RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_strg, Rsigma_crip, Rsigma_col]=BeamStress_calc_dev(x);

c(1)=max(RVon_skn/5e8)-1; 
c(2)=max(RVon_spr/5e8)-1;
c(3)=max(RVon_skn./Rsigmab_skn)-1;
c(4)=max(Rsigma_strg./Rsigma_crip)-1;
disp(c)

ceq=[];

end


















