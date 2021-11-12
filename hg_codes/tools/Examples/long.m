
Param=eval('A321');

d = 10; % Number of random variables 
% n = 50; % Number of emulations in surrogate model (arbitrarily large number)
% % choose 50 considering the time spent on Nastran as it takes above half an
% % hour to calculate 50 results with Nastran

 

N =1; % Number of training samples

 

%--------------------------------------------------------------------------
%% lower and upper bounds of variables
lower = [1/1000   1/1000   1/1000  1/1000   1/1000  0.1/1000   0.1/1000   0.1/1000  0.1/1000   0.1/1000 ]; %lower bounds 
upper = [30/1000  30/1000 30/1000  30/1000  30/1000 10/1000  10/1000 10/1000  10/1000  10/1000]; %upper bounds
param = [lower;upper]; % input required for sampling

 

xdd=[0.0473207320000000,0.359121570000000,0.909435710000000,0.919022650000000,0.414010470000000,0.376680380000000,0.463936930000000,0.261764570000000,0.543042370000000,0.769543070000000];

 

xxu=xdd.*repmat(param(2,:) - param(1,:),N,1) + repmat(param(1,:),N,1);

 

conshape=ones(1,36);
linearshape=1-(1/36:1/36:1);
mlinear=-1*linearshape+1;%another linear shape
% quadraticshape=((0.1:0.1:1)-1).*(2*(0.1:0.1:1)-1);
quadraticshape=-1*(1/36:1/36:1).*(1/36:1/36:1)+1;%quadratic shape
mquadratic=-1*quadraticshape+1;%another quadratic shape

 

iiii=1;
    %spar thickness dis
    shape{iiii}=xxu(iiii,1)*conshape+xxu(iiii,2)*linearshape+xxu(iiii,3)*quadraticshape+xxu(iiii,4)*mlinear+xxu(iiii,5)*mquadratic;
shapetotal=shape{iiii};
%skin thickness dis
shape2{iiii}=xxu(iiii,6)*conshape+xxu(iiii,7)*linearshape+xxu(iiii,8)*quadraticshape+xxu(iiii,9)*mlinear+xxu(iiii,10)*mquadratic;

 

shapetotal2=shape2{iiii};  

 

% shape3{iiii}=xxu(iiii,11)*conshape+xxu(iiii,12)*linearshape+xxu(iiii,13)*quadraticshape+xxu(iiii,14)*mlinear+xxu(iiii,15)*mquadratic;
% 
% shapetotal3=shape3{iiii}; 

 

    xxuse=[shapetotal(1:25),shapetotal2(1:25),5.3*10^(-5)*ones(1,25)];
Param.Wing.Thickness=xxuse;

 

fwt_x=[shapetotal(26:36),shapetotal2(26:36),5.3*10^(-5)*ones(1,11)];
Param.FWT.Thickness=fwt_x;


run_folder='C:\Git\bin';

[Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis_Fast_v1(Param, run_folder)
    
