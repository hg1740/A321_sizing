
figure 

plot(Y_all,Load_distribution.Moment_P2.LC1,'b-s')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC2,'r-s')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC3,'m-s')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC6,'b-v')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC6,'r-v')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC6,'m-v')
hold on 
plot(Y_all,Load_distribution.Moment_P2.Cruise,'k-s')



figure 

plot(Y_all,Load_distribution.Shear_P2.LC1,'b-s')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC2,'r-s')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC3,'m-s')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC6,'b-v')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC6,'r-v')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC6,'m-v')
hold on 
plot(Y_all,Load_distribution.Shear_P2.Cruise,'k-s')





Load_AR10_Eta0=load('Res_AR10_Eta_100_Loads');
Load_AR10_Eta10=load('Res_AR10_Eta_90_Loads');
Load_AR10_Eta20=load('Res_AR10_Eta_80_Loads');
Load_AR10_Eta30=load('Res_AR10_Eta_70_Loads');
Load_AR10_Eta40=load('Res_AR10_Eta_60_Loads');

figure 
plot(Load_AR10_Eta0.Loads.Y, Load_AR10_Eta0.Loads.Moment.pullup_cruise,'bs-')
hold on
plot(Load_AR10_Eta10.Loads.Y, Load_AR10_Eta10.Loads.Moment.pullup_cruise,'ks-')
hold on
plot(Load_AR10_Eta20.Loads.Y, Load_AR10_Eta20.Loads.Moment.pullup_cruise,'rs-')
hold on
plot(Load_AR10_Eta30.Loads.Y, Load_AR10_Eta30.Loads.Moment.pullup_cruise,'gs-')
hold on
plot(Load_AR10_Eta40.Loads.Y, Load_AR10_Eta40.Loads.Moment.pullup_cruise,'ms-')


figure 
plot(Load_AR10_Eta0.Loads.Y, Load_AR10_Eta0.Loads.Shear.pullup_cruise,'bs-')
hold on
plot(Load_AR10_Eta10.Loads.Y, Load_AR10_Eta10.Loads.Shear.pullup_cruise,'ks-')
hold on
plot(Load_AR10_Eta20.Loads.Y, Load_AR10_Eta20.Loads.Shear.pullup_cruise,'rs-')
hold on
plot(Load_AR10_Eta30.Loads.Y, Load_AR10_Eta30.Loads.Shear.pullup_cruise,'gs-')
hold on
plot(Load_AR10_Eta40.Loads.Y, Load_AR10_Eta40.Loads.Shear.pullup_cruise,'ms-')


%% tail wing 

% nodes ID

%% Loads distributions from different load cases

% Model 1
Param_AR10_Eta0=load('Res_AR10_Eta_100_Model');

Param_AR10_Eta0=Param_AR10_Eta0.Param;

% Mode 2 
Param_AR10_Eta10=load('Res_AR10_Eta_90_Model');

Param_AR10_Eta10=Param_AR10_Eta10.Param;

% Mode 3
Param_AR10_Eta20=load('Res_AR10_Eta_80_Model');

Param_AR10_Eta20=Param_AR10_Eta20.Param;

% Mode 4
Param_AR10_Eta30=load('Res_AR10_Eta_70_Model');

Param_AR10_Eta30=Param_AR10_Eta30.Param;

% Mode 5
Param_AR10_Eta40=load('Res_AR10_Eta_60_Model');

Param_AR10_Eta40=Param_AR10_Eta40.Param;


run_folder1='C:\Git\A321_sizing\hg_codes\results\AR10_test\eta0';
run_folder2='C:\Git\A321_sizing\hg_codes\results\AR10_test\eta10';
run_folder3='C:\Git\A321_sizing\hg_codes\results\AR10_test\eta20';
run_folder4='C:\Git\A321_sizing\hg_codes\results\AR10_test\eta30';
run_folder5='C:\Git\A321_sizing\hg_codes\results\AR10_test\eta40';

[Aircraft1, Wingbox_right1, FEM_full1, Y_all1, Box_dimensions1, Box_CrossSec1, Load_distribution1, Wing_Delta1, Root_Delta1, Internal_Stresses1]=Stress_Analysis_Fast_v1(Param_AR10_Eta0, run_folder1);
            
tic
[Aircraft2, Wingbox_right2, FEM_full2, Y_all2, Box_dimensions2, Box_CrossSec2, Load_distribution2, Wing_Delta2, Root_Delta2, Internal_Stresses2]=Stress_Analysis_Fast_v1(Param_AR10_Eta10, run_folder2);
toc    

[Aircraft3, Wingbox_right3, FEM_full3, Y_all3, Box_dimensions3, Box_CrossSec3, Load_distribution3, Wing_Delta3, Root_Delta3, Internal_Stresses3]=Stress_Analysis_Fast_v1(Param_AR10_Eta20, run_folder3);
    

[Aircraft4, Wingbox_right4, FEM_full4, Y_all4, Box_dimensions4, Box_CrossSec4, Load_distribution4, Wing_Delta4, Root_Delta4, Internal_Stresses4]=Stress_Analysis_Fast_v1(Param_AR10_Eta30, run_folder4);
    

[Aircraft5, Wingbox_right5, FEM_full5, Y_all5, Box_dimensions5, Box_CrossSec5, Load_distribution5, Wing_Delta5, Root_Delta5, Internal_Stresses5]=Stress_Analysis_Fast_v1(Param_AR10_Eta40, run_folder5);

figure % Moment eta=0

plot(Y_all1,Load_distribution1.Moment_P2.LC1,'rs-')
hold on 
plot(Y_all1,Load_distribution1.Moment_P2.LC2,'bs-')
hold on 
plot(Y_all1,Load_distribution1.Moment_P2.LC3,'ks-')
hold on 
plot(Y_all1,Load_distribution1.Moment_P2.LC6,'ms-')
hold on 
plot(Y_all1,Load_distribution1.Moment_P2.LC9,'gs-')

figure % Shear eta=0

plot(Y_all1,Load_distribution1.Shear_P2.LC1,'rs-')
hold on 
plot(Y_all1,Load_distribution1.Shear_P2.LC2,'bs-')
hold on 
plot(Y_all1,Load_distribution1.Shear_P2.LC3,'ks-')
hold on 
plot(Y_all1,Load_distribution1.Shear_P2.LC6,'ms-')
hold on 
plot(Y_all1,Load_distribution1.Shear_P2.LC9,'gs-')




figure % Moment eta=10

plot(Y_all2,Load_distribution2.Moment_P2.LC1,'rs-')
hold on 
plot(Y_all2,Load_distribution2.Moment_P2.LC2,'bs-')
hold on 
plot(Y_all2,Load_distribution2.Moment_P2.LC3,'ks-')
hold on 
plot(Y_all2,Load_distribution2.Moment_P2.LC6,'ms-')
hold on 
plot(Y_all2,Load_distribution2.Moment_P2.LC9,'gs-')
hold on 
plot(Y_all2,Load_distribution2.Moment_P2.Cruise,'k--')

figure % Shear eta=10

plot(Y_all2,Load_distribution2.Shear_P2.LC1,'rs-')
hold on 
plot(Y_all2,Load_distribution2.Shear_P2.LC2,'bs-')
hold on 
plot(Y_all2,Load_distribution2.Shear_P2.LC3,'ks-')
hold on 
plot(Y_all2,Load_distribution2.Shear_P2.LC6,'ms-')
hold on 
plot(Y_all2,Load_distribution2.Shear_P2.LC9,'gs-')
hold on 
plot(Y_all2,Load_distribution2.Shear_P2.Cruise,'k--')


figure % Moment eta=20

plot(Y_all3,Load_distribution3.Moment_P2.LC1,'rs-')
hold on 
plot(Y_all3,Load_distribution3.Moment_P2.LC2,'bs-')
hold on 
plot(Y_all3,Load_distribution3.Moment_P2.LC3,'ks-')
hold on 
plot(Y_all3,Load_distribution3.Moment_P2.LC6,'ms-')
hold on 
plot(Y_all3,Load_distribution3.Moment_P2.LC9,'gs-')
hold on 
plot(Y_all3,Load_distribution3.Moment_P2.Cruise,'k--')

figure % Shear eta=20

plot(Y_all3,Load_distribution3.Shear_P2.LC1,'rs-')
hold on 
plot(Y_all3,Load_distribution3.Shear_P2.LC2,'bs-')
hold on 
plot(Y_all3,Load_distribution3.Shear_P2.LC3,'ks-')
hold on 
plot(Y_all3,Load_distribution3.Shear_P2.LC6,'ms-')
hold on 
plot(Y_all3,Load_distribution3.Shear_P2.LC9,'gs-')
hold on 
plot(Y_all3,Load_distribution3.Shear_P2.Cruise,'k--')



figure % Moment eta=30

plot(Y_all4,Load_distribution4.Moment_P2.LC1,'rs-')
hold on 
plot(Y_all4,Load_distribution4.Moment_P2.LC2,'bs-')
hold on 
plot(Y_all4,Load_distribution4.Moment_P2.LC3,'ks-')
hold on 
plot(Y_all4,Load_distribution4.Moment_P2.LC6,'ms-')
hold on 
plot(Y_all4,Load_distribution4.Moment_P2.LC9,'gs-')
hold on 
plot(Y_all4,Load_distribution4.Moment_P2.Cruise,'k--')

figure % Shear eta=30

plot(Y_all4,Load_distribution4.Shear_P2.LC1,'rs-')
hold on 
plot(Y_all4,Load_distribution4.Shear_P2.LC2,'bs-')
hold on 
plot(Y_all4,Load_distribution4.Shear_P2.LC3,'ks-')
hold on 
plot(Y_all4,Load_distribution4.Shear_P2.LC6,'ms-')
hold on 
plot(Y_all4,Load_distribution4.Shear_P2.LC9,'gs-')
hold on 
plot(Y_all4,Load_distribution4.Shear_P2.Cruise,'k--')


figure % Moment eta=40

plot(Y_all5,Load_distribution5.Moment_P2.LC1,'rs-')
hold on 
plot(Y_all5,Load_distribution5.Moment_P2.LC2,'bs-')
hold on 
plot(Y_all5,Load_distribution5.Moment_P2.LC3,'ks-')
hold on 
plot(Y_all5,Load_distribution5.Moment_P2.LC6,'ms-')
hold on 
plot(Y_all5,Load_distribution5.Moment_P2.LC9,'gs-')
hold on 
plot(Y_all5,0.67*Load_distribution5.Moment_P2.Cruise,'k--')

figure % Shear eta=40

plot(Y_all5,Load_distribution5.Shear_P2.LC1,'rs-')
hold on 
plot(Y_all5,Load_distribution5.Shear_P2.LC2,'bs-')
hold on 
plot(Y_all5,Load_distribution5.Shear_P2.LC3,'ks-')
hold on 
plot(Y_all5,Load_distribution5.Shear_P2.LC6,'ms-')
hold on 
plot(Y_all5,Load_distribution5.Shear_P2.LC9,'gs-')
hold on 
plot(Y_all5,Load_distribution5.Shear_P2.Cruise,'k--')

