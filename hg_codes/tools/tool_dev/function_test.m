Param=eval('A321');

% update FWT span 
Param.FWT.Fold_eta=0.8;

run_folder='C:\Git\A321_sizing\hg_codes\results\test_temp';


%% HL test 
[Static_Loads,Delta,Upper_Bound]=Stress_Analysis_HL(Param, run_folder);


figure 
plot(Static_Loads.Y,Static_Loads.Moment_P2,'s-')
hold on 
plot(Static_Loads.Y,Upper_Bound.Moment_P2,'v-')


figure 
plot(Static_Loads.Y(1:end-1),Delta.Wing_Delta.Wing_threshold.Max_Moment,'s-')
hold on 
plot(Static_Loads.Y(1:end-1),Delta.Wing_Delta.Wing_threshold.Min_Moment,'s-')

hold on 

plot(Static_Loads.Y(1:end-1),Delta.Wing_Delta.Wing_hinge_failure.Max_Moment,'v-')
hold on 
plot(Static_Loads.Y(1:end-1),Delta.Wing_Delta.Wing_hinge_failure.Min_Moment,'v-')


figure 
Moment1=Static_Loads.Moment_P2 + [Delta.Wing_Delta.Wing_threshold.Max_Moment;0];

Moment2=(Static_Loads.Moment_P2 + [Delta.Wing_Delta.Wing_hinge_failure.Max_Moment;0])*(2/3);

plot(Static_Loads.Y,Moment1,'s-')
hold on 
plot(Static_Loads.Y,Moment2,'v-')



%% ZF test

[ZF_Static_Loads,ZF_Upper_Bound]=Stress_Analysis_ZF(Param, run_folder);

%% sizing 


tic
[Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis_Fast_v1(Param, run_folder)
toc  


figure 

plot(Y_all,Load_distribution.Moment_P2.LC1,'s-')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC2,'s-')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC3,'s-')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC6,'s-')
hold on 
plot(Y_all,Load_distribution.Moment_P2.LC9,'s-')
hold on 
plot(Y_all,Load_distribution.Moment_P2.HL,'v-')
hold on 
plot(Y_all,Load_distribution.Moment_P2.ZF,'o-')

legend('2.5g cruise','2.5g sea level','-g dive','g+1mc','g-1mc','HL','ZF')


figure 

plot(Y_all,Load_distribution.Shear_P2.LC1,'s-')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC2,'s-')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC3,'s-')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC6,'s-')
hold on 
plot(Y_all,Load_distribution.Shear_P2.LC9,'s-')
hold on 
plot(Y_all,Load_distribution.Shear_P2.HL,'v-')
hold on 
plot(Y_all,Load_distribution.Shear_P2.ZF,'o-')

legend('2.5g cruise','2.5g sea level','-g dive','g+1mc','g-1mc','HL','ZF')



figure 

plot(Static_Loads.Y,1.5*Static_Loads.Moment_P2,'s-')

hold on 

plot(Static_Loads.Y,Upper_Bound.Moment_P2,'v-')

hold on 

plot(Static_Loads.Y,Static_Loads.Moment_P2+[Delta.Wing.Max_Moment;0],'o-')



figure 

plot(Static_Loads.Y,1.5*Static_Loads.Shear_P2,'s-')

hold on 

plot(Static_Loads.Y,Upper_Bound.Shear_P2,'v-')

hold on 

plot(Static_Loads.Y,Static_Loads.Shear_P2+[Delta.Wing.Max_Shear;0],'o-')


figure 

plot(Static_Loads.Y(1:end-1),Delta.Wing.Max_Moment,'v-')

figure 

plot(Static_Loads.Y(1:end-1),Delta.Wing.Max_Shear,'v-')



Loads_gust1=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test1','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off','Gust_Eta',1);

Loads_gust2=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test1','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off','Gust_Eta',0.25);

Loads_gust3=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test1','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off','Gust_Eta',1);

figure 

% plot(Loads_gust1.Y(1:end-1),Loads_gust1.Wing_Delta.Max_Moment,'s-')
% hold on 

plot(Loads_gust2.Y(1:end-1),Loads_gust2.Wing_Delta.Max_Moment,'v-')
hold on 
plot(Loads_gust3.Y(1:end-1),Loads_gust3.Wing_Delta.Max_Moment,'s-')


% if ~isfolder(fullfile(pwd,'bin'))
%    mkdir(fullfile(pwd,'bin'))
% else
%    delete(fullfile('bin','*'))
% end
% run_folder = fullfile(pwd,'bin');
% 
% [CDi,CD0,CL,k,Distribution,Load_distribution, Displacements_Res]=Static_Trim_v1(Param, run_folder, 'Load_Factor',2.5,'File_Name','pull_up_test','Hinge_Lock','off');
% 
% Loads=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');
% Loads_on=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','on');
% 
% figure 
% plot(Load_distribution.Y,Load_distribution.Moment_P2,'s-')
% 
% 
% figure 
% plot(Loads.Y(1:end-1),Loads.Wing_Delta.Max_Moment,'s-')
% hold on 
% plot(Loads_on.Y(1:end-1),Loads_on.Wing_Delta.Max_Moment,'v-')