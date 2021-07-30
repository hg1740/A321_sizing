
% Create Aircraft Model

run_folder='C:\Git\A321_sizing\hg_codes\results\hinge_angle_test';

Wing_Model=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR16\Res_AR16_Eta_70_Model.mat');


% Hinge angle 40 model 

Param_H40=Wing_Model.Param;

% update hinge angle 

Param_H40.FWT.Flare_angle=40/2;


% Hinge angle 35 model 

Param_H35=Wing_Model.Param;

% update hinge angle 

Param_H35.FWT.Flare_angle=35/2;



% Hinge angle 30 model 

Param_H30=Wing_Model.Param;

% update hinge angle 

Param_H30.FWT.Flare_angle=30/2;


% Hinge angle 25 model 

Param_H25=Wing_Model.Param;

% update hinge angle 

Param_H25.FWT.Flare_angle=25/2;


% Hinge angle 20 model 

Param_H20=Wing_Model.Param;

% update hinge angle 

Param_H20.FWT.Flare_angle=20/2;



% Hinge angle 15 model 

Param_H15=Wing_Model.Param;

% update hinge angle 

Param_H15.FWT.Flare_angle=15/2;



% Hinge angle 10 model 

Param_H10=Wing_Model.Param;

% update hinge angle 

Param_H10.FWT.Flare_angle=10/2;


% Hinge angle 5 model 

Param_H5=Wing_Model.Param;

% update hinge angle 

Param_H5.FWT.Flare_angle=5/2;



% Hinge angle 0 model 

Param_H0=Wing_Model.Param;

% update hinge angle 

Param_H0.FWT.Flare_angle=0/2;


%% Run static trim 

[FEM_full_H40,~,~,~,~,Distribution_H40,Load_distribution_H40, Displacements_Res_H40]=Static_Trim_v1(Param_H40, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_H35,~,~,~,~,Distribution_H35,Load_distribution_H35, Displacements_Res_H35]=Static_Trim_v1(Param_H35, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_H30,~,~,~,~,Distribution_H30,Load_distribution_H30, Displacements_Res_30]=Static_Trim_v1(Param_H30, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_H25,~,~,~,~,Distribution_H25,Load_distribution_H25, Displacements_Res_H25]=Static_Trim_v1(Param_H25, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_H20,~,~,~,~,Distribution_H20,Load_distribution_H20, Displacements_Res_H20]=Static_Trim_v1(Param_H20, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_H15,~,~,~,~,Distribution_H15,Load_distribution_H15, Displacements_Res_H15]=Static_Trim_v1(Param_H15, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_H10,~,~,~,~,Distribution_H10,Load_distribution_H10, Displacements_Res_H10]=Static_Trim_v1(Param_H10, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_H5,~,~,~,~,Distribution_H5,Load_distribution_H5, Displacements_Res_H5]=Static_Trim_v1(Param_H5, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_H0,~,~,~,~,Distribution_H0,Load_distribution_H0, Displacements_Res_H0]=Static_Trim_v1(Param_H0, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

[FEM_full_HL,~,~,~,~,Distribution_HL,Load_distribution_HL, Displacements_Res_HL]=Static_Trim_v1(Param_H0, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','on');

%% Run Gust 

Loads_H40=Gust_Analysis_v1(Param_H40,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');

Loads_H35=Gust_Analysis_v1(Param_H35,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');

Loads_H30=Gust_Analysis_v1(Param_H30,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');

Loads_H25=Gust_Analysis_v1(Param_H25,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');

Loads_H20=Gust_Analysis_v1(Param_H20,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');

Loads_H15=Gust_Analysis_v1(Param_H15,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');



%% Static Results

% Moment

figure

plot(Load_distribution_H40.Y,Load_distribution_H40.Moment_P2,'rs-')
hold on 
plot(Load_distribution_H35.Y,Load_distribution_H35.Moment_P2,'ks-')
hold on 
plot(Load_distribution_H30.Y,Load_distribution_H30.Moment_P2,'ms-')
hold on 
plot(Load_distribution_H25.Y,Load_distribution_H25.Moment_P2,'o-')
hold on 
plot(Load_distribution_H20.Y,Load_distribution_H20.Moment_P2,'gs-')
hold on 
plot(Load_distribution_H15.Y,Load_distribution_H15.Moment_P2,'bs-')
hold on 
plot(Load_distribution_H10.Y,Load_distribution_H10.Moment_P2,'r--','LineWidth',1.2)
% hold on 
% plot(Load_distribution_H5.Y,Load_distribution_H5.Moment_P2,'b--','LineWidth',1.2)
hold on 
plot(Load_distribution_H0.Y,Load_distribution_H0.Moment_P2,'k--','LineWidth',1.2)
hold on 
plot(Load_distribution_HL.Y,Load_distribution_HL.Moment_P2,'k--','LineWidth',1.2)

% Shear

figure

plot(Load_distribution_H40.Y,Load_distribution_H40.Shear_P2,'rs-')
hold on 
plot(Load_distribution_H35.Y,Load_distribution_H35.Shear_P2,'ks-')
hold on 
plot(Load_distribution_H30.Y,Load_distribution_H30.Shear_P2,'ms-')
hold on 
plot(Load_distribution_H25.Y,Load_distribution_H25.Shear_P2,'o-')
hold on 
plot(Load_distribution_H20.Y,Load_distribution_H20.Shear_P2,'gs-')
hold on 
plot(Load_distribution_H15.Y,Load_distribution_H15.Shear_P2,'bs-')
hold on 
plot(Load_distribution_H10.Y,Load_distribution_H10.Shear_P2,'rs-')
% hold on 
% plot(Load_distribution_H5.Y,Load_distribution_H5.Shear_P2,'bs-')
hold on 
plot(Load_distribution_H0.Y,Load_distribution_H0.Shear_P2,'k--','LineWidth',1.2)
hold on 
plot(Load_distribution_HL.Y,Load_distribution_HL.Shear_P2,'k--','LineWidth',1.2)





% Lift distribution 

figure 

plot(Distribution_H0.Y,Distribution_H0.Lift_var,'b-')
hold on 
plot(Distribution_H5.Y,Distribution_H10.Lift_var,'r-')
hold on 
plot(Distribution_H15.Y,Distribution_H15.Lift_var,'k-')
hold on 
plot(Distribution_H25.Y,Distribution_H25.Lift_var,'g-')
hold on 
plot(Distribution_H35.Y,Distribution_H35.Lift_var,'m-')
hold on 
plot(Distribution_H40.Y,Distribution_H40.Lift_var,'r--')

%% Gust Results


% Moment

figure

plot(Loads_H40.Y(1:end-1),Loads_H40.Wing_Delta.Max_Moment,'rs-')
hold on 
plot(Loads_H35.Y(1:end-1),Loads_H35.Wing_Delta.Max_Moment,'ks-')
hold on 
plot(Loads_H30.Y(1:end-1),Loads_H30.Wing_Delta.Max_Moment,'ms-')
hold on 
plot(Loads_H25.Y(1:end-1),Loads_H25.Wing_Delta.Max_Moment,'o-')
hold on 
plot(Loads_H20.Y(1:end-1),Loads_H20.Wing_Delta.Max_Moment,'gs-')
hold on 
plot(Loads_H15.Y(1:end-1),Loads_H15.Wing_Delta.Max_Moment,'bs-')

hold on 

plot(Loads_H40.Y(1:end-1),Loads_H40.Wing_Delta.Min_Moment,'rs-')
hold on 
plot(Loads_H35.Y(1:end-1),Loads_H35.Wing_Delta.Min_Moment,'ks-')
hold on 
plot(Loads_H30.Y(1:end-1),Loads_H30.Wing_Delta.Min_Moment,'ms-')
hold on 
plot(Loads_H25.Y(1:end-1),Loads_H25.Wing_Delta.Min_Moment,'o-')
hold on 
plot(Loads_H20.Y(1:end-1),Loads_H20.Wing_Delta.Min_Moment,'gs-')
hold on 
plot(Loads_H15.Y(1:end-1),Loads_H15.Wing_Delta.Min_Moment,'bs-')





% Shear

figure

plot(Loads_H40.Y(1:end-1),Loads_H40.Wing_Delta.Max_Shear,'rs-')
hold on 
plot(Loads_H35.Y(1:end-1),Loads_H35.Wing_Delta.Max_Shear,'ks-')
hold on 
plot(Loads_H30.Y(1:end-1),Loads_H30.Wing_Delta.Max_Shear,'ms-')
hold on 
plot(Loads_H25.Y(1:end-1),Loads_H25.Wing_Delta.Max_Shear,'o-')
hold on 
plot(Loads_H20.Y(1:end-1),Loads_H20.Wing_Delta.Max_Shear,'gs-')
hold on 
plot(Loads_H15.Y(1:end-1),Loads_H15.Wing_Delta.Max_Shear,'bs-')

hold on 


plot(Loads_H40.Y(1:end-1),Loads_H40.Wing_Delta.Min_Shear,'rs-')
hold on 
plot(Loads_H35.Y(1:end-1),Loads_H35.Wing_Delta.Min_Shear,'ks-')
hold on 
plot(Loads_H30.Y(1:end-1),Loads_H30.Wing_Delta.Min_Shear,'ms-')
hold on 
plot(Loads_H25.Y(1:end-1),Loads_H25.Wing_Delta.Min_Shear,'o-')
hold on 
plot(Loads_H20.Y(1:end-1),Loads_H20.Wing_Delta.Min_Shear,'gs-')
hold on 
plot(Loads_H15.Y(1:end-1),Loads_H15.Wing_Delta.Min_Shear,'bs-')
