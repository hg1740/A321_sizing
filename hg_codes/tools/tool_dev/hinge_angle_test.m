run_folder='C:\Git\A321_sizing\hg_codes\results\hinge_angle_test';


Wing_Model=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR13\Res_AR13_Eta_80_Model.mat');

Param=Wing_Model.Param;

% update hinge angle 

Param.FWT.Flare_angle=-13;

[FEM_full,CDi,CD0,CL,k,Aerodynamic_distribution,Load_distribution,Displacements_Res]=Static_Trim_v1(Param, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

%% Result plot 

model = mni.import_matran(fullfile(run_folder,'Trim_analysis.dat'));
model.draw

%extract the data
f06 = mni.result.f06(fullfile(run_folder,'Trim_analysis.f06'));
res_disp =  f06.read_disp;
res_aeroP = f06.read_aeroP;
res_aeroF = f06.read_aeroF;

% apply deformation result
[~,i] = ismember(model.GRID.GID,res_disp.GP);
model.GRID.Deformation = [res_disp.dX(:,i);res_disp.dY(:,i);res_disp.dZ(:,i)];

% apply aero pressure
model.CAERO1.PanelPressure = res_aeroP.Cp;

%apply aero forces
f = [res_aeroF.aeroFx;res_aeroF.aeroFy;res_aeroF.aeroFz;...
    res_aeroF.aeroMx;res_aeroF.aeroMy;res_aeroF.aeroMz];
model.CAERO1.PanelForce = f';

% update the plot to apply deformations and aero pressures + forces
model.update('Scale',1)
% 