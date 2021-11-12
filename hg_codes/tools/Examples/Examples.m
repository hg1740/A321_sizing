
%% Note %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Request for .OP4 out-put is located at:   ALENA\+awi\+mthods\Nastran Line
% 455 - 459 currently set for NASTRAN 2019, need to be adjusted for
% different versions. 

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load model parameters 
Param=eval('A321_v2');

% Update Param 
Param.FWT.Fold_eta=0.8;

if Param.FWT.Fold_eta==1
    
    Param=rmfield(Param,'FWT');
    
end

%% create run_folder
if ~isfolder(fullfile(pwd,'bin'))
   mkdir(fullfile(pwd,'bin'))
else
   delete(fullfile('bin','*'))
end
run_folder = fullfile(pwd,'bin');
%% Run analysis

% Hinge_Lock : on / off

% Run static trim analysis --------------------------------------
[FEM_full,CDi,CD0,CL,k,Aerodynamic_distribution,Load_distribution,Displacements_Res,Box_dimensions, Box_CrossSec]=Static_Trim_v1(Param, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

% Result plotting (e.g. Bending_moemnt in plane-2 (out of plane))
figure 
plot(Load_distribution.Y,Load_distribution.Moment_P2,'s-')

xlabel('Spanwise distance')
ylabel('Out of plane bending moment (Nm)')


% Run gust analysis (1mc)-----------------------------------------
Loads=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test','Mach_Number',0.78,'Altitude',36000,'Hinge_Lock','off');

% Result plotting: incremental moment along the wingspan
figure 
plot(Loads.Y(1:end-1),Loads.Wing_Delta.Max_Moment,'s-')
xlabel('Spanwise distance')
ylabel('Incremental bending moment (Nm)')

% Result plotting: incremental moment at the wing root
figure 
Gust_length=linspace(18,214,7);
plot(Gust_length,Loads.Root_Delta.Max_Moment,'s-')
hold on 
plot(Gust_length,Loads.Root_Delta.Min_Moment,'s-')
xlabel('Gust length (m)')
ylabel('Incremental bending moment (Nm)')


% Result plotting: time history incremental moment at the wing root
figure 
Time=linspace(0,2.5,201);
plot(Time,Loads.Time_Response.Root_Moment,'-','LineWidth',1)

xlabel('Time (s)')
ylabel('Incremental bending moment (Nm)')
set(gcf,'Color','w')
