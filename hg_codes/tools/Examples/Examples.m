
%% Load model parameters 
Param=eval('A321');

% Update Param 
Param.FWT.Fold_eta=0.8;

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
[CDi,CD0,CL,k,Distribution,Load_distribution, Displacements_Res]=Static_Trim_v1(Param, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');

% Result plotting (e.g. Bending_moemnt in plane-2 (out of plane))
figure 
plot(Load_distribution.Y,Load_distribution.Moment_P2,'s-')

xlabel('Spanwise distance')
ylabel('Out of plane bending moment (Nm)')


% Run gust analysis (1mc)-----------------------------------------
Loads=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');

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



