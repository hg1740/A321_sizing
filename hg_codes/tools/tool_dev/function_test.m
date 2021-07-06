Param=eval('A321');

% update FWT span 
Param.FWT.Fold_eta=0.8;

% run_folder='C:\Git\A321_sizing\hg_codes\results\test_temp';

if ~isfolder(fullfile(pwd,'bin'))
   mkdir(fullfile(pwd,'bin'))
else
   delete(fullfile('bin','*'))
end
run_folder = fullfile(pwd,'bin');

[CDi,CD0,CL,k,Distribution,Load_distribution, Displacements_Res]=Static_Trim_v1(Param, run_folder, 'Load_Factor',2.5,'File_Name','pull_up_test','Hinge_Lock','off');

Loads=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','off');
Loads_on=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_test','March_Number',0.78,'Altitude',36000,'Hinge_Lock','on');

figure 
plot(Load_distribution.Y,Load_distribution.Moment_P2,'s-')


figure 
plot(Loads.Y(1:end-1),Loads.Wing_Delta.Max_Moment,'s-')
hold on 
plot(Loads_on.Y(1:end-1),Loads_on.Wing_Delta.Max_Moment,'v-')