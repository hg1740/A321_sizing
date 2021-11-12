AR10_eta100=load('NRes_A126_AR10_Eta_100_Model');
AR10_eta90=load('NRes_A126_AR10_Eta_90_Model');
AR10_eta80=load('NRes_A126_AR10_Eta_80_Model');
AR10_eta70=load('NRes_A126_AR10_Eta_70_Model');
AR10_eta60=load('NRes_A126_AR10_Eta_60_Model');

AR10_eta100_Param=AR10_eta100.Param;
AR10_eta90_Param=AR10_eta90.Param;
AR10_eta80_Param=AR10_eta80.Param;
AR10_eta70_Param=AR10_eta70.Param;
AR10_eta60_Param=AR10_eta60.Param;


[Load_Distribution_eta100, Sizing_Loads_eta100, Box_dimensions_eta100, Box_CrossSec_eta100]=Sizing_Evelope(AR10_eta100_Param,run_folder);
[Load_Distribution_eta90, Sizing_Loads_eta90, Box_dimensions_eta90, Box_CrossSec_eta90]=Sizing_Evelope(AR10_eta90_Param,run_folder);
[Load_Distribution_eta80, Sizing_Loads_eta80, Box_dimensions_eta80, Box_CrossSec_eta80]=Sizing_Evelope(AR10_eta80_Param,run_folder);
[Load_Distribution_eta70, Sizing_Loads_eta70, Box_dimensions_eta70, Box_CrossSec_eta70]=Sizing_Evelope(AR10_eta70_Param,run_folder);
[Load_Distribution_eta60, Sizing_Loads_eta60, Box_dimensions_eta60, Box_CrossSec_eta60]=Sizing_Evelope(AR10_eta60_Param,run_folder);

figure 

% plot(Sizing_Loads_eta100.Y, Sizing_Loads_eta100.Moment_P2,'s-')
% hold on 
% plot(Sizing_Loads_eta90.Y, Sizing_Loads_eta90.Moment_P2,'s-')
% hold on 
% plot(Sizing_Loads_eta80.Y, Sizing_Loads_eta80.Moment_P2,'s-')
% hold on 
plot(Sizing_Loads_eta70.Y, Sizing_Loads_eta70.Moment_P2,'s-')
hold on 
plot(Sizing_Loads_eta60.Y, Sizing_Loads_eta60.Moment_P2,'v-')


figure 
plot(Load_Distribution_eta60.Y, Load_Distribution_eta60.Case1.Moment_P2,'s-')
hold on 
plot(Load_Distribution_eta60.Y, Load_Distribution_eta60.Case3.Moment_P2,'s-')
hold on 

plot(Load_Distribution_eta70.Y, Load_Distribution_eta70.Case1.Moment_P2,'v-')
hold on 
plot(Load_Distribution_eta70.Y, Load_Distribution_eta70.Case3.Moment_P2,'v-')


figure 
plot(Load_Distribution_eta60.Y, Load_Distribution_eta60.Case1.Shear_P2,'s-')
hold on 
plot(Load_Distribution_eta60.Y, Load_Distribution_eta60.Case3.Shear_P2,'s-')
hold on 

plot(Load_Distribution_eta70.Y, Load_Distribution_eta70.Case1.Shear_P2,'v-')
hold on 
plot(Load_Distribution_eta70.Y, Load_Distribution_eta70.Case3.Shear_P2,'v-')




figure 

plot(AR10_eta70_Param.Y,[AR10_eta70_Param.Wing.SparCap_Thickness,AR10_eta70_Param.FWT.SparCap_Thickness(2:end)],'s-');
hold on
plot(AR10_eta60_Param.Y,[AR10_eta60_Param.Wing.SparCap_Thickness,AR10_eta60_Param.FWT.SparCap_Thickness(2:end)],'v-');
















