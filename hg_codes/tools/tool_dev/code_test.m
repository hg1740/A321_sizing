
run_folder = 'C:\Git\A321_sizing\hg_codes\results\test_temp'; %[-], folder for exporting the NASTRAN model

Param=eval('A321_v2');

Param=rmfield(Param,'FWT');

Safty_Factor=1.5;

Yield_strength=5.2e8;


[~, Sizing_Loads, Box_dimensions, Box_CrossSec]=Sizing_Evelope(Param,run_folder);

Param_Initial=Spar_Cap_Sizing(Param, Sizing_Loads, Box_dimensions, Safty_Factor, Yield_strength);

[~, Sizing_Loads1, Box_dimensions1, Box_CrossSec1]=Sizing_Evelope(Param_Initial,run_folder);


[Area,Iyy, Izz, J]=Beam_Proeprties_v1(Box_dimensions.Inboard.Width,Box_dimensions.Inboard.Height,Param);


[Area1,Iyy1, Izz1, J1]=Beam_Proeprties_v1(Box_dimensions1.Inboard.Width,Box_dimensions1.Inboard.Height,Param_Initial);





% 
% Param1=Param;
% 
% Param1.Wing.SparCap_Thickness=Param.Wing.SparCap_Thickness*0.1;