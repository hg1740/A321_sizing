% 
% 
% 
% [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, varargout]=Aircraft_Models_v2(Param)
% 
% 
% a_b=[1,2,3,4,5,6];
% 
% Ks=[8.2,5.8,5.2,5.1,5,5];
% 
% figure 
% plot(a_b,Ks)


Param=eval('A321_v1');

run_folder='C:\Git\A321_sizing\hg_codes\results\example_test1';

[FEM_full,CDi,CD0,CL,k,Aerodynamic_distribution,Load_distribution,Displacements_Res,Box_dimensions, Box_CrossSec]=Static_Trim_v1(Param, run_folder, 'Load_Factor',2.5,'File_Name','Trim_analysis','Hinge_Lock','off');



Stresses=Internal_Stress_Calc(Param, Box_CrossSec, Box_dimensions, Load_distribution)