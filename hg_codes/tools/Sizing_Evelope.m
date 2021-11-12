function [Load_Distribution, Sizing_Loads, Box_dimensions, Box_CrossSec]=Sizing_Evelope(Param,run_folder)


%% Load evelope


if isfield(Param,'FWT')
    % Zero Fuel 36000 ft pull up (Hinge free)
    Param_ZF=Param;
    
    % Update Fuel mass
    Param_ZF.Masses.Fuel_Mass=1;
    
    [~,~,~,~,~,~,TrimLoad_HF1,~,Box_dimensions, Box_CrossSec]=Static_Trim_v1(Param_ZF, run_folder, 'Load_Factor',2.5,'File_Name','Pull_up36000ft','Hinge_Lock','off','Altitude',36000,'Mach_Num',0.78);
    
        
%     % 1g + severe gust (Hinge free)
%     [~,~,~,~,~,~,TrimLoad_HF2,~,~, ~]=Static_Trim_v1(Param, run_folder, 'Load_Factor',1,'File_Name','Level_3000ft','Hinge_Lock','off','Altitude',3000,'Mach_Num',0.48);
%     
%     GustLoads_HF=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_HF','Mach_Number',0.78,'Altitude',36000,'Hinge_Lock','off','Altitude',3000, 'Gust_Eta', 1);
    
    
    % Gust threshold (Hinge locked)
    [~,~,~,~,~,~,TrimLoad_HL,~,~, ~]=Static_Trim_v1(Param, run_folder, 'Load_Factor',1,'File_Name','Level_3000ft','Hinge_Lock','on','Altitude',3000,'Mach_Num',0.48);
    
    GustLoads_HL=Gust_Analysis_v1(Param,run_folder,'File_Name','gust_HL','Mach_Number',0.48,'Altitude',3000,'Hinge_Lock','on','Altitude',3000, 'Gust_Eta', 0.33);
    
    
else
    
    
    % Zero Fuel 36000 ft pull up (Hinge free)
    Param_ZF=Param;
    
    % Update Fuel mass
    Param_ZF.Masses.Fuel_Mass=1;
    
    [~,~,~,~,~,~,TrimLoad_HF1,~,Box_dimensions, Box_CrossSec]=Static_Trim_v1(Param_ZF, run_folder, 'Load_Factor',2.5,'File_Name','Pull_up36000ft','Hinge_Lock','off','Altitude',36000,'Mach_Num',0.78);
    
  
end


%% sizing load

if isfield(Param,'FWT')
    
    % ZF static
    Moment_P2_Case1=abs(TrimLoad_HF1.Moment_P2);
    Shear_P2_Case1=abs(TrimLoad_HF1.Shear_P2);
    Torque_Case1=abs(TrimLoad_HF1.Torque);
    
    % 1g + gust
%     Moment_P2_Case2=abs(TrimLoad_HF2.Moment_P2) + [GustLoads_HF.Wing_Delta.Max_Moment; 0];
%     Shear_P2_Case2=abs(TrimLoad_HF2.Shear_P2) + [GustLoads_HF.Wing_Delta.Max_Shear; 0];
%     Torque_Case2=abs(TrimLoad_HF2.Torque) + [GustLoads_HF.Wing_Delta.Max_Torque; 0];
    
    % gust threshold
    Moment_P2_Case3=abs(TrimLoad_HL.Moment_P2) + [GustLoads_HL.Wing_Delta.Max_Moment; 0];
    Shear_P2_Case3=abs(TrimLoad_HL.Shear_P2) + [GustLoads_HL.Wing_Delta.Max_Shear; 0];
    Torque_Case3=abs(TrimLoad_HL.Torque) + [GustLoads_HL.Wing_Delta.Max_Torque; 0];
    
    
    
    %% Result output
    
    Load_Distribution.Case1.Moment_P2=Moment_P2_Case1;
    Load_Distribution.Case1.Shear_P2=Shear_P2_Case1;
    Load_Distribution.Case1.Torque=Torque_Case1;
    
%     Load_Distribution.Case2.Moment_P2=Moment_P2_Case2;
%     Load_Distribution.Case2.Shear_P2=Shear_P2_Case2;
%     Load_Distribution.Case2.Torque=Torque_Case2;
    
    Load_Distribution.Case3.Moment_P2=Moment_P2_Case3;
    Load_Distribution.Case3.Shear_P2=Shear_P2_Case3;
    Load_Distribution.Case3.Torque=Torque_Case3;
    
    Load_Distribution.Y=TrimLoad_HF1.Y;
    
    Sizing_Loads.Moment_P2=max([Moment_P2_Case1,Moment_P2_Case3],[],2);
    Sizing_Loads.Shear_P2=max([Shear_P2_Case1,Shear_P2_Case3],[],2);
    Sizing_Loads.Torque=max([Torque_Case1,Torque_Case3],[],2);
    Sizing_Loads.Y=TrimLoad_HF1.Y;
    
    
else
    
    % ZF static
    Moment_P2_Case1=abs(TrimLoad_HF1.Moment_P2);
    Shear_P2_Case1=abs(TrimLoad_HF1.Shear_P2);
    Torque_Case1=abs(TrimLoad_HF1.Torque);
    

    %% Result output
    
    Load_Distribution.Case1.Moment_P2=Moment_P2_Case1;
    Load_Distribution.Case1.Shear_P2=Shear_P2_Case1;
    Load_Distribution.Case1.Torque=Torque_Case1;
    
    
    Load_Distribution.Y=TrimLoad_HF1.Y;
    
    Sizing_Loads.Moment_P2=Moment_P2_Case1;
    Sizing_Loads.Shear_P2=Shear_P2_Case1;
    Sizing_Loads.Torque=Torque_Case1;
    Sizing_Loads.Y=TrimLoad_HF1.Y;
    
end





end

