function Param_Update=Spar_Cap_Sizing(Param, Loads, Box_dimensions, Safty_Factor, Yield_strength)


%% Material properties

    E=70e9;
    sigma_Y=Yield_strength;
    
    
 %% stress calculation   
 
 % loads extraction
 
 M_P2=Safty_Factor*Loads.Moment_P2;
 
 % box dimensions
 h=Box_dimensions.Inboard.Height';
 w=Box_dimensions.Inboard.Width';
 
 
 % spar cap 
 lcw=0.1*w; % spar cap top length
 lch=0.1*h; % spar cap web length
 
 
 % Cap thickness
 I_ref=M_P2.*(0.5*h)./(4*sigma_Y);
 
 Spar_Cap_Thickness= I_ref./(lcw.*(0.5*h).^2 + lch.^3/12 + lch.*(0.5*h-0.5*lch).^2);
 
 
 % Update Param
 
 if isfield(Param,'FWT')
     
     Param_Update=Param;
     
     Param_Update.Wing.SparCap_Thickness= Spar_Cap_Thickness(1:25)';
     
     Param_Update.Wing.SparWeb_Thickness=1e-5*ones(1,25);
     
     
     Param_Update.Wing.Skin_String.Skin_Thickness=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.Effective_Width=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.Stg_Pitch=Param.Wing.Skin_String.Effective_Width;
     
     Param_Update.Wing.Skin_String.Strg_Depth=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgFlange_Width=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgGround_Width=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgThickness_Ground=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgThickness_Web=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgThickness_Flange=1e-5*ones(1,25);
     
     
     
     Param_Update.FWT.SparCap_Thickness= Spar_Cap_Thickness(25:35)';
     
     Param_Update.FWT.SparWeb_Thickness=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.Skin_Thickness=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.Effective_Width=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.Stg_Pitch=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.Strg_Depth=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.StrgFlange_Width=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.StrgGround_Width=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.StrgThickness_Ground=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.StrgThickness_Web=1e-5*ones(1,11);
     
     Param_Update.FWT.Skin_String.StrgThickness_Flange=1e-5*ones(1,11);
      
     
 else
     
     
     Param_Update=Param;
     
     Param_Update.Wing.SparCap_Thickness= Spar_Cap_Thickness(1:25)';
     
     Param_Update.Wing.SparWeb_Thickness=0.002*ones(1,25);
     
     
     Param_Update.Wing.Skin_String.Skin_Thickness=0.002*ones(1,25);
     
     Param_Update.Wing.Skin_String.Effective_Width=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.Stg_Pitch=Param.Wing.Skin_String.Effective_Width;
     
     Param_Update.Wing.Skin_String.Strg_Depth=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgFlange_Width=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgGround_Width=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgThickness_Ground=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgThickness_Web=1e-5*ones(1,25);
     
     Param_Update.Wing.Skin_String.StrgThickness_Flange=1e-5*ones(1,25);
       
     
 end
 
 

end