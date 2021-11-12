function Param_update=Thickness_Update(Param,Constraint_uppers,Internal_Stresse, Safty_Factor, Yield_strength)


skin_coefficient=step_size(Constraint_uppers.Skin_Constraints_Upper);

spar_cap_coefficient=step_size(Constraint_uppers.Spar_cap_Constraints_Upper);

spar_web_coefficient=step_size(Constraint_uppers.Spar_web_Constraints_Upper);

strg_coefficient=step_size(Constraint_uppers.Strg_Constraints_Upper);


if numel(skin_coefficient) > 25
    
    % update spar cap thickness
    Param.Wing.SparCap_Thickness=Param.Wing.SparCap_Thickness.*spar_cap_coefficient(1:25);
    Param.FWT.SparCap_Thickness=Param.FWT.SparCap_Thickness.*spar_cap_coefficient(26:35);
            
    % update spar web thickness
    Param.Wing.SparWeb_Thickness=Param.Wing.SparWeb_Thickness.*spar_web_coefficient(1:25);
    Param.FWT.SparWeb_Thickness=Param.FWT.SparWeb_Thickness.*spar_web_coefficient(26:35);
    
    % update skin thickness
    Param.Wing.Skin_Thickness=Param.Wing.Skin_Thickness.*skin_coefficient(1:25);
    Param.FWT.Skin_Thickness=Param.FWT.Skin_Thickness.*skin_coefficient(26:35);
    
    % update stringer area
    Param.Wing.Stringer_Area=Param.Wing.Stringer_Area.*strg_coefficient(1:25);
    Param.FWT.Stringer_Area=Param.FWT.Stringer_Area.*strg_coefficient(26:35);
    
    % update rib thickness
    sigma=Internal_Stresse.sparweb_Crush;
    
    Y=Internal_Stresse.Y;
    y1=Y(1:end-1);
    y2=Y(2:end);
    l=y2-y1;
    
    Rib_Y=(y1+y2)/2;
    
    crush1=sigma(1:end-1);
    crush2=sigma(2:end);
    
    Rib_crush=(crush1+crush2)/2;
    
    Rib_Thickness=Safty_Factor.*Rib_crush.*l'./ Yield_strength;
    
    Param.Wing.Rib_Thickness=Rib_Thickness;
    
    Param.Wing.Rib_Y=Rib_Y;
    
    Param.Wing.Y=Y;
       
    
else
    
    % update spar cap thickness
    Param.Wing.SparCap_Thickness=Param.Wing.SparCap_Thickness.*spar_cap_coefficient(1:25);
    
    % update spar web thickness
    Param.Wing.SparWeb_Thickness=Param.Wing.SparWeb_Thickness.*spar_web_coefficient(1:25);
    
    % update skin thickness
    Param.Wing.Skin_Thickness=Param.Wing.Skin_Thickness.*skin_coefficient(1:25);
    
    % update stringer area
    Param.Wing.Stringer_Area=Param.Wing.Stringer_Area.*strg_coefficient(1:25);
    
    % update rib thickness
    sigma=Internal_Stresse.sparweb_Crush;
    
    Y=Internal_Stresse.Y;
    y1=Y(1:end-1);
    y2=Y(2:end);
    l=y2-y1;
    
    Rib_Y=(y1+y2)/2;
    
    crush1=sigma(1:end-1);
    crush2=sigma(2:end);
    
    Rib_crush=(crush1+crush2)/2;
    
    Rib_Thickness=Safty_Factor.*Rib_crush.*l'./ Yield_strength;
    
    Param.Wing.Rib_Thickness=Rib_Thickness;
    
    Param.Wing.Rib_Y=Rib_Y;
    
    Param.Wing.Y=Y;
    
   
end


Param_update = Param;



end
