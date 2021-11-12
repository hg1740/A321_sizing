function Param_Update=Spar_Sizing(Param, Box_dimensions, Loads, Safty_Factor, Yield_strength)


% This function aims to adjust spar thickness based on the initial guess to
% ensure a correct overall second moment of area.


%% Material properties

E=70e9;
sigma_Y=Yield_strength;

shear_strength=sigma_Y/sqrt(3);


%% stress calculation

% loads extraction

M_P2=Safty_Factor*abs(Loads.Moment_P2);
T=Safty_Factor*abs(Loads.Torque);
V=Safty_Factor*abs(Loads.Shear_P2);

% box dimensions
h=Box_dimensions.Inboard.Height';
w=Box_dimensions.Inboard.Width';

% spar cap dimensions
hs=Param.Wing.CapEta_height*h;
ws=Param.Wing.CapEta_width*w;


% calculating the required second moment of area
I_ref=M_P2.*(0.5*h)/sigma_Y;


% initiate the loop
indicator=1;

while indicator > 0.05
    

    
    % calculate current stiffness 
    [~,Iyy, ~, ~]=Beam_Proeprties_v1(w',h',Param);
    
    
    % skin and spar cap thickness
    t_skin=Param.Wing.Skin_String.Skin_Thickness;
    
    t_spr=Param.Wing.SparCap_Thickness;
    

    % fist moment of areas   
    Q_skn=w'.*t_skin.*(0.5*h');
    Q_spar= (ws'.*t_spr.*(0.5*h') + hs'.*t_spr.*(0.5*h'-0.5*hs'))*2;
    
    Q0=Q_skn+Q_spar;
    
    % spar web thickness calculation 
    t_web=(V'.*Q0./(2*Iyy) + T'./(2*h'.*w'))./(shear_strength - 0.125*V'.*h'.^2./(2*Iyy));
    
    
    % update spar web
    Param.Wing.SparWeb_Thickness=t_web;
    
    % update stiffness
    [~,Iyy_web, ~, ~]=Beam_Proeprties_v1(Box_dimensions.Inboard.Width,Box_dimensions.Inboard.Height,Param);
    
    % adjusting spar cap thickness based on the delta (4 spar caps)
    delta_Iyy=0.25*(I_ref'-Iyy_web);
    
    indicator=max(abs(delta_Iyy./Iyy_web));
    
    cons=ws'.*(0.5*h').^2 + hs'.*(0.5*h' - 0.5*hs').^2;
    
    dt=delta_Iyy./cons;
    
    
    % update spar cap thickness
    
    Param.Wing.SparCap_Thickness=Param.Wing.SparCap_Thickness + dt;
    
  
    
end

Param_Update=Param;



end




