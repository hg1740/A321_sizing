function Param_Update=Spar_Sizing(Param, Box_dimensions, Loads, Safty_Factor, Yield_strength)


% This function aims to determine the spar area to ensure sufficient second
% moment of area to withstand the bending loads.


%% Material properties
E=70e9;
sigma_Y=Yield_strength;
shear_strength=sigma_Y/sqrt(3);


%% short hand variables
if isfield(Param,'FWT')
    
    tsw=[Param.Wing.SparWeb_Thickness,Param.FWT.SparWeb_Thickness(2:end)];
    
else
    
    tsw=Param.Wing.SparWeb_Thickness;
    
end

%% stress calculation

% loads extraction
M_P2=Safty_Factor*abs(Loads.Moment_P2);
T=Safty_Factor*abs(Loads.Torque);
V=Safty_Factor*abs(Loads.Shear_P2);

% box dimensions
h=Box_dimensions.Inboard.Height';
w=Box_dimensions.Inboard.Width';

% spar cap dimensions
ws=Param.Wing.CapEta_width*w;


Cap_area_Y=0.5*M_P2'./(h'*sigma_Y);

Cap_Thickness_Y=Cap_area_Y./ws';

Cap_Thickness_B=(0.5*M_P2'.*ws'./(3.6*E*h')).^(1/3);

Cap_Thickness=max([Cap_Thickness_Y;Cap_Thickness_B]);


% Update cap thickness
Param_Update=Param;

if isfield(Param, 'FWT')
    
    Param_Update.Wing.SparCap_Thickness=Cap_Thickness(1:25);
    Param_Update.FWT.SparCap_Thickness=Cap_Thickness(25:end);
    
else
    
    Param_Update.Wing.SparCap_Thickness=Cap_Thickness(1:25);
    
end

% Current second moment of area
if isfield(Param, 'FWT')
    
    [~,Iyy, ~, ~]=Beam_Proeprties_v1(w',h',Param_Update,'FWT',true,'Wing',true);
    
else
    
    [~,Iyy, ~, ~]=Beam_Proeprties_v1(w',h',Param_Update,'FWT',false,'Wing',true);
    
end


% check bending stress on the top web
Bending_Stress=M_P2'.*(0.5*h')./Iyy;


% short hand variable for spar cap 

if isfield(Param_Update,'FWT')
    
    tsc=[Param_Update.Wing.SparCap_Thickness,Param_Update.FWT.SparCap_Thickness(2:end)];
    
else
    
    tsc=Param_Update.Wing.SparCap_Thickness;
    
end


% spar web stresses
Q_skn=w'.*tsc.*(0.5*h');
Q_spar= ws'.*tsc.*(0.5*h')*2;
Q_web=(0.5*tsw.*h'.*(0.25*h'))*2;

Q=Q_skn + Q_spar + Q_web;

Shear_stress=V'.*Q./(2*Iyy.*tsw) + T'./(2*h'.*w'.*tsw);


% web critical buckling stresses
Sigma_buckling_web=21*E*(tsw./h);


% web constraint check
Constraint_web1=Shear_stress./shear_strength;

Constraint_web2=Bending_Stress./Sigma_buckling_web;

Constraint_web=max([Constraint_web1;Constraint_web2]);


% Update thickness

SparWeb_adjust=step_size(Constraint_web);

if isfield(Param,'FWT')
    
    Param_Update.Wing.SparWeb_Thickness=Param.Wing.SparWeb_Thickness.*SparWeb_adjust(1:25);
    
    Param_Update.FWT.SparWeb_Thickness=Param.FWT.SparWeb_Thickness.*SparWeb_adjust(25:end);
    
else
    
    Param_Update.Wing.SparWeb_Thickness=Param.Wing.SparWeb_Thickness.*SparWeb_adjust;
    
end


end




