function  Param_Update=SkinStringer_Sizing(Param, Box_dimensions, Sizing_Loads, Safty_Factor, Yield_strength)



%% Material properties

E=70e9;
sigma_Y=Yield_strength;

%% Wing-box cross sectional parameters

if isfield(Param,'FWT')
    
    NumSec=35;
    
else
    
    NumSec=25;
    
end



% cross sectional properties
if isfield(Param,'FWT')
    [~,IYY, ~, ~]=Beam_Proeprties_v1(Box_dimensions.Inboard.Width,Box_dimensions.Inboard.Height,Param,'FWT',true,'Wing',true);
else
    [~,IYY, ~, ~]=Beam_Proeprties_v1(Box_dimensions.Inboard.Width,Box_dimensions.Inboard.Height,Param,'FWT',false,'Wing',true);
end

% initialise param

t_skin=zeros(1,NumSec);
be=zeros(1,NumSec);
b=zeros(1,NumSec);

ba=zeros(1,NumSec);
bw=zeros(1,NumSec);
bf=zeros(1,NumSec);

ta=zeros(1,NumSec);
tw=zeros(1,NumSec);
tf=zeros(1,NumSec);


for jj=1:NumSec
    
    % External loads
    M_P2=Safty_Factor*Sizing_Loads.Moment_P2(jj);
    T=Safty_Factor*Sizing_Loads.Torque(jj);
    S_P2=Safty_Factor*Sizing_Loads.Shear_P2(jj);
    
    % box dimensions
    hs=Box_dimensions.Inboard.Height(jj);
    ws=Box_dimensions.Inboard.Width(jj);
    
    % bending stiffness
    Iyy=IYY(jj);
    
    % rib pitch
    Lr=Param.Wing.Rib_Pitch;
    
    % effective length
    c=1.5;
    L=Lr/sqrt(c);
    
    L_inch=39.37*L; % (inch)
    
    % bending stress
    sigma=0.5*M_P2*hs/Iyy;
    
    sigma_psi=sigma*0.0001450377; % (psi)
    
    % calculating skin thickness
    % Note factor of 1.5 applied as the area ratio Ast/Askin=0.5; 
    % Total force act on the skin-stringer panel = sigma * (Ask+ 0.5*Askin)
    t_skin(jj)=1.5*sigma_psi*L_inch/(3000^2);
    
    % N load intensity (load per unit lenth)
    N=t_skin(jj).*sigma_psi;
 
    Fe=2000*(N/L_inch)^0.5;
    
    % effective width
    if Fe < 0.0001450377*sigma_Y
        be_t=-7.743e-6*(Fe/1000)^4 + 0.0006387*(Fe/1000)^3 + 0.007084*(Fe/1000)^2 - 1.966*(Fe/1000) + 76.83;
        
    else
        
        % TODO: this part need to update
        % should pop up a warning message instead ??
        
        be_t=10;
        
    end
    
    % skin-stringer panel effective width
    be(jj)=be_t*t_skin(jj); %(inch)
    
    % stringer pitch = effective width
    b(jj)=be(jj); %(inch)
    
    % ground thickness
    ta(jj)=0.7*t_skin(jj); %(inch)
    
    % ground width
    ba(jj)=ta(jj)*9.35; %(inch)
    
    % stringer area
    A_st=0.5*b(jj)*t_skin(jj); %(inch^2)
    
    % stringer depth
    bw(jj)=(be_t*((A_st-2*ba(jj)*ta(jj))/1.327))^0.5; %(inch)
    
    % stringer web thickness
    tw(jj)=bw(jj)/be_t; %(inch)
    
    % stringer flange width
    bf(jj)=0.327*bw(jj); %(inch)
    
    % stringer flange thickness
    tf(jj)=tw(jj); %(inch)
    
    
end



% Update Results


Param_Update=Param;


if isfield(Param,'FWT')
    
    
    % In board
    Param_Update.Wing.Skin_String.Skin_Thickness=t_skin(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.Effective_Width=be(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.Stg_Pitch=b(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.Strg_Depth=bw(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.StrgFlange_Width=bf(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.StrgGround_Width=ba(1:25)*2*0.0254; 
    
    Param_Update.Wing.Skin_String.StrgThickness_Ground=ta(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.StrgThickness_Web=tw(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.StrgThickness_Flange=tf(1:25)*0.0254;
    
    
    % FWT
    Param_Update.FWT.Skin_String.Skin_Thickness=t_skin(25:35)*0.0254;
    
    Param_Update.FWT.Skin_String.Effective_Width=be(25:35)*0.0254;
    
    Param_Update.FWT.Skin_String.Stg_Pitch=b(25:35)*0.0254;
    
    Param_Update.FWT.Skin_String.Strg_Depth=bw(25:35)*0.0254;
    
    Param_Update.FWT.Skin_String.StrgFlange_Width=bf(25:35)*0.0254;
    
    Param_Update.FWT.Skin_String.StrgGround_Width=ba(25:35)*0.0254;
    
    Param_Update.FWT.Skin_String.StrgThickness_Ground=ta(25:35)*0.0254;
    
    Param_Update.FWT.Skin_String.StrgThickness_Web=tw(25:35)*0.0254;
    
    Param_Update.FWT.Skin_String.StrgThickness_Flange=tf(25:35)*0.0254;
    
    
    
else
    
    
    Param_Update.Wing.Skin_String.Skin_Thickness=t_skin(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.Effective_Width=be(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.Stg_Pitch=b(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.Strg_Depth=bw(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.StrgFlange_Width=bf(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.StrgGround_Width=ba(1:25)*2*0.0254; 
    
    Param_Update.Wing.Skin_String.StrgThickness_Ground=ta(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.StrgThickness_Web=tw(1:25)*0.0254;
    
    Param_Update.Wing.Skin_String.StrgThickness_Flange=tf(1:25)*0.0254;
    
    
end


end





%     % Spar
%     t_cap1=Param.Wing.SparCap_Thickness;
%     t_web1=Param.Wing.SparWeb_Thickness;
%     
%     % skin-stringer panel
%     t_skin1=Param.Wing.Skin_String.Skin_Thickness;
%     be1=Param.Wing.Skin_String.Effective_Width;
%     b1=Param.Wing.Skin_String.Stg_Pitch;
%     
%     ba1=Param.Wing.Skin_String.StrgGround_Width;
%     bw1=Param.Wing.Skin_String.Strg_Depth;
%     bf1=Param.Wing.Skin_String.StrgFlange_Width;
%     
%     ta1=Param.Wing.Skin_String.StrgThickness_Ground;
%     tw1=Param.Wing.Skin_String.StrgThickness_Web;
%     tf1=Param.Wing.Skin_String.StrgThickness_Flange;
% 
%     
% 
%     % FWT spar
%     t_cap2=Param.FWT.SparCap_Thickness;
%     t_web2=Param.FWT.SparWeb_Thickness;
%     
%     % FWT skin-stringer panel
%     t_skin2=Param.FWT.Skin_String.Skin_Thickness;
%     be2=Param.FWT.Skin_String.Effective_Width;
%     b2=Param.FWT.Skin_String.Stg_Pitch;
%     
%     ba2=Param.FWT.Skin_String.StrgGround_Width;
%     bw2=Param.FWT.Skin_String.Strg_Depth;
%     bf2=Param.FWT.Skin_String.StrgFlange_Width;
%     
%     ta2=Param.FWT.Skin_String.StrgThickness_Ground;
%     tw2=Param.FWT.Skin_String.StrgThickness_Web;
%     tf2=Param.FWT.Skin_String.StrgThickness_Flange;
%     
%     % Combine
%     t_cap=[t_cap1,t_cap2(2:end)];
%     t_web=[t_web1,t_web2(2:end)];
%     
%     
%     t_skin=[t_skin1,t_skin2(2:end)];
%     be=[be1,be2(2:end)];
%     b=[b1,b2(2:end)];
%     
%     ba=[ba1,ba2(2:end)];
%     bw=[bw1,bw2(2:end)];
%     bf=[bf1,bf2(2:end)];
%     
%     ta=[ta1,ta2(2:end)];
%     tw=[tw1,tw2(2:end)];
%     tf=[tf1,tf2(2:end)];