function [Area,Iyy, Izz, J]=Beam_Proeprties_v1(Width,Height,Param,varargin)


p=inputParser;
addParameter(p, 'FWT', false, @(x)validateattributes(x, {'logical'},  {'scalar'}));
addParameter(p, 'Wing', false, @(x)validateattributes(x, {'logical'},  {'scalar'}));
parse(p,varargin{:})

%Shorthand variables
w=Width;
h=Height;


if p.Results.Wing==true && p.Results.FWT==false
    
    lcw=Param.Wing.CapEta_width*w; % spar cap top length
    lch=Param.Wing.CapEta_height*h; % spar cap web length
    
    t_cap=Param.Wing.SparCap_Thickness;
    t_web=Param.Wing.SparWeb_Thickness;
    
    
    % skin-stringer panel params
    % Note stringer pitch is assumed to be equal to the effective width of the ss panel.
    % Hence be=b.
    t_skin=Param.Wing.Skin_String.Skin_Thickness;
    be=Param.Wing.Skin_String.Effective_Width;
    b=Param.Wing.Skin_String.Stg_Pitch;
    
    ba=Param.Wing.Skin_String.StrgGround_Width;
    bw=Param.Wing.Skin_String.Strg_Depth;
    bf=Param.Wing.Skin_String.StrgFlange_Width;
    
    ta=Param.Wing.Skin_String.StrgThickness_Ground;
    tw=Param.Wing.Skin_String.StrgThickness_Web;
    tf=Param.Wing.Skin_String.StrgThickness_Flange;
    
elseif p.Results.Wing==false && p.Results.FWT==true
    
    
    lcw=Param.FWT.CapEta_width*w; % spar cap top length
    lch=Param.FWT.CapEta_height*h; % spar cap web length
    
    t_cap=Param.FWT.SparCap_Thickness;
    t_web=Param.FWT.SparWeb_Thickness;
    
    
    % skin-stringer panel params
    % Note stringer pitch is assumed to be equal to the effective width of the ss panel.
    % Hence be=b.
    t_skin=Param.FWT.Skin_String.Skin_Thickness;
    be=Param.FWT.Skin_String.Effective_Width;
    b=Param.FWT.Skin_String.Stg_Pitch;
    
    ba=Param.FWT.Skin_String.StrgGround_Width;
    bw=Param.FWT.Skin_String.Strg_Depth;
    bf=Param.FWT.Skin_String.StrgFlange_Width;
    
    ta=Param.FWT.Skin_String.StrgThickness_Ground;
    tw=Param.FWT.Skin_String.StrgThickness_Web;
    tf=Param.FWT.Skin_String.StrgThickness_Flange;
    
    
elseif p.Results.FWT==true && p.Results.Wing==true
    
    
    lcw=Param.Wing.CapEta_width*w; % spar cap top length
    lch=Param.Wing.CapEta_height*h; % spar cap web length
    
    t_cap=[Param.Wing.SparCap_Thickness,Param.FWT.SparCap_Thickness(2:end)];
    t_web=[Param.Wing.SparWeb_Thickness,Param.FWT.SparWeb_Thickness(2:end)];
    
    
    % skin-stringer panel params
    % Note stringer pitch is assumed to be equal to the effective width of the ss panel.
    % Hence be=b.
    t_skin=[Param.Wing.Skin_String.Skin_Thickness, Param.FWT.Skin_String.Skin_Thickness(2:end)];
    be=[Param.Wing.Skin_String.Effective_Width, Param.FWT.Skin_String.Effective_Width(2:end)];
    b=[Param.Wing.Skin_String.Stg_Pitch,Param.FWT.Skin_String.Stg_Pitch(2:end)];
    
    ba=[Param.Wing.Skin_String.StrgGround_Width, Param.FWT.Skin_String.StrgGround_Width(2:end)];
    bw=[Param.Wing.Skin_String.Strg_Depth, Param.FWT.Skin_String.Strg_Depth(2:end)];
    bf=[Param.Wing.Skin_String.StrgFlange_Width,Param.FWT.Skin_String.StrgFlange_Width(2:end)];
    
    ta=[Param.Wing.Skin_String.StrgThickness_Ground,Param.FWT.Skin_String.StrgThickness_Ground(2:end)];
    tw=[Param.Wing.Skin_String.StrgThickness_Web,Param.FWT.Skin_String.StrgThickness_Web(2:end)];
    tf=[Param.Wing.Skin_String.StrgThickness_Flange,Param.FWT.Skin_String.StrgThickness_Flange(2:end)];
    
 
    
end


% Area of induvidual stringer at each section.
Strg_Area=ba.*ta + bf.*tf + bw.*tw;


%% out of plane moment of inertia Iyy

% cap
% Iyy_sparcap_=lcw.*t_cap.^3/12 + lcw.*t_cap.*(0.5*h).^2 + t_cap.*lch.^3/12 + t_cap.*lch.*(0.5*h - 0.5*lch).^2;
Iyy_sparcap_= lcw.*t_cap.^3/12 + lcw.*t_cap.*(0.5*h).^2 + t_cap.*lch.^3/12 + t_cap.*lch.*(0.5*h - 0.5*lch).^2;

Iyy_sparcap=Iyy_sparcap_*4;


% web
Iyy_web_=t_web.*h.^3/12;

Iyy_web=Iyy_web_*2;


% skin-stringer panel 

% skin
Iyy_skin=w.*t_skin.^3/12 + w.*t_skin.*(0.5*h).^2;

% stringer
seg1=ba.*ta.^3/12 + ba.*ta.*(0.5*h).^2;

seg2=tw.*bw.^3/12 + tw.*bw.*(0.5*h-0.5*bw).^2;

seg3=bf.*tf.^3/12 + bf.*tf.*(0.5*h-bf).^2;

Iyy_strg_= seg1+seg2+seg3;

Num_strg=round(w./be);

Iyy_SkinStrg=Iyy_strg_.*Num_strg*2 + Iyy_skin*2;

% total Iyy

Iyy=Iyy_sparcap + Iyy_web + Iyy_SkinStrg;


%% in-plane moment of inertia Izz

% cap
Izz_sparcap_=lch.*t_cap.^3/12 + t_cap.*lch.*(0.5*w).^2 + t_cap.*lcw.^3/12 + t_cap.*lcw.*(0.5*w - 0.5*lcw).^2;

Izz_sparcap=Izz_sparcap_*4;


% web
Izz_web_=h.*t_web.^3/12 + t_web.*h.*(0.5*w).^2;

Izz_web=Izz_web_*2;


% skin-stringer panel 

% skin
Izz_skin_=t_skin.*w.^3/12;

Izz_skin=Izz_skin_*2;


% stringer 
seg1_zz=ta.*ba.^3/12;

seg2_zz=bw.*tw.^3/12;

seg3_zz=tf.*bf.^3/12 + bf.*tf.*(0.5*bf).^2;

Izz_strg_= seg1_zz+seg2_zz+seg3_zz;


Izz_strg=zeros(1,length(w));

for i=1:length(length(w))
    
    if mod(Num_strg(i),2)==0
        
        sp=be(i);
        
        offset=0.5*sp:sp:w/2;
        Izz_strg(i)=(Num_strg(i)*Izz_strg_(i) + Strg_Area(i)*(sum(offset.^2)))*2;
        
    elseif mod(Num_strg(i),2)==1
        
        sp=be(i);
        
        offset=0:sp:w/2;
        Izz_strg(i)=(Num_strg(i)*Izz_strg_(i) + Strg_Area(i)*(sum(offset.^2)))*2;
        
    end
    
end

Izz=Izz_strg+Izz_web+Izz_skin+Izz_sparcap;


%% Torsional moment of area 

% sapr web + skin 
J=4*(h.*w).^2./(4*lcw./t_cap + 4*lch./t_cap + 2*w./t_skin + 2*h./t_web);

% J=4*(h.*w).^2./(4*lcw./t_cap + 4*lch./t_cap + 2*w./t_skin + 2*h./t_web);

%% Areas of beam cross section 

Area_cap=(lcw + lch).*t_cap;
Area_skin=w.*t_skin;
Area_web=h.*t_web;

Area.spar_cap=4*Area_cap;
Area.spar_web=2*Area_web;
Area.skin=2*Area_skin;
Area.Strg=2*Strg_Area.*Num_strg;

Area.cross_section=Area.spar_cap + Area.spar_web + Area.skin + Area.Strg;



end