function [Area,Iyy, Izz, J]=Beam_proeprties(Width,Height,SparCap_Thickness,SparWeb_Thickness, Skin_Thickness, Strg_Area, CapEta_w, CapEta_h)


% stringer pitch 
sp=0.24;

%Shorthand variables
w=Width;
h=Height;

lcw=CapEta_w*w; % spar cap top length
lch=CapEta_h*h; % spar cap web length 

t_cap=SparCap_Thickness;
t_web=SparWeb_Thickness;
t_skin=Skin_Thickness;
A_strg=Strg_Area;


%% out of plane moment of inertia Iyy

% cap
Iyy_sparcap_=lcw.*t_cap.^3/12 + lcw.*t_cap.*(0.5*h).^2 + t_cap.*lch.^3/12 + t_cap.*lch.*(0.5*h - 0.5*lch).^2;

Iyy_sparcap=Iyy_sparcap_*4;

% skin
Iyy_skin_=w.*t_skin.^3/12 + t_skin.*w.*(0.5*h).^2;

Iyy_skin=Iyy_skin_*2;

% web
Iyy_web_=t_web.*h.^3/12;

Iyy_web=Iyy_web_*2;

% stringer 
d_strg=sqrt(A_strg/0.36);
t_strg=0.12*d_strg;

seg1=d_strg.*t_strg.^3/12 + d_strg.*t_strg.*(0.5*h).^2;

seg2=t_strg.*d_strg.^3/12 + d_strg.*t_strg.*(0.5*h-0.5*d_strg).^2;

seg3=d_strg.*t_strg.^3/12 + d_strg.*t_strg.*(0.5*h-d_strg).^2;

Iyy_strg_= seg1+seg2+seg3;

Num_strg=w/sp;

Iyy_strg=Iyy_strg_.*Num_strg*2;

% total Iyy

Iyy=Iyy_sparcap + Iyy_web + Iyy_skin + Iyy_strg;


%% in-plane moment of inertia Izz

% cap
Izz_sparcap_=lch.*t_cap.^3/12 + t_cap.*lcw.*(0.5*w).^2 + t_cap.*lcw.^3/12 + t_cap.*lcw.*(0.5*w - 0.5*lcw).^2;

Izz_sparcap=Izz_sparcap_*4;

% skin
Izz_skin_=t_skin.*w.^3/12;

Izz_skin=Izz_skin_*2;

% web
Izz_web_=h.*t_web.^3/12 + t_web.*h.*(0.5*w).^2;

Izz_web=Izz_web_*2;

% stringer 
seg1_zz=t_strg.*d_strg.^3/12 + d_strg.*t_strg.*(0.5*d_strg).^2;

seg2_zz=d_strg.*t_strg.^3/12;

seg3_zz=t_strg.*d_strg.^3/12 + d_strg.*t_strg.*(0.5*d_strg).^2;

Izz_strg_= seg1_zz+seg2_zz+seg3_zz;


Izz_strg=zeros(1,length(w));

for i=1:length(length(w))
    
    if mod(Num_strg(i),2)==0
        offset=0.5*sp:sp:w/2;
        Izz_strg(i)=(Num_strg(i)*Izz_strg_(i)+3*d_strg(i)*t_strg(i)*(sum(offset.^2)))*2;
        
    elseif mod(Num_strg(i),2)==1
        offset=0:sp:w/2;
        Izz_strg(i)=(Num_strg(i)*Izz_strg_(i)+3*d_strg(i)*t_strg(i)*(sum(offset.^2)))*2;
        
    end
    
end

Izz=Izz_strg+Izz_web+Izz_skin+Izz_sparcap;



%% Torsional moment of area 

% spar carp + sapr web + skin 
J=4*(h.*w).^2./(4*lcw./t_cap + 4*lch./t_cap + 2*w./t_skin + 2*h./t_web);




%% Areas of beam cross section 

Area_cap=(lcw + lch).*t_cap;
Area_skin=w.*t_skin;
Area_web=h.*t_web;

Area.spar_cap=Area_cap;
Area.spar_Web=Area_web;
Area.skin=Area_skin;
Area.Strg=Strg_Area.*Num_strg;

Area.cross_section=Area.spar_cap*4 + Area.spar_Web*2 + Area.skin*2 + Area.Strg*2;



end