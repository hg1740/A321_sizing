function [Rib_Thickness,Rib_Weight]=Rib_Sizing(Param,Sizing_Loads,Safty_Factor)

E=70e9;
yield_strength=5.2e8;
rho=2810;

if isfield(Param,'FWT')
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT_R]=Aircraft_Models_v3(Param);
    
else
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v3(Param);
end



% box dimensions
hs=Box_dimensions.Inboard.Height;
ws=Box_dimensions.Inboard.Width;

% second moment of area
Iyy=Box_CrossSec.Iyy;

% loads
M_P2=Safty_Factor*Sizing_Loads.Moment_P2;
T=Safty_Factor*Sizing_Loads.Torque;
S_P2=Safty_Factor*Sizing_Loads.Shear_P2;


% spar cap thickness
if isfield(Param,'FWT')
    t_sc=[Param.Wing.SparCap_Thickness,Param.FWT.SparCap_Thickness(2:end)];
else
    t_sc=Param.Wing.SparCap_Thickness;   
end

% crushing stress
N_crush=0.5*M_P2'.^2.*hs.*t_sc./(Iyy.^2*E);

N1=N_crush(1:end-1);
N2=N_crush(2:end);

N=(N1+N2)/2;

% rib pitch 

% Y=Y_eta;
% Y=Sizing_Loads.Y;

if isfield(Param,'FWT')
    
    Y_fwt=Y_eta(end) + linspace(FWT_R.RData(1),FWT_R.RData(2),11);
    
    Y=[Y_eta,Y_fwt(2:end)];
    
else
    
    Y=Y_eta;
    
end


Y1=Y(1:end-1);
Y2=Y(2:end);

Spacing=Y2-Y1;

% Chord width 
w1=ws(1:end-1);
w2=ws(2:end);

wm=(w1+w2)/2;

% rib depth
d1=hs(1:end-1);
d2=hs(2:end);

dm=(d1+d2)/2;

% crushing force
Crush_F=N.*Spacing.*wm;

% critical thickness for yielding 
t_y=Crush_F./(wm*yield_strength);

% critical thickness for wide column buckling
A=pi^2*E/(12*(1-0.33^2));

t_cb=(Crush_F.*dm./(A*wm)).^(1/3);


% k=10;
% t_b=(Crush_F.*wm./(k*E)).^(1/3);
% 
Rib_Thickness=max([t_y;t_cb]);

% Rib weight 
Rib_Weight=sum(1.6*wm.*dm.*Rib_Thickness*rho);


end






