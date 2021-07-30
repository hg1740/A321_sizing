% x=[-17.9,-6.29,-2,2,6.29,17.9];
% 
% Chord=[1.35,3.8,6,6,3.8,1.35];
% 
% Chord_Square=Chord.^2;
% 
% A1 = trapz(x, Chord);
% 
% A2 = trapz(x, Chord_Square);
% 
% Chord_Position=interp1([6,3.8],[0 4.29],4.3);


%% Load baseline AC properties

Param=eval('A321');

% Update Wing properties
Param.Wing.AR=10.172; % 10.172 for A321

Param.Wing.TotalArea=126; 

% Update layuot properties
Param.Layout.Horizontal_Tail_Position=41;

if isfield(Param,'FWT')
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v1(Param);
    
    
else
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
    
end

draw(FEM_full);

set(gcf,'Color','w')