
    
Param=eval('A321');

% Select fold length 
Param.FWT.Fold_eta=1;


% Update Wing properties
Param.Wing.AR=10.172; % 10.172 for A321

Param.Wing.TotalArea=126; % 126m^2 for A321

Param.Wing.Dihedral=0;

Param.FWT.fold_angle=0;

Param.FWT.Flare_angle=30;


%---------------------------------------------------------

if Param.FWT.Fold_eta==1
    
    Param=rmfield(Param,'FWT');
    
end

% Generate wing planform properties 

if isfield(Param,'FWT') 
    [Geo_Wing, Geo_FWT]= Wing_Gen_V1(Param);
else
    [Geo_Wing]= Wing_Gen_V1(Param);
end

% Update Wing planform 

Param.Wing.Root_Chord=Geo_Wing.Root_Chord;

Param.Wing.TE_Sweep1=Geo_Wing.TE_Sweep1;

Param.Wing.TE_Sweep2=Geo_Wing.TE_Sweep2;

Param.Wing.Span=Geo_Wing.Span;

Param.Wing.Semi_Span=(Param.Wing.Span-4)/2;

Param.Wing.HalfArea=(Param.Wing.TotalArea-Param.Wing.Root_Chord*4)/2;


if isfield(Param,'FWT') 
    
    % update FWT planeform 
    
    Param.FWT.Root_Chord=Geo_FWT.Root_Chord;
    
    Param.FWT.Root_Height=Geo_FWT.Root_Height;
    
    Param.FWT.Tip_Chord=Geo_FWT.Tip_Chord;
    
    Param.FWT.Tip_Height=Geo_FWT.Tip_Height;
    
    
end


%% Model check

if isfield(Param,'FWT')
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v1(Param);
    
    
else
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
    
end

draw(FEM_full);

set(gcf,'Color','w')