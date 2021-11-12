%A321

% In this version of the Param, spar wiill be sbilitted into spar cap and
% web. The skin and stringer panel are sized as a whole (dependent on each other)

function Param=A321_v2

 
Param.Wing.AR=10.172;

Param.Wing.Root_Chord=6;

% spar
Param.Wing.CapEta_width=0.1;

Param.Wing.CapEta_height=0;

Param.Wing.SparCap_Thickness=0.02*ones(1,25);

Param.Wing.SparWeb_Thickness=0.005*ones(1,25);

 
% skin-stringer panels 
Param.Wing.Skin_String.Skin_Thickness=0.01*ones(1,25);

Param.Wing.Skin_String.Effective_Width=0.01*ones(1,25);

Param.Wing.Skin_String.Stg_Pitch=Param.Wing.Skin_String.Effective_Width;

Param.Wing.Skin_String.Strg_Depth=0.01*ones(1,25);

Param.Wing.Skin_String.StrgFlange_Width=0.01*ones(1,25);

Param.Wing.Skin_String.StrgGround_Width=0.01*ones(1,25);

Param.Wing.Skin_String.StrgThickness_Ground=0.001*ones(1,25);

Param.Wing.Skin_String.StrgThickness_Web=0.001*ones(1,25);

Param.Wing.Skin_String.StrgThickness_Flange=0.001*ones(1,25);






Param.Wing.LE_Sweep=27; 

Param.Wing.TE_Sweep1=0;

Param.Wing.TE_Sweep2=16.5;

Param.Wing.BeamLoc=0.4;

Param.Wing.Kink=0.27;

Param.Wing.AeroPanel_Number=10;

Param.Wing.AeroPanel_AR=2.5;

Param.Wing.Span=35.8;

Param.Wing.Semi_Span=(Param.Wing.Span-4)/2;

Param.Wing.Dihedral=5;

Param.Wing.TotalArea=126;

Param.Wing.HalfArea=(Param.Wing.TotalArea-Param.Wing.Root_Chord*4)/2;

Param.Wing.ThicknessToChord_Root=0.15;

Param.Wing.ThicknessToChord_kink=0.12;

Param.Wing.ThicknessToChord_tip=0.11;

Param.Wing.Rib_Pitch=0.6;


% Jig shape

Param.Wing.Jig_Twist=deg2rad([0,0]);

Param.Wing.Jig_Eta=[0,1];

Param.Connector.Jig_Twist=deg2rad([0,0]);

Param.Connector.Jig_Eta=[0,1];


% Folding wingtip 

Param.FWT.Fold_angle=10;

Param.FWT.Flare_angle=25;

Param.FWT.Fold_eta=0.75;

Param.FWT.Root_Chord=0.5;

Param.FWT.Root_Height=0.5;

Param.FWT.Tip_Chord=0.3;

Param.FWT.Tip_Height=0.1;

Param.FWT.Hinge_Stiffness=1e-4;

% FWT sizing variables 

% spar
Param.FWT.CapEta_width=0.1;

Param.FWT.CapEta_height=0.1;

Param.FWT.SparCap_Thickness=0.02*ones(1,11);

Param.FWT.SparWeb_Thickness=0.005*ones(1,11);

 
% skin-stringer panels 
Param.FWT.Skin_String.Skin_Thickness=0.01*ones(1,11);

Param.FWT.Skin_String.Effective_Width=0.01*ones(1,11);

Param.FWT.Skin_String.Stg_Pitch=Param.FWT.Skin_String.Effective_Width;

Param.FWT.Skin_String.Strg_Depth=0.01*ones(1,11);

Param.FWT.Skin_String.StrgFlange_Width=0.01*ones(1,11);

Param.FWT.Skin_String.StrgGround_Width=0.01*ones(1,11);

Param.FWT.Skin_String.StrgThickness_Ground=0.001*ones(1,11);

Param.FWT.Skin_String.StrgThickness_Web=0.001*ones(1,11);

Param.FWT.Skin_String.StrgThickness_Flange=0.001*ones(1,11);


% Jig shape

Param.FWT.Jig_Twist=deg2rad([0,0]);

Param.FWT.Jig_Eta=[0,1];


% Layout: Position parameters

Param.Layout.Fuselage_Length=45;

Param.Layout.Fuselage_Width=4;

Param.Layout.Wing_Position=20;

Param.Layout.Engine_Position=4.29;

Param.Layout.Horizontal_Tail_Position=42;

Param.Layout.Vertical_Tail_Position=41;

% Mass case

Param.Masses.Hinge=200; %kg

Param.Masses.Secondary_Mass=1200;

Param.Masses.Fuselage_Structure=25000;

Param.Masses.Payload=25000;

Param.Masses.Engine=7362/2;

Param.Masses.Pylon=1239/2;

Param.Masses.MTOW=93500;

Param.Masses.OWE=48500; 

Param.Masses.Fuselage_shell_mass = 2*pi*2*0.004*2800*44.5;

Param.Masses.Horizontal_tail=682; 

Param.Masses.Vertical_tail=522;

Param.Masses.Fuel_Capacity=32940; % L

Param.Masses.Fuel_Fraction=0.723;

Param.Masses.Fuel_Density=840; %g/L

Param.Masses.Fuel_Mass=Param.Masses.Fuel_Capacity * Param.Masses.Fuel_Fraction * Param.Masses.Fuel_Density/1000; 


Param.Material.Modulus=70e9;

Param.Material.Poisson=0.33;

Param.Material.Density=2800;


end


