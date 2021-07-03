    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[0,0,0]; %15
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';
    
    % Num of element 
    Wingbox_right.NumBeamElem = 23;

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = 16;   
    Wingbox_right.LESweep     = [27, 27];
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [1, 16.5, 16.5];
    Wingbox_right.TESweep_eta = [0, 0.27, 1];
    Wingbox_right.RootChord   = 6;
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_right.Eta_;
    Wingbox_right.BeamLoc     = repmat(0.4, size(all_eta));

    Wingbox_right.BeamLoc_eta = all_eta;

    %Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_right = awi.model.Spar;
    FrontSpar_right.XLoc = [0.15, 0.15];
    FrontSpar_right.Eta  = [0   , 1];
    RearSpar_right = awi.model.Spar;
    RearSpar_right.XLoc = [0.65, 0.65];
    RearSpar_right.Eta  = [0   , 1];

    Wingbox_right.add([FrontSpar_right, RearSpar_right]);

    %Define internal layout
    Wingbox_right.RibPitch      = 0.65;
    Wingbox_right.StringerPitch = 0.15;

    %Make the material
    E_wing  = 70e9; %[N/m^2], typical YM of aluminium
    nu_wing = 0.333;
    rho_wing=2810;
    Mat_wing = awi.model.Material;
    Mat_wing.E  = E_wing;
    Mat_wing.Nu = nu_wing;
    Mat_wing.G  = E_wing / (2 * (1 + nu_wing));
    Mat_wing.Rho=rho_wing;
    
    Wingbox_right.Material_eta = [0, 1];
    Wingbox_right.Material     = [Mat_wing, Mat_wing];
    
    build(Wingbox_right)
    
    % Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
    Wingbox_right.AeroPanelLength=[];
    Wingbox_right.NumAeroPanel=10;
    Wingbox_right.AeroPanelAR=1;
    
    %     Wingbox_right.AeroPanelLength=0.4;
    
  
    eta_=[0,1];
    
    Wingbox_right.A   =  0.1;
    Wingbox_right.A_eta=eta_;

    Wingbox_right.I11 = 0.01;
    Wingbox_right.I11_eta=eta_;

    Wingbox_right.I22 = 0.01;
    Wingbox_right.I22_eta = eta_;

    Wingbox_right.J   = 0.01;
    Wingbox_right.J_eta= eta_;


    % Jig shape
    
    Wingbox_right.Twist = [0,0];
    Wingbox_right.Twist_eta = [0,1];
    
    build(Wingbox_right)
    
    FWT_R = insertWingFold(Wingbox_right, 'FlareAngle', 25, 'FoldAngle', 10,'EtaFold',0.4);
    FWT_R.HingeStiffness = [1e14 1e14 1e14 1e14 1e-4 1e14];
    
    
    FEM_test=convertToFE(Wingbox_right);
    
    draw(FEM_test)
    
    
    
    
    
    
    Wingbox_left = awi.model.LiftingSurface;
    Wingbox_left.Name = 'A320Wing_left';
    Wingbox_left.Origin=[Param.Layout.Wing_Position,-2,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_left.ActiveSet = 'sSet';

    %Wing dimensions
    Wingbox_left.SpanVector  = 'Y';
    Wingbox_left.Span        = -Param.Wing.Semi_Span;  
    Wingbox_left.LESweep     = [-Param.Wing.LE_Sweep, -Param.Wing.LE_Sweep];
    Wingbox_left.LESweep_eta = [0, 1];
    Wingbox_left.TESweep     = [-Param.Wing.TE_Sweep1, -Param.Wing.TE_Sweep2, -Param.Wing.TE_Sweep2];
    Wingbox_left.TESweep_eta = [0, Param.Wing.Kink, 1];
    Wingbox_left.RootChord   = Param.Wing.Root_Chord;   
    
    %Dihedral 
    Wingbox_left.Dihedral=[Param.Wing.Dihedral,Param.Wing.Dihedral];
    Wingbox_left.Dihedral_eta=[0,1];
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_left.Eta_;
    Wingbox_left.BeamLoc     = repmat(Param.Wing.BeamLoc, size(all_eta));
%     Wingbox_right.BeamLoc     = [0.34,0.4,0.4];
    Wingbox_left.BeamLoc_eta = all_eta;

    %Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_left = awi.model.Spar;
    FrontSpar_left.XLoc = [0.15, 0.15];
    FrontSpar_left.Eta  = [0,  1];
    RearSpar_left = awi.model.Spar;
    RearSpar_left.XLoc = [0.65, 0.65];
    RearSpar_left.Eta  = [0, 1];

    Wingbox_left.add([FrontSpar_left, RearSpar_left]);

    %Define internal layout
    Wingbox_left.RibPitch      = 0.65;
    Wingbox_left.StringerPitch = 0.15;

    Wingbox_left.Material_eta = [0, 1];
    Wingbox_left.Material     = [Mat_wing, Mat_wing];
    
    build(Wingbox_left)
    
    % Aeropanel definition
    
    Wingbox_left.AeroPanelLength=[];
    Wingbox_left.NumAeroPanel=Param.Wing.AeroPanel_Number;
    Wingbox_left.AeroPanelAR=Param.Wing.AeroPanel_AR;
    
    
    

