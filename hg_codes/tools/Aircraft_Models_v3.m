function [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, varargout]=Aircraft_Models_v3(Param)

    %% Right root connector

    Connector_right = awi.model.LiftingSurface;
    Connector_right.Name = 'Connector_Right';
    Connector_right.Origin=[Param.Layout.Wing_Position,0,0]; %15

    %Use the Leading/Trailing edge sweep to define the planform
    Connector_right.ActiveSet = 'sSet';

    %Tail wing dimensions
    Connector_right.SpanVector  = 'Y';
    Connector_right.Span        = 2;
    Connector_right.LESweep     = [0, 0];
    Connector_right.LESweep_eta = [0, 1];
    Connector_right.TESweep     = [0,  0];
    Connector_right.TESweep_eta = [0,  1];
    Connector_right.RootChord   = Param.Wing.Root_Chord;

    %Make sure the beam is at the midchord
    all_eta           = Connector_right.Eta_;
    Connector_right.BeamLoc     = repmat(Param.Wing.BeamLoc, size(all_eta));
    Connector_right.BeamLoc_eta = all_eta;


    % Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_root_right = awi.model.Spar;
    FrontSpar_root_right.XLoc = [0.15, 0.15];
    FrontSpar_root_right.Eta  = [0   , 1];
    RearSpar_root_right = awi.model.Spar;
    RearSpar_root_right.XLoc = [0.65, 0.65];
    RearSpar_root_right.Eta  = [0   , 1];
    Connector_right.add([FrontSpar_root_right, RearSpar_root_right]);

    %Define internal layout
    Connector_right.RibPitch      = 0.65;
    Connector_right.StringerPitch = 0.15;
    
    %Make the connector material
    E0  = 15.4e10; %[N/m^2], typical YM of IM7 composite
    nu0 = 0.333;
    rho0=1550; 
    Mat_conn = awi.model.Material;
    Mat_conn.E  = E0;
    Mat_conn.Nu = nu0;
    Mat_conn.G  = E0 / (2 * (1 + nu0));
    Mat_conn.Rho=rho0;

    % material properties
    Connector_right.Material_eta = [0, 1];
    Connector_right.Material     = [Mat_conn, Mat_conn];

    % Define box beam corss section
    Connector_box_right=awi.model.BoxBeam;
    Connector_box_right.BoxType='SymmetricBox';
    Connector_box_right.Height=1;
    Connector_box_right.Width=3;
    Connector_box_right.CoverThickness=0.08;
    Connector_box_right.SparThickness=0.08;
    getGeometricProps(Connector_box_right)
    
    Connector_right.BoxBeam = Connector_box_right;
    Connector_right.A   = Connector_box_right.Abb;
    Connector_right.I11 = Connector_box_right.Ixx;
    Connector_right.I22 = Connector_box_right.Izz;
    Connector_right.J   = Connector_box_right.Jbb;
%     Connector_right.NSM = Connector_box_right.NSM;
%     Connector_right.NSI = Connector_box_right.NSI;

    for i=1:1:3
        handle_connectorR=strcat('PM_tail_R','i');
        handle_connectorR=awi.model.PointMass;
        handle_connectorR.SOffset=-0.1+i*0.2;
        handle_connectorR.Mass=1;
        handle_connectorR.MassGroup='Group3';
        Connector_right.add(handle_connectorR);

    end
    
    % Aeropanel definition
    Connector_right.AeroPanelLength=[];
    Connector_right.NumAeroPanel=Param.Wing.AeroPanel_Number;
    
    % Jig shape
    
    Connector_right.Twist=Param.Connector.Jig_Twist;
    
    Connector_right.Twist_eta=Param.Connector.Jig_Eta;

    build(Connector_right);
    
    
    %% Left root connector
    
    Connector_left = awi.model.LiftingSurface;
    Connector_left.Name = 'Connector_Right';
    Connector_left.Origin=[Param.Layout.Wing_Position,0,0]; %15
    
    %Use the Leading/Trailing edge sweep to define the planform
    Connector_left.ActiveSet = 'sSet';

    %Tail wing dimensions
    Connector_left.SpanVector  = 'Y';
    Connector_left.Span        = -2;
    Connector_left.LESweep     = [0, 0];
    Connector_left.LESweep_eta = [0, 1];
    Connector_left.TESweep     = [0,  0];
    Connector_left.TESweep_eta = [0,  1];
    Connector_left.RootChord   = Param.Wing.Root_Chord;

    %Make sure the beam is at the midchord
    all_eta           = Connector_left.Eta_;
    Connector_left.BeamLoc     = repmat(Param.Wing.BeamLoc, size(all_eta));
    Connector_left.BeamLoc_eta = all_eta;
%     Connector_right.XOffset=35;
%     Tailwing_right.YOffset=1;

    % Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_root_left = awi.model.Spar;
    FrontSpar_root_left.XLoc = [0.15,  0.15];
    FrontSpar_root_left.Eta  = [0   ,  1];
    RearSpar_root_left = awi.model.Spar;
    RearSpar_root_left.XLoc = [0.65, 0.65];
    RearSpar_root_left.Eta  = [0   , 1];
    Connector_left.add([FrontSpar_root_left, RearSpar_root_left]);

    %Define internal layout
    Connector_left.RibPitch      = 0.65;
    Connector_left.StringerPitch = 0.15;
    
    % material properties
    Connector_left.Material_eta = [0, 1];
    Connector_left.Material     = [Mat_conn, Mat_conn];

    % Define box beam corss section
    Connector_box_left=awi.model.BoxBeam;
    Connector_box_left.BoxType='SymmetricBox';
    Connector_box_left.Height=1;
    Connector_box_left.Width=3;
    Connector_box_left.CoverThickness=0.08;
    Connector_box_left.SparThickness=0.08;
    getGeometricProps(Connector_box_left)
    
    Connector_left.BoxBeam = Connector_box_left;
    Connector_left.A   = Connector_box_left.Abb;
    Connector_left.I11 = Connector_box_left.Ixx;
    Connector_left.I22 = Connector_box_left.Izz;
    Connector_left.J   = Connector_box_left.Jbb;
%     Connector_left.NSM = Connector_box_left.NSM;
%     Connector_left.NSI = Connector_box_left.NSI;

    for i=1:1:3
        handle_connectorL=strcat('PM_tail_R','i');
        handle_connectorL=awi.model.PointMass;
        handle_connectorL.SOffset=-0.1+i*0.2;
        handle_connectorL.Mass=1;
        handle_connectorL.MassGroup='Group3';
        Connector_left.add(handle_connectorL);

    end
    
    % Aeropanel definition
    Connector_left.AeroPanelLength=[];
    Connector_left.NumAeroPanel=Param.Wing.AeroPanel_Number;
    
    % Jig shape
    
    Connector_left.Twist=-Param.Connector.Jig_Twist;
    
    Connector_left.Twist_eta=Param.Connector.Jig_Eta;

    build(Connector_left);

%% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[Param.Layout.Wing_Position,2,0]; %15
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';
    
    % Num of element 
    Wingbox_right.NumBeamElem = 23;

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = Param.Wing.Semi_Span;   
    Wingbox_right.LESweep     = [Param.Wing.LE_Sweep, Param.Wing.LE_Sweep];
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [Param.Wing.TE_Sweep1, Param.Wing.TE_Sweep2, Param.Wing.TE_Sweep2];
    Wingbox_right.TESweep_eta = [0, Param.Wing.Kink, 1];
    Wingbox_right.RootChord   = Param.Wing.Root_Chord;
    
    
    %Dihedral 
    Wingbox_right.Dihedral=[Param.Wing.Dihedral,Param.Wing.Dihedral];
    Wingbox_right.Dihedral_eta=[0,1];
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_right.Eta_;
    Wingbox_right.BeamLoc     = repmat(Param.Wing.BeamLoc, size(all_eta));

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
    Wingbox_right.NumAeroPanel=Param.Wing.AeroPanel_Number;
    Wingbox_right.AeroPanelAR=Param.Wing.AeroPanel_AR;
    
    %     Wingbox_right.AeroPanelLength=0.4;
    
    Kink=Param.Wing.Kink;
    
    
    %% FWT
    
    if isfield(Param,'FWT')
        
        %% initialise properties
        
        Wingbox_right.A   =  [1,1];
        Wingbox_right.A_eta=[0,1];
        
        Wingbox_right.I11 = [1,1];
        Wingbox_right.I11_eta=[0,1];
        
        Wingbox_right.I22 = [1,1];
        Wingbox_right.I22_eta = [0,1];
        
        Wingbox_right.J   = [1,1];
        Wingbox_right.J_eta= [0,1];
        
        Wingbox_right.Twist = deg2rad([0, 0]);
        Wingbox_right.Twist_eta = [0, 1];
        
        %% FWT right definition
        
        FWT_R = insertWingFold(Wingbox_right, 'FlareAngle', Param.FWT.Flare_angle, 'FoldAngle', Param.FWT.Fold_angle,'EtaFold',Param.FWT.Fold_eta);
        FWT_R.HingeStiffness = [1e14 1e14 1e14 1e14 Param.FWT.Hinge_Stiffness 1e14];
        
        FWT_R.AeroPanelLength=[];
        FWT_R.NumAeroPanel=Param.Wing.AeroPanel_Number;
        FWT_R.AeroPanelAR=Param.Wing.AeroPanel_AR;
        
        FWT_R.NumBeamElem = 10;
        
        FWT_eta_= 0:0.1:1;
        FWT_Bheight=interp1([0,1],[Param.FWT.Root_Height,Param.FWT.Tip_Height],FWT_eta_);
        FWT_Bwidth=interp1([0,1],0.5*[Param.FWT.Root_Chord,Param.FWT.Tip_Chord],FWT_eta_);
        
  
        % computing beam properties ------------------------------------
        
        [FWT_A_val,FWT_Iyy_val, FWT_Izz_val, FWT_J_val]=Beam_Proeprties_v1(FWT_Bwidth, FWT_Bheight, Param,'FWT',true,'Wing',false);
        
        % --------------------------------------------------------------
        
        % FWT structural properties
        
        FWT_R.A   =  FWT_A_val.cross_section;
        FWT_R.A_eta=FWT_eta_;
        
        FWT_R.I11 = FWT_Izz_val;
        FWT_R.I11_eta=FWT_eta_;
        
        FWT_R.I22 = FWT_Iyy_val;
        FWT_R.I22_eta = FWT_eta_;
        
        FWT_R.J   = FWT_J_val;
        FWT_R.J_eta= FWT_eta_;
        
        
        % FWT Jig Twist
        FWT_R.Twist = Param.FWT.Jig_Twist;
        FWT_R.Twist_eta = Param.FWT.Jig_Eta;
        
        
        %Make the material for FWT
        E_fwt  = 70e9; %[N/m^2], typical YM of aluminium
        nu_fwt = 0.333;
        rho_fwt=2810;
        Mat_fwt = awi.model.Material;
        Mat_fwt.E  = E_fwt;
        Mat_fwt.Nu = nu_fwt;
        Mat_fwt.G  = E_fwt / (2 * (1 + nu_fwt));
        Mat_fwt.Rho=rho_fwt;
        FWT_R.Material_eta = [0, 1];
        FWT_R.Material     = [Mat_fwt, Mat_fwt];
        
        
        % add point masses- distributed
        for i=1:1:2
            handle=strcat('PM_right','i');
            handle=awi.model.PointMass;
            handle.SOffset=0.2+i*0.2;
            handle.Mass=1;
            handle.Inertia11 = 0.1;
            handle.Inertia22 =  0.1;
            handle.Inertia33 =  0.1;
            handle.Inertia12 =  0.1;
            handle.Inertia23 =  0.1;
            handle.Inertia13 =  0.1;
            handle.MassGroup='Group1';
            FWT_R.add(handle);
            
        end
        
        % hinge weight
        hinge_weight=awi.model.PointMass;
        hinge_weight.SOffset=0;
        hinge_weight.Mass=Param.Masses.Hinge;
        FWT_R.add(hinge_weight);
        
        varargout{1}=FWT_R;
        
        Kink=Param.Wing.Kink*(1/Param.FWT.Fold_eta);
        
        
    end
   
     
  %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=Wingbox_right.NumBeamElem+2;
    

    % -------------------------------------------------------------

    % set width and height array 
    YData=Wingbox_right.YData;
    SparWidth=Wingbox_right.Chord*0.5;

    RootH=Wingbox_right.Chord(1)*Param.Wing.ThicknessToChord_Root; % root thickness/chord = 0.15
    MidH=Wingbox_right.Chord(2)*Param.Wing.ThicknessToChord_kink;  % middle thickness/chord = 0.12
    TipH=Wingbox_right.Chord(end)*Param.Wing.ThicknessToChord_tip;% tip thickness/chord = 0.11


    % set up eta values
    elnum=Wingbox_right.NumBeamElem + 1; % total number of beam elements along the wing
    
    if isfield(Param,'FWT')
        Num_seg1=round(elnum*Kink); % number of elements in the inboard section
        
    else
        
        Num_seg1=ceil(elnum*Kink); % number of elements in the inboard section
        
    end
          
    Num_seg2=elnum - Num_seg1; % number of elements in the outboard section
    
    Num_sec1=Num_seg1+1;    % number of sections in the inboard section
    Num_sec2=Num_seg2+1;    % number of sections in the outboard section
    
    eta1_=linspace(0,Kink, Num_sec1);
    eta2_=linspace(Kink,1,Num_sec2);
    etaS=[eta1_(1:end-1),eta2_(1:end)];

    RData=Wingbox_right.RData;
    eta_R=RData/RData(end);
    eta_Y=YData/YData(end);
    etaRS=interp1(eta_Y,eta_R,etaS);

    Bwidth=interp1(RData/RData(end),SparWidth,etaRS);
    Bheight=interp1(RData/RData(end),0.79*[RootH,MidH,TipH],etaRS);
    
    % computing beam properties -----------------------------------------
    
    [A_val,Iyy_val, Izz_val, J_val]=Beam_Proeprties_v1(Bwidth, Bheight, Param,'FWT',false,'Wing',true);
    
    % -------------------------------------------------------------------

    eta_=etaRS;
    Y_eta=etaRS*RData(end);
    
    Wingbox_right.A   =  A_val.cross_section;
    Wingbox_right.A_eta=eta_;

    Wingbox_right.I11 = Izz_val;
    Wingbox_right.I11_eta=eta_;

    Wingbox_right.I22 = Iyy_val;
    Wingbox_right.I22_eta = eta_;

    Wingbox_right.J   = J_val;
    Wingbox_right.J_eta= eta_;


    % Jig shape
    
    Wingbox_right.Twist = Param.Wing.Jig_Twist;
    Wingbox_right.Twist_eta = Param.Wing.Jig_Eta;
    
    build(Wingbox_right)
    
    if isfield(Param,'FWT')
        
        Box_dimensions.Inboard.Height = [Bheight,FWT_Bheight(2:end)];
        Box_dimensions.Inboard.Width = [Bwidth,FWT_Bwidth(2:end)];

        
        Box_CrossSec.Izz=[Izz_val,FWT_Izz_val(2:end)];
        Box_CrossSec.Iyy=[Iyy_val,FWT_Iyy_val(2:end)];
        Box_CrossSec.J=[J_val,FWT_J_val(2:end)];
        Box_CrossSec.A=[A_val.cross_section,FWT_A_val.cross_section(2:end)];
        
    else
        
        Box_dimensions.Inboard.Height = Bheight;
        Box_dimensions.Inboard.Width = Bwidth;

        
        Box_CrossSec.Izz=Izz_val;
        Box_CrossSec.Iyy=Iyy_val;
        Box_CrossSec.J=J_val;
        Box_CrossSec.A=A_val.cross_section;
        
    end
    
    %% Mass definition
    

    wingmass_eta=etaRS;
    Mwidth=interp1(RData/RData(end),SparWidth,wingmass_eta);
    
    mass_set=(Param.Masses.Secondary_Mass+0.5*Param.Masses.Fuel_Mass)*(Mwidth)/sum(Mwidth);
      
%     m=total_mass/19;
    
    for i=1:1:25
        handle=strcat('PM_right','i');
        handle=awi.model.PointMass;
        handle.SOffset=wingmass_eta(i);
%         handle.SOffset=0+i*0.2;
        handle.Mass=mass_set(i);
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle.MassGroup='Group1';
        Wingbox_right.add(handle);

    end
    
    
 
    %% attachments - engine
    
    Engine=awi.model.BluffBody;
    Engine.Name='Engine';
    
    % cylinder body
    Engine.Radius=[1.4, 1.4, 1];
    Engine.Eta =  [0, 0.6, 1];
    Engine.Length = 3.5;    

    % Engine location - user defined
    Y_Engine=Param.Layout.Engine_Position;
    
    X_Engine=interp1(Wingbox_right.YData,Wingbox_right.XData,Y_Engine);
    Z_Engine=interp1(Wingbox_right.YData,Wingbox_right.ZData,Y_Engine);
    
    Engine.Origin = [X_Engine-Engine.Length+Param.Layout.Wing_Position, Y_Engine + 2, Z_Engine];

    
    %Make engine material
    E1  = 70e9; %[N/m^2],set as a rigid body
    nu = 0.333;
    Engine_Mat = awi.model.Material;
    Engine_Mat.E  = E1;
    Engine_Mat.Nu = nu;
    Engine_Mat.G  = E1 / (2 * (1 + nu));
    Engine_Mat.Rho=1; % using lumped mass instead
    
    
    % use the strong material
    Engine.Material_eta = [0, 1];
    Engine.Material     = [Engine_Mat, Engine_Mat];
    
    % Engine stiffness
    Engine_radius=1;
    Engine_thickness=0.015;
    Engine_A=2*pi*Engine_radius*Engine_thickness;
    
    Engine_I11=pi*Engine_radius^3*Engine_thickness;
    Engine_I22=pi*Engine_radius^3*Engine_thickness;
    Engine_J=2*pi*Engine_radius^3*Engine_thickness;

    
    Engine.A   = Engine_A;
    Engine.I11 = Engine_I11;
    Engine.I22 = Engine_I22;
    Engine.J   = Engine_J;

    
    %Aeropanel althoufh it is useless now
    Engine.AeroPanelLength=0.5;
    
    % add engine mass
    engine_mass=awi.model.PointMass;   
    engine_mass.SOffset=0.1;
    engine_mass.Mass=Param.Masses.Engine;
    Engine.add(engine_mass);
    
    % add pylon
    pylon_mass=awi.model.PointMass;   
    pylon_mass.SOffset=0.9;
    pylon_mass.Mass=Param.Masses.Pylon;
    Engine.add(pylon_mass);

    build(Engine)
      
    Wingbox_right.add(Engine)

    
%     %Control surfaces - flaps
%     flap_R=awi.model.ControlSurface;
%     flap_R.Eta=[0, 0.24];
%     flap_R.xLE=[0.8,0.8];
%     flap_R.xTE=[1,1];
%     flap_R.Max_def=0.1;
%     flap_R.Max_rate=0.1;
%     flap_R.HingeLine='LE';
%     flap_R.Label='FlapR';
%     flap_R.FaceColor='m';
%     
% %     flap_R.NumAeroPanel=10;
%     flap_R.AeroPanelLength=0.4;
%     
%     build(flap_R)
%     Wingbox_right.add(flap_R);
%     
%     Wingbox_right.ModelControlSurf = 1;
    
    
    build(Wingbox_right);
    
 
    
    %% Wingbox 2 - left and control surf.

    Wingbox_left = awi.model.LiftingSurface;
    Wingbox_left.Name = 'A320Wing_left';
    Wingbox_left.Origin=[Param.Layout.Wing_Position,-2,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_left.ActiveSet = 'sSet';
    
    % Num of element
    Wingbox_left.NumBeamElem = 23;

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
    
    % AeroPanelLength
    %     NumAeroPanel
    %     Wingbox_left.NumAeroPanel=20;
    %     Wingbox_left.AeroPanelLength=0.4;
    Wingbox_left.AeroPanelLength=[];
    Wingbox_left.NumAeroPanel=Param.Wing.AeroPanel_Number;
    Wingbox_left.AeroPanelAR=Param.Wing.AeroPanel_AR;
    
    
     %% FWT
    
     if isfield(Param,'FWT')
         
         Wingbox_left.A   =  1;
         Wingbox_left.A_eta=[0,1];
         
         Wingbox_left.I11 = 1;
         Wingbox_left.I11_eta=[0,1];
         
         Wingbox_left.I22 = 1;
         Wingbox_left.I22_eta = [0,1];
         
         Wingbox_left.J   = 1;
         Wingbox_left.J_eta= [0,1];
         
         %% FWT left definition
         
         FWT_L = insertWingFold(Wingbox_left, 'FlareAngle', -Param.FWT.Flare_angle, 'FoldAngle', Param.FWT.Fold_angle,'EtaFold',Param.FWT.Fold_eta);
         FWT_L.HingeStiffness = [1e14 1e14 1e14 1e14 Param.FWT.Hinge_Stiffness 1e14];
         
         FWT_L.AeroPanelLength=[];
         FWT_L.NumAeroPanel=Param.Wing.AeroPanel_Number;
         FWT_L.AeroPanelAR=Param.Wing.AeroPanel_AR;
         
         FWT_L.NumBeamElem = 10;
         
         
         % FWT structural properties
         FWT_eta_= 0:0.1:1;
         
         FWT_L.A   =  FWT_A_val.cross_section;
         FWT_L.A_eta=FWT_eta_;
         
         FWT_L.I11 = FWT_Izz_val;
         FWT_L.I11_eta=FWT_eta_;
         
         FWT_L.I22 = FWT_Iyy_val;
         FWT_L.I22_eta = FWT_eta_;
         
         FWT_L.J   = FWT_J_val;
         FWT_L.J_eta= FWT_eta_;
         
         FWT_L.Material_eta = [0, 1];
         FWT_L.Material     = [Mat_fwt, Mat_fwt];
         
         % Jig shape
         FWT_L.Twist = -Param.FWT.Jig_Twist;
         FWT_L.Twist_eta = Param.FWT.Jig_Eta;
         
         % add point masses
         for i=1:1:2
             handle=strcat('PM_right','i');
             handle=awi.model.PointMass;
             handle.SOffset=0.2+i*0.2;
             handle.Mass=1;
             handle.Inertia11 = 0.1;
             handle.Inertia22 =  0.1;
             handle.Inertia33 =  0.1;
             handle.Inertia12 =  0.1;
             handle.Inertia23 =  0.1;
             handle.Inertia13 =  0.1;
             handle.MassGroup='Group1';
             FWT_L.add(handle);
             
         end
         
         % hinge mass 
         hinge_weight_L=awi.model.PointMass;
         hinge_weight_L.SOffset=0;
         hinge_weight_L.Mass=Param.Masses.Hinge;
         FWT_L.add(hinge_weight_L);
         
     end
    
       
    %% Create discretised boxbeam with varied cross section prperties along the span 

    eta_=etaRS;
    Wingbox_left.A   =  A_val.cross_section;
    Wingbox_left.A_eta=eta_;

    Wingbox_left.I11 = Izz_val;
    Wingbox_left.I11_eta=eta_;

    Wingbox_left.I22 = Iyy_val;
    Wingbox_left.I22_eta = eta_;

    Wingbox_left.J   = J_val;
    Wingbox_left.J_eta= eta_;
    
    % Jig shape
    
    Wingbox_left.Twist = -Param.Wing.Jig_Twist;
    Wingbox_left.Twist_eta = Param.Wing.Jig_Eta;

    build(Wingbox_left)
   

    %% Mass definition
    
    
    for i=1:1:25
        handle=strcat('PM_left','i');
        handle=awi.model.PointMass;
        handle.SOffset=wingmass_eta(i);
%         handle.SOffset=0+i*0.2;
        handle.Mass=mass_set(i);
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle.MassGroup='Group1';
        Wingbox_left.add(handle);

    end
    

    %% attachments 2  - engine_left
    
    Engine2=awi.model.BluffBody;
    Engine2.Name='Engine_left';
    
    % cylinder body
    Engine2.Radius=[1.4, 1.4, 1];
    Engine2.Eta =  [0, 0.6, 1];
    Engine2.Length = 3.5;
    
%     %Engine location 
%     Y_eng=Wingbox_left.PanelCoords.LE.Y(2);
%     X_eng=(Wingbox_left.PanelCoords.LE.X(2)+Wingbox_left.PanelCoords.TE.X(2))/2;
%     Z_eng=(Wingbox_left.PanelCoords.LE.Z(2)+Wingbox_left.PanelCoords.TE.Z(2))/2;
%     
%     Engine2.Origin = [X_eng-3.5, Y_eng, Z_eng];

    Engine2.Origin = [X_Engine-Engine.Length+Param.Layout.Wing_Position, -(Y_Engine + 2), Z_Engine];
   
%     Engine.XOffset=16.471008-3.5;
%     Engine.YOffset=5.8170588;
%     Engine.ZOffset=-2;
    
%     %Make the material
%     E1  = 76e9; %[N/m^2],set as a rigid body
%     nu = 0.333;
%     Mat1 = awi.model.Material;
%     Mat1.E  = E1;
%     Mat1.Nu = nu;
%     Mat1.G  = E1 / (2 * (1 + nu));
%     Mat1.Rho=2800;
    
    
    % use the strong material
    Engine2.Material_eta = [0, 1];
    Engine2.Material     = [Engine_Mat, Engine_Mat];
    Engine2.A   = Engine_A;
    Engine2.I11 = Engine_I11;
    Engine2.I22 = Engine_I22;
    Engine2.J   = Engine_J;

%     Engine2.A   = 0.04432;
%     Engine2.I11 = 0.002;
%     Engine2.I22 = 0.002;
%     Engine2.J   = 0.001636;
    
    %Aeropanel althoufh it is useless now
    Engine2.AeroPanelLength=0.5;
    
    % add engine mass
    engine2_mass=awi.model.PointMass;   
    engine2_mass.SOffset=0.1;
    engine2_mass.Mass=Param.Masses.Engine;
    Engine2.add(engine2_mass);
    
    % add pylon
    pylon2_mass=awi.model.PointMass;   
    pylon2_mass.SOffset=0.9;
    pylon2_mass.Mass=Param.Masses.Pylon;
    Engine2.add(pylon2_mass);

    build(Engine2)
      
    Wingbox_left.add(Engine2)
    
%     %Control surfaces - flaps
%     flap_L=awi.model.ControlSurface;
%     flap_L.Eta=[0, 0.24];
%     flap_L.xLE=[0.8,0.8];
%     flap_L.xTE=[1,1];
%     flap_L.Max_def=0.1;
%     flap_L.Max_rate=0.1;
%     flap_L.HingeLine='LE';
%     flap_L.Label='FlapL';
%     flap_L.FaceColor='m';
%     
%     flap_L.AeroPanelLength=0.4;
%     
%     build(flap_R)
%     Wingbox_left.add(flap_L);
%     
%     Wingbox_left.ModelControlSurf = 1;
     
    build(Wingbox_left);
    
    
    
       %% Create a BluffBody

    Body=awi.model.BluffBody;
    Body.Name='Fuselage';
    % cylinder body
    % Body.Radius=[2,2];
    % Body.Eta=[0,1];

    % real body
    Body.Eta = [0;0.005268;0.010536;0.015805;0.021073;...
        0.026342;0.03161;0.036879;0.042147;0.047415;0.052684;...
        0.057952;0.063221;0.0684890;0.073758;0.079026;0.084294;0.089563;0.094831;0.1001;0.411022;...
        0.721944;0.736578;0.751213;0.765847;0.780482;0.795117;0.809751;0.824386;0.83902;0.853655;...
        0.868289;0.882924;0.897558;0.912193;0.926827;0.941462;0.956096;0.970731;0.985365;1]';

    Body.Radius = [0.01;0.3844030;0.565081;0.707928;0.830682;0.940375;...
        1.04067;1.13377;1.22112;1.30374;1.38237;1.45758;1.52981;1.59941;1.66667;...
        1.73182;1.79508;1.8566;1.91653;1.975;2.11455;2.11455;2.11455;2.11455;2.11455;...
        2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;1.9;...
        1.75;1.6;1.4;1.2;1.0;0.01]';
    
    Body.Origin = [0, 0, 0];
    Body.Length=Param.Layout.Fuselage_Length;
    % Body.XOffset=-15;

    %Make the material
    E1  = 70e9; %[N/m^2],set as a rigid body
    nu = 0.333;
    Body_Mat = awi.model.Material;
    Body_Mat.E  = E1;
    Body_Mat.Nu = nu;
    Body_Mat.G  = E1 / (2 * (1 + nu));
    Body_Mat.Nu = nu;
    Body_Mat.Rho = 2800; 

    
    % use the strong material

    Body.Material_eta = [0, 1];
    Body.Material     = [Body_Mat, Body_Mat];

    %define  panel size
    % Body.NumAeroPanel=5;
    Body.AeroPanelLength=0.5;

    % define cross sectional definition
%     Bodybox=awi.model.BoxBeam;
%     Bodybox.BoxType='SymmetricBox';
%     Bodybox.Height=2;
%     Bodybox.Width=2;
%     Bodybox.CoverThickness=0.01;
%     Bodybox.SparThickness=0.01;
%     getGeometricProps(Bodybox)
%     
%     
%     Body.BoxBeam = Bodybox;

    Body_radius=2; 
    Body_thickness=0.004;
    CS_A=2*pi*Body_radius*Body_thickness;
    
    CS_I11=pi*Body_radius^3*Body_thickness;
    CS_I22=pi*Body_radius^3*Body_thickness;
    CS_J=2*pi*Body_radius^3*Body_thickness;
%     CS_Inertia = Mat1.Rho*CS_A*Body_radius^2;
    
    Body.A   = CS_A;
    Body.I11 = CS_I11;
    Body.I22 = CS_I22;
    Body.J   = CS_J;
    %     Body.NSM = Bodybox.NSM;
%     Body.NSI = CS_Inertia;
    
    % add point masses- fuselage
    %     M_fs=MWE-Wing_mass;
    %     m_fs=M_fs/11;
    %
    %     %payload
    %     M_p=Payload_max;
    %     m_p= M_p/11;
    
    % mass including structure mass + payload
    Mass_val=(Param.Masses.Fuselage_Structure + Param.Masses.Payload)/11;

    
    for i=1:1:11
        handle=strcat('PM_body','i');
        handle=awi.model.PointMass;
        handle.SOffset=-0.1+i*0.1;
        handle.Mass=Mass_val;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle.MassGroup='Group2';
        Body.add(handle);

    end

    build(Body)


    %% Generate tailwing Right and control surf.

    Tailwing_right = awi.model.LiftingSurface;
    Tailwing_right.Name = 'Tail_Wing_Right';

    %Use the Leading/Trailing edge sweep to define the planform
    Tailwing_right.ActiveSet = 'sSet';

    %Tail wing dimensions
    Tailwing_right.SpanVector  = 'Y';
    Tailwing_right.Span        = 12.45/2;
    Tailwing_right.LESweep     = [32, 32];
    Tailwing_right.LESweep_eta = [0, 1];
    Tailwing_right.TESweep     = [15,  15];
    Tailwing_right.TESweep_eta = [0,  1];
    Tailwing_right.RootChord   = 3.31;
    
    %Dihedral
    Tailwing_right.Dihedral=[Param.Wing.Dihedral,Param.Wing.Dihedral];
    Tailwing_right.Dihedral_eta=[0,1];

    %Make sure the beam is at the midchord
    all_eta           = Tailwing_right.Eta_;
    Tailwing_right.BeamLoc     = repmat(0.4, size(all_eta));
    Tailwing_right.BeamLoc_eta = all_eta;
    Tailwing_right.XOffset = Param.Layout.Horizontal_Tail_Position;
%     Tailwing_right.YOffset=1;

    % Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_tail_right = awi.model.Spar;
    FrontSpar_tail_right.XLoc = [0.15, 0.15];
    FrontSpar_tail_right.Eta  = [0   , 1];
    RearSpar_tail_right = awi.model.Spar;
    RearSpar_tail_right.XLoc = [0.65, 0.65];
    RearSpar_tail_right.Eta  = [0   , 1];
    Tailwing_right.add([FrontSpar_tail_right, RearSpar_tail_right]);

    %Define internal layout
    Tailwing_right.RibPitch      = 0.65;
    Tailwing_right.StringerPitch = 0.15;
    
    %Make the material
    Et  = 70e9; %[N/m^2],set as a rigid body
    nut = 0.333;
    Tail_Mat = awi.model.Material;
    Tail_Mat.E  = Et;
    Tail_Mat.Nu = nut;
    Tail_Mat.G  = Et / (2 * (1 + nut));
    Tail_Mat.Nu = nut;
    Tail_Mat.Rho = 2800;

    % material properties
    Tailwing_right.Material_eta = [0, 1];
    Tailwing_right.Material     = [Tail_Mat, Tail_Mat];

    % Define box beam corss section
    tailbox_right=awi.model.BoxBeam;
    tailbox_right.BoxType='SymmetricBox';
    tailbox_right.Height=0.5;
    tailbox_right.Width=1;
    tailbox_right.CoverThickness=0.006;
    tailbox_right.SparThickness=0.006;
    getGeometricProps(tailbox_right)
    Tailwing_right.BoxBeam = tailbox_right;
    Tailwing_right.A   = tailbox_right.Abb;
    Tailwing_right.I11 = tailbox_right.Ixx;
    Tailwing_right.I22 = tailbox_right.Izz;
    Tailwing_right.J   = tailbox_right.Jbb;
%     Tailwing_right.NSM = tailbox_right.NSM;
%     Tailwing_right.NSI = tailbox_right.NSI;

%     Horiaontal_tail_mass_val=0.5*Horizontal_tail/5;
%     
%     for i=1:1:5
%         handle_tailR=strcat('PM_tail_R','i');
%         handle_tailR=awi.model.PointMass;
%         handle_tailR.SOffset=-0.1+i*0.2;
%         handle_tailR.Mass=Horiaontal_tail_mass_val;
% %         handle.Inertia11 =  0;
% %         handle.Inertia22 =  0;
% %         handle.Inertia33 =  0;
% %         handle.Inertia12 =  0;
% %         handle.Inertia23 =  0;
% %         handle.Inertia13 =  0;
%         handle.MassGroup='Group3';
%         Tailwing_right.add(handle_tailR);
% 
%     end
    

    % Aeropanel definition
    Tailwing_right.AeroPanelLength=0.5;

    %Control surfaces - elevators
    myelevator_right=awi.model.ControlSurface;
    myelevator_right.Eta=[0, 1];
    myelevator_right.xLE=[0.6,0.6];
    myelevator_right.xTE=[1,1];
    myelevator_right.Max_def=0.1;
    myelevator_right.Max_rate=0.1;
    myelevator_right.HingeLine='LE';
    myelevator_right.Label='elevatR';
    myelevator_right.FaceColor='m';
    myelevator_right.AeroPanelLength=0.5;
    build(myelevator_right)
    Tailwing_right.add(myelevator_right);

    Tailwing_right.ModelControlSurf = 1;


    build(Tailwing_right);
    
    
      %% Generate tailwing Left and control surf.

    Tailwing_left = awi.model.LiftingSurface;
    Tailwing_left.Name = 'Tail_Wing_Left';

    %Use the Leading/Trailing edge sweep to define the planform
    Tailwing_left.ActiveSet = 'sSet';

    %Tail wing dimensions
    Tailwing_left.SpanVector  = 'Y';
    Tailwing_left.Span        = -12.45/2;
    Tailwing_left.LESweep     = [-32, -32];
    Tailwing_left.LESweep_eta = [0, 1];
    Tailwing_left.TESweep     = [-15,  -15];
    Tailwing_left.TESweep_eta = [0,  1];
    Tailwing_left.RootChord   = 3.31;
    
    %Dihedral
    Tailwing_left.Dihedral=[Param.Wing.Dihedral,Param.Wing.Dihedral];
    Tailwing_left.Dihedral_eta=[0,1];

    %Make sure the beam is at the midchord
    all_eta           = Tailwing_left.Eta_;
    Tailwing_left.BeamLoc     = repmat(0.4, size(all_eta));
    Tailwing_left.BeamLoc_eta = all_eta;
    Tailwing_left.XOffset=Param.Layout.Horizontal_Tail_Position;
%     Tailwing_left.YOffset=-1;

    % Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_tail_left = awi.model.Spar;
    FrontSpar_tail_left.XLoc = [0.15, 0.15];
    FrontSpar_tail_left.Eta  = [0   , 1];
    RearSpar_tail_left = awi.model.Spar;
    RearSpar_tail_left.XLoc = [0.65, 0.65];
    RearSpar_tail_left.Eta  = [0   , 1];
    Tailwing_left.add([FrontSpar_tail_left, RearSpar_tail_left]);

    %Define internal layout
    Tailwing_left.RibPitch      = 0.65;
    Tailwing_left.StringerPitch = 0.15;

    % material properties
    Tailwing_left.Material_eta = [0, 1];
    Tailwing_left.Material     = [Tail_Mat, Tail_Mat];

    % Define box beam corss section
    tailbox_left=awi.model.BoxBeam;
    tailbox_left.BoxType='SymmetricBox';
    tailbox_left.Height=0.5;
    tailbox_left.Width=1;
    tailbox_left.CoverThickness=0.006;
    tailbox_left.SparThickness=0.006;
    getGeometricProps(tailbox_left)
    Tailwing_left.BoxBeam = tailbox_left;
    Tailwing_left.A   = tailbox_left.Abb;
    Tailwing_left.I11 = tailbox_left.Ixx;
    Tailwing_left.I22 = tailbox_left.Izz;
    Tailwing_left.J   = tailbox_left.Jbb;
%     Tailwing_right.NSM = tailbox_right.NSM;
%     Tailwing_right.NSI = tailbox_right.NSI;

%     for i=1:1:5
%         handle_tailL=strcat('PM_tail_L','i');
%         handle_tailL=awi.model.PointMass;
%         handle_tailL.SOffset=-0.1+i*0.2;
%         handle_tailL.Mass=Horiaontal_tail_mass_val;
% %         handle.Inertia11 =  0;
% %         handle.Inertia22 =  0;
% %         handle.Inertia33 =  0;
% %         handle.Inertia12 =  0;
% %         handle.Inertia23 =  0;
% %         handle.Inertia13 =  0;
%         handle_tailL.MassGroup='Group4';
%         Tailwing_left.add(handle_tailL);
% 
%     end

    % Aeropanel definition
    Tailwing_left.AeroPanelLength=0.5;

    %Control surfaces - elevators
    myelevator_left=awi.model.ControlSurface;
    myelevator_left.Eta=[0, 1];
    myelevator_left.xLE=[0.6,0.6];
    myelevator_left.xTE=[1,1];
    myelevator_left.Max_def=0.1;
    myelevator_left.Max_rate=0.1;
    myelevator_left.HingeLine='LE';
    myelevator_left.Label='elevatL';
    myelevator_left.FaceColor='m';
    
    myelevator_left.AeroPanelLength=0.5;
    
    build(myelevator_left)
    Tailwing_left.add(myelevator_left);

    Tailwing_left.ModelControlSurf = 1;


    build(Tailwing_left);


    %% Generate vertical wing and rudder

    Verticalwing=awi.model.LiftingSurface;
    Verticalwing.Name = 'Vertical_wing';
    Verticalwing.ActiveSet = 'pSet';
    Verticalwing.Chord     = [3.31, 1.5];
    Verticalwing.Chord_eta = [0, 1];
    Verticalwing.Span      = 12.45/2;

    Verticalwing.SpanVector = 'Z';
    Verticalwing.Sweep = [30, 30];
    Verticalwing.Dihedral = [0,0];


    all_eta           = Verticalwing.Eta_;
    Verticalwing.BeamLoc     = repmat(0.5, size(all_eta));
    Verticalwing.BeamLoc_eta = all_eta;
    Verticalwing.XOffset=Param.Layout.Vertical_Tail_Position;


    Verticalwing.Material_eta = [0, 1];
    Verticalwing.Material     = [Tail_Mat, Tail_Mat];

    % % Aeropanel definition
    Verticalwing.NumAeroPanel=8;

    % Define box beam corss section
    Verticalbox=awi.model.BoxBeam;
    Verticalbox.BoxType='SymmetricBox';
    Verticalbox.Height=0.5;
    Verticalbox.Width=1;
    Verticalbox.CoverThickness=0.005;
    Verticalbox.SparThickness=0.005;
    getGeometricProps(Verticalbox)
    Verticalwing.BoxBeam = Verticalbox;
    Verticalwing.A   = Verticalbox.Abb;
    Verticalwing.I11 = Verticalbox.Ixx;
    Verticalwing.I22 = Verticalbox.Izz;
    Verticalwing.J   = Verticalbox.Jbb;
%     Verticalwing.NSM = Verticalbox.NSM;
%     Verticalwing.NSI = Verticalbox.NSI;

%     Vertical_tail_mass_val=Vertical_tail/5;
% 
%     for i=1:1:5
%         handle_vertical=strcat('PM_vertical_L','i');
%         handle_vertical=awi.model.PointMass;
%         handle_vertical.SOffset=-0.1+i*0.2;
%         handle_vertical.Mass=Vertical_tail_mass_val;
% %         handle.Inertia11 =  0;
% %         handle.Inertia22 =  0;
% %         handle.Inertia33 =  0;
% %         handle.Inertia12 =  0;
% %         handle.Inertia23 =  0;
% %         handle.Inertia13 =  0;
%         handle_vertical.MassGroup='Group5';
%         Verticalwing.add(handle_vertical);
% 
%     end

    build(Verticalwing);


    %% Build aircraft model
    Aircraft = awi.model.Aircraft;

    Aircraft.add(Body);
    
    
    Body.add(Connector_right)
    Body.add(Connector_left)
    
    Connector_right.add(Wingbox_right)
    Connector_left.add(Wingbox_left)
    
    Body.add(Tailwing_right)
    Body.add(Tailwing_left)
    Body.add(Verticalwing)


    %The analysis methods require an 'awi.model.Aircraft' object
    % This is because some information is only known at the aircraft level,
    % e.g. all-up mass, reference span, area, etc.
    % Aircraft = awi.model.Aircraft;
    % Aircraft.add(LS);
    
    if isfield(Param,'FWT')
        
        Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea, Wingbox_left.SurfaceArea,...
            Connector_right.SurfaceArea,  Connector_left.SurfaceArea, FWT_R.SurfaceArea, FWT_L.SurfaceArea]);
        
        
        Taper_ratio=FWT_R.Chord(end)/Wingbox_right.Chord(1);
        Mean_cord_coefficient=(2/3)*(1+Taper_ratio+Taper_ratio^2)/(1+Taper_ratio);
        
        
        Aircraft.RefSpan  = (Wingbox_right.Span + Connector_right.Span + FWT_R.Span)*2;
        Aircraft.RefChord = Wingbox_right.RootChord*Mean_cord_coefficient; %mean aerodynamic chord = 0.697 for A321 wing;
        
    else
        
        Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea, Wingbox_left.SurfaceArea,...
            Connector_right.SurfaceArea,  Connector_left.SurfaceArea]);
        
        
        Taper_ratio=Wingbox_right.Chord(end)/Wingbox_right.Chord(1);
        Mean_cord_coefficient=(2/3)*(1+Taper_ratio+Taper_ratio^2)/(1+Taper_ratio);
        
        
        Aircraft.RefSpan  = Wingbox_right.Span*2+Connector_right.Span*2;
        Aircraft.RefChord = Wingbox_right.RootChord*Mean_cord_coefficient; %mean aerodynamic chord = 0.697 for A321 wing;
        
    end

    %% Generate the FEM 
% 
%     run_folder = [
%         'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\full_A321_sizing_test3']; %[-], folder for exporting the NASTRAN model

    % Convert to a finite element model and draw it
    FEM_full = convertToFE(Aircraft);
    FEM_full.updateDmiEntry();
    
    
    

end