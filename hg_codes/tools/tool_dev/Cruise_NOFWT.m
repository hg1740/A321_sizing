

function [Aircraft, FEM_full, Y_Data, Y_all, Load_distribution, Lift_all, Cl, Cdi]=Cruise_NOFWT(run_folder,x,Semi_span,Root_chord,LE_sweep,TE_sweep1,TE_sweep2,BeamLoc,Fuselage_length, Wing_position, Engine_position, Horizontal_tail_position, Vertical_tail_position, Secondary_mass,Fuel_mass,Fuselage_total_mass,fold_angle,flare_angle,fold_eta,hinge_stiffness,FWT_root_width,FWT_root_height,FWT_tip_width,FWT_tip_height, FWT_thickness)

    
    % calculate mean aerodynamic chord
    Mid_chord=0.63685*Root_chord;
    Tip_chord=0.2248*Root_chord;
    Taper_ratio=Tip_chord/Root_chord;
    Mean_cord_coefficient=(2/3)*(1+Taper_ratio+Taper_ratio^2)/(1+Taper_ratio);
    
    % fixed masses
    Engine_mass=7362/2; % kg
    Pylon=1239/2; % kg
    Horizontal_tail=682; % kg
    Vertical_tail=522; % kg
   
     %% Right root connector

    Connector_right = awi.model.LiftingSurface;
    Connector_right.Name = 'Connector_Right';
    Connector_right.Origin=[Wing_position,0,0]; %15

    %Use the Leading/Trailing edge sweep to define the planform
    Connector_right.ActiveSet = 'sSet';

    %Tail wing dimensions
    Connector_right.SpanVector  = 'Y';
    Connector_right.Span        = 2;
    Connector_right.LESweep     = [0, 0];
    Connector_right.LESweep_eta = [0, 1];
    Connector_right.TESweep     = [0,  0];
    Connector_right.TESweep_eta = [0,  1];
    Connector_right.RootChord   = Root_chord;

    %Make sure the beam is at the midchord
    all_eta           = Connector_right.Eta_;
    Connector_right.BeamLoc     = repmat(BeamLoc, size(all_eta));
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
%     Connector_right.AeroPanelLength=0.5;
    Connector_right.NumAeroPanel=10;

    build(Connector_right);
    
    
    %% Left root connector
    
    Connector_left = awi.model.LiftingSurface;
    Connector_left.Name = 'Connector_Left';
    Connector_left.Origin=[Wing_position,0,0]; %15
    
    %Use the Leading/Trailing edge sweep to define the planform
    Connector_left.ActiveSet = 'sSet';

    %Tail wing dimensions
    Connector_left.SpanVector  = 'Y';
    Connector_left.Span        = -2;
    Connector_left.LESweep     = [0, 0];
    Connector_left.LESweep_eta = [0, 1];
    Connector_left.TESweep     = [0,  0];
    Connector_left.TESweep_eta = [0,  1];
    Connector_left.RootChord   = Root_chord;

    %Make sure the beam is at the midchord
    all_eta           = Connector_left.Eta_;
    Connector_left.BeamLoc     = repmat(BeamLoc, size(all_eta));
    Connector_left.BeamLoc_eta = all_eta;
%     Connector_right.XOffset=35;
%     Tailwing_right.YOffset=1;

    % Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_root_left = awi.model.Spar;
    FrontSpar_root_left.XLoc = [0.15, 0.15];
    FrontSpar_root_left.Eta  = [0   , 1];
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
%     Connector_left.AeroPanelLength=0.5;
    Connector_left.NumAeroPanel=10;

    build(Connector_left);
   
    
  %% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[Wing_position,2,0]; %15
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';
    
    % Num of element 
    Wingbox_right.NumBeamElem = 23;

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = Semi_span;   %34.1/2;
    Wingbox_right.LESweep     = [LE_sweep, LE_sweep];
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [TE_sweep1, TE_sweep2, TE_sweep2];
    Wingbox_right.TESweep_eta = [0, 0.27, 1];
    Wingbox_right.RootChord   = Root_chord;
    
    
    %Dihedral 
    Wingbox_right.Dihedral=[5,5];
    Wingbox_right.Dihedral_eta=[0,1];
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_right.Eta_;
    Wingbox_right.BeamLoc     = repmat(BeamLoc, size(all_eta));

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
    
    %% Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
    Wingbox_right.NumAeroPanel=10;
    Wingbox_right.AeroPanelAR=2.5;
    
    %     Wingbox_right.AeroPanelLength=0.4;
    
    %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=Wingbox_right.NumBeamElem+2;
    kink_eta=0.27;
    
    %%sizing variables ---------------------------------------------
    
    thickness1=x(1:NumSec);
    thickness2=x(NumSec+1:NumSec*2);
    Astrg=x(NumSec*2+1:NumSec*3);
        
    d_strg=sqrt(Astrg/0.36);
    t_strg=0.12*d_strg;
    % -------------------------------------------------------------

    % etaS=linspace(0,Wingbox_right.Span,NumSec);

    % set width and height array 
    YData=Wingbox_right.YData;
    SparWidth=Wingbox_right.Chord*0.5;

    RootH=Wingbox_right.Chord(1)*0.15; % root thickness/chord = 0.15
    MidH=Wingbox_right.Chord(2)*0.12;  % middle thickness/chord = 0.12
    TipH=Wingbox_right.Chord(end)*0.11;% tip thickness/chord = 0.11


    % set up eta values
    elnum=Wingbox_right.NumBeamElem + 1; % total number of beam elements along the wing
    Num_seg1=ceil(elnum*kink_eta); % number of elements in the inboard section
    Num_seg2=elnum - Num_seg1; % number of elements in the outboard section
    
    Num_sec1=Num_seg1+1;    % number of sections in the inboard section
    Num_sec2=Num_seg2+1;    % number of sections in the outboard section
    
    eta1_=linspace(0,kink_eta, Num_sec1);
    eta2_=linspace(kink_eta,1,Num_sec2);
    etaS=[eta1_(1:end-1),eta2_(1:end)];

    RData=Wingbox_right.RData;
    eta_R=RData/RData(end);
    eta_Y=YData/YData(end);
    etaRS=interp1(eta_Y,eta_R,etaS);

    Bwidth=interp1(RData/RData(end),SparWidth,etaRS);
    Bheight=interp1(RData/RData(end),0.79*[RootH,MidH,TipH],etaRS);
%      Bheight=interp1(RData/RData(end),0.79*[RootH,RootH,RootH],etaRS);

    % stringer pitch 
    strg_n=0.24;

    %intialise data array
    A_val=zeros(1,NumSec);
    Ixx_val=zeros(1,NumSec);
    Izz_val=zeros(1,NumSec);
    J_val=zeros(1,NumSec);
    
    %NSM
    NSM_val=zeros(1,NumSec);
    NSI_val=zeros(1,NumSec);
    
    %offset from shear center
    SCy_val=zeros(1,NumSec);
    SCz_val=zeros(1,NumSec);
    NAy_val=zeros(1,NumSec);
    NAz_val=zeros(1,NumSec);
    CMy_val=zeros(1,NumSec);
    CMz_val=zeros(1,NumSec);
    
    %offset from CoG
    xOff=interp1(YData/YData(end),Wingbox_right.XData,etaS);
    xOff_1=xOff(1:end-1);
    xOff_2=xOff(2:end);
    xOff_val=[0,xOff_2-xOff_1];
    

    for ii=1:NumSec

        boxname=strcat('Box',string(ii));
        boxname=awi.model.BoxBeam;
        boxname.BoxType='SymmetricBox';
        boxname.Height=Bheight(ii);
        boxname.Width=Bwidth(ii);
        boxname.CoverThickness=thickness2(ii);
        boxname.SparThickness=thickness1(ii);

        NumStrg=floor(Bwidth(ii)/strg_n);

        ts=t_strg(ii);
        ds=d_strg(ii);
        hs=Bheight(ii);
        ws=Bwidth(ii);

        Istrg_xx=((ts*ds^3/12)+(ds*ts^3/12 + ds*ts*(ds/2)^2)*2 + 3*ds*ts*(hs/2-ds/2)^2)*NumStrg*2;
        Istrg_zz_=(ds*ts^3/12 + (ts*ds^3/12 + ts*ds*(ds/2)^2)*2);

        if mod(NumStrg,2)==0
            offset=0.12:strg_n:ws/2;
            Istrg_zz=(NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;

        elseif mod(NumStrg,2)==1
            offset=0:strg_n:ws/2;
            Istrg_zz=(NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;

        end

        getGeometricProps(boxname)
        A_val(ii)=boxname.Abb+0;
        Ixx_val(ii)=boxname.Ixx+Istrg_xx;
        Izz_val(ii)=boxname.Izz+Istrg_zz;
        J_val(ii)=boxname.Jbb;
        
        % NSM
       
        NSM_val(ii)=boxname.NSM;
        NSI_val(ii)=boxname.NSI;
        
        % offset
        SCy_val(ii)=boxname.xSC;
        SCz_val(ii)=boxname.zSC;
        NAy_val(ii)=boxname.xNA;
        NAz_val(ii)=boxname.zNA;
        CMy_val(ii)=boxname.xCM;
        CMz_val(ii)=boxname.zCM;

    end


    eta_=etaRS;
    Wingbox_right.A   =  A_val;
    Wingbox_right.A_eta=eta_;

    Wingbox_right.I11 = Izz_val;
    Wingbox_right.I11_eta=eta_;

    Wingbox_right.I22 = Ixx_val;
    Wingbox_right.I22_eta = eta_;

    Wingbox_right.J   = J_val;
    Wingbox_right.J_eta= eta_;

    % NSM and NSI
%     Wingbox_right.NSM=NSM_val;
%     Wingbox_right.NSM_eta= eta_;
%     
%     Wingbox_right.NSI=NSI_val;
%     Wingbox_right.NSI_eta= eta_;
    

    build(Wingbox_right)


    %% Mass definition
    

    wingmass_eta=etaRS;
    Mwidth=interp1(RData/RData(end),SparWidth,wingmass_eta);
    
    mass_set=(Secondary_mass+Fuel_mass)*(Mwidth)/sum(Mwidth);
      
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
    
    
%     %Engine location 
%     Y_eng=Wingbox_right.PanelCoords.LE.Y(2);
%     X_eng=(Wingbox_right.PanelCoords.LE.X(2)+Wingbox_right.PanelCoords.TE.X(2))/2;
%     Z_eng=(Wingbox_right.PanelCoords.LE.Z(2)+Wingbox_right.PanelCoords.TE.Z(2))/2;
%     
%     Engine.Origin = [X_eng-3.5, Y_eng, Z_eng];


    % Engine location - user defined
    Y_Engine=Engine_position;
    
    X_Engine=interp1(Wingbox_right.YData,Wingbox_right.XData,Y_Engine);
    Z_Engine=interp1(Wingbox_right.YData,Wingbox_right.ZData,Y_Engine);
    
    Engine.Origin = [X_Engine-Engine.Length+Wing_position, Y_Engine + 2, Z_Engine];

    
    %Make engine material
    E1  = 76e9; %[N/m^2],set as a rigid body
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
%     Engine_Inertia = Mat1.Rho*Engine_A*Engine_radius^2;
    
    Engine.A   = Engine_A;
    Engine.I11 = Engine_I11;
    Engine.I22 = Engine_I22;
    Engine.J   = Engine_J;
    %     Body.NSM = Bodybox.NSM;
%     Engine.NSI = Engine_Inertia;

%     Engine.A   = 0.04432;
%     Engine.I11 = 0.002;
%     Engine.I22 = 0.002;
%     Engine.J   = 0.001636;
    
    %Aeropanel althoufh it is useless now
    Engine.AeroPanelLength=0.5;
    
    % add engine mass
    engine_mass=awi.model.PointMass;   
    engine_mass.SOffset=0.1;
    engine_mass.Mass=Engine_mass;
    Engine.add(engine_mass);
    
    % add pylon
    pylon_mass=awi.model.PointMass;   
    pylon_mass.SOffset=0.9;
    pylon_mass.Mass=Pylon;
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
    Wingbox_left.Origin=[Wing_position,-2,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_left.ActiveSet = 'sSet';

    %Wing dimensions
    Wingbox_left.SpanVector  = 'Y';
    Wingbox_left.Span        = -Semi_span;  
    Wingbox_left.LESweep     = [-LE_sweep, -LE_sweep];
    Wingbox_left.LESweep_eta = [0, 1];
    Wingbox_left.TESweep     = [-TE_sweep1, -TE_sweep2, -TE_sweep2];
    Wingbox_left.TESweep_eta = [0, 0.27, 1];
    Wingbox_left.RootChord   = Root_chord;   
    
    %Dihedral 
    Wingbox_left.Dihedral=[5,5];
    Wingbox_left.Dihedral_eta=[0,1];
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_left.Eta_;
    Wingbox_left.BeamLoc     = repmat(BeamLoc, size(all_eta));
%     Wingbox_right.BeamLoc     = [0.34,0.4,0.4];
    Wingbox_left.BeamLoc_eta = all_eta;

    %Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_left = awi.model.Spar;
    FrontSpar_left.XLoc = [0.15, 0.15];
    FrontSpar_left.Eta  = [0   , 1];
    RearSpar_left = awi.model.Spar;
    RearSpar_left.XLoc = [0.65, 0.65];
    RearSpar_left.Eta  = [0   , 1];

    Wingbox_left.add([FrontSpar_left, RearSpar_left]);

    %Define internal layout
    Wingbox_left.RibPitch      = 0.65;
    Wingbox_left.StringerPitch = 0.15;

%     %Make the material
%     E  = 70e9; %[N/m^2], typical YM of aluminium
%     nu = 0.333;
%     rho=2810;
%     Mat = awi.model.Material;
%     Mat.E  = E;
%     Mat.Nu = nu;
%     Mat.G  = E / (2 * (1 + nu));
%     Mat.Rho=rho;
    Wingbox_left.Material_eta = [0, 1];
    Wingbox_left.Material     = [Mat_wing, Mat_wing];
    
    build(Wingbox_left)
    
    %% Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
    Wingbox_left.NumAeroPanel=10;
    Wingbox_left.AeroPanelAR=2.5;
    %     Wingbox_left.AeroPanelLength=0.4;
    
   
    %% Create discretised boxbeam with varied cross section prperties along the span 

    eta_=etaRS;
    Wingbox_left.A   =  A_val;
    Wingbox_left.A_eta=eta_;

    Wingbox_left.I11 = Izz_val;
    Wingbox_left.I11_eta=eta_;

    Wingbox_left.I22 = Ixx_val;
    Wingbox_left.I22_eta = eta_;

    Wingbox_left.J   = J_val;
    Wingbox_left.J_eta= eta_;

%     % NSM and NSI
%     Wingbox_left.NSM=NSM_val;
%     Wingbox_left.NSM_eta= eta_;
%     
%     Wingbox_left.NSI=NSI_val;
%     Wingbox_left.NSI_eta= eta_;
    
  
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

    Engine2.Origin = [X_Engine-Engine.Length+Wing_position, -(Y_Engine + 2), Z_Engine];
   
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
    engine2_mass.Mass=Engine_mass;
    Engine2.add(engine2_mass);
    
    % add pylon
    pylon2_mass=awi.model.PointMass;   
    pylon2_mass.SOffset=0.9;
    pylon2_mass.Mass=Pylon;
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
    Body.Length=Fuselage_length;
    % Body.XOffset=-15;

    
    %Make the material
    E1  = 76e9; %[N/m^2],set as a rigid body
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
    Mass_val=Fuselage_total_mass/11;

    
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
    Tailwing_right.Dihedral=[5,5];
    Tailwing_right.Dihedral_eta=[0,1];

    %Make sure the beam is at the midchord
    all_eta           = Tailwing_right.Eta_;
    Tailwing_right.BeamLoc     = repmat(0.5, size(all_eta));
    Tailwing_right.BeamLoc_eta = all_eta;
    Tailwing_right.XOffset = Horizontal_tail_position;
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
    Et  = 76e9; %[N/m^2],set as a rigid body
    nut = 0.333;
    Tail_Mat = awi.model.Material;
    Tail_Mat.E  = Et;
    Tail_Mat.Nu = nut;
    Tail_Mat.G  = E1 / (2 * (1 + nut));
    Tail_Mat.Nu = nu;
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

    %Make sure the beam is at the midchord
    all_eta           = Tailwing_left.Eta_;
    Tailwing_left.BeamLoc     = repmat(0.5, size(all_eta));
    Tailwing_left.BeamLoc_eta = all_eta;
    Tailwing_left.XOffset=Horizontal_tail_position;
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
    Verticalwing.XOffset=Vertical_tail_position;


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

    Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea, Wingbox_left.SurfaceArea,...
        Connector_right.SurfaceArea,  Connector_left.SurfaceArea, FWT_R.SurfaceArea, FWT_L.SurfaceArea]);
    
    
    Aircraft.RefSpan  = Wingbox_right.Span*2 + Connector_right.Span*2 + FWT_R.Span*2;
    Aircraft.RefChord = Wingbox_right.RootChord*Mean_cord_coefficient; %mean aerodynamic chord = 0.697 for A321 wing;
%     Aircraft.RefChord = Aircraft.RefArea/Aircraft.RefSpan; 



    %% Generate the FEM 
% 
%     run_folder = [
%         'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\full_A321_sizing_test3']; %[-], folder for exporting the NASTRAN model

    % Convert to a finite element model and draw it
    FEM_full = convertToFE(Aircraft);

    % %Export it to a file
    export(FEM_full, run_folder);



    %% LoadCase 1: cruising altitude level flight
    TrimLoadcase1 = awi.model.LoadCase;
    
    acMass = 94000;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;
    
%     FlightPoint=awi.model.FlightPoint;
%     FlightPoint.Mach=0.78;
%     FlightPoint.AcVelocity=FlightPoint.Mach*340;
%     FlightPoint.Altitude = 36000;
%     getFlightPointData(FlightPoint,'ISA');

    TrimLoadcase1.Name = 'A321_cruise_2p5g';
    TrimLoadcase1.Altitude   = altitude;
    TrimLoadcase1.Mach       = mach_number;
    TrimLoadcase1.AcVelocity = aircraft_velocity;
    TrimLoadcase1.AcMass = acMass;
    

    TrimLoadcase1.PitchAngle=0;
    TrimLoadcase1.RollAngle =0;
    TrimLoadcase1.ID = 1020;
    TrimLoadcase1.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase1.CsDeflection=flap_angle*pi/180;
    
    
    build(TrimLoadcase1)
    

  %% NASTRAN method - RUN SOL 144

     NastranMethods1 = awi.methods.Nastran;
     NastranMethods1.AnalysisModel = FEM_full;
     MassCases=awi.model.MassCases.empty;
     ID0=200;    
     
    % sizing case 1
    trimFile1 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase1, MassCases,run_folder,'DatFilename','A321_36000ft_1g');
    
   
   
    delete(strcat(run_folder, '\*.xdb'));
    delete(strcat(run_folder, '\*.h5'));
    delete(strcat(run_folder, '\*.log'));
    delete(strcat(run_folder, '\*.f06'));
    delete(strcat(run_folder, '\*.f04'));
    delete(strcat(run_folder, '\*.op4'));
    
    delete(strcat(run_folder, '\*.xdb.*'));
    delete(strcat(run_folder, '\*.h5.*'));
    delete(strcat(run_folder, '\*.log.*'));
    delete(strcat(run_folder, '\*.f06.*'));
    delete(strcat(run_folder, '\*.f04.*'));
         
    NastranMethods1.runNastran(trimFile1);
 
    
    %% extract results
 
    Trim1_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_36000ft_1g.f06'),'ReadF06',true,'ReadHDF5',false);

    
    WingNodes=NastranMethods1.WingNode;
    
    index1=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[WingNodes(1:end-1).GID]);    
    index2=ismember([Trim1_res.f06data.Bendingmoment.LGRID],[WingNodes(end).GID]);
    
    % extract forces:
    
    Trim1.M_P1=[Trim1_res.f06data.Bendingmoment.UMPLN1(index1),Trim1_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Trim1.M_P2=[Trim1_res.f06data.Bendingmoment.UMPLN2(index1),Trim1_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Trim1.T=[Trim1_res.f06data.Bendingmoment.UTORQUE1(index1),Trim1_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque

    Trim1.S_P1=[Trim1_res.f06data.Bendingmoment.USPLN1(index1),Trim1_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Trim1.S_P2=[Trim1_res.f06data.Bendingmoment.USPLN2(index1),Trim1_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
   
    % Y-Data
    WingNodes=NastranMethods1.WingNode;
    
    X=[WingNodes.X];

    Y_Data=X(2,:);
    

    %% Internal Load Summery 
    
    Moment_P1_Loadcase1=Trim1.M_P1; % internal load M1 for 2.5 g pull up
    Moment_P2_Loadcase1=Trim1.M_P2; % internal load M2 for 2.5 g pull up
    Torque_Loadcase1=Trim1.T;       % internal load T for 2.5 g pull up
    Shear_P1_Loadcase1=Trim1.S_P1;
    Shear_P2_Loadcase1=Trim1.S_P2;
    
    %% Result output
    
    % Out of plane bending moment distributions for all load cases
    Load_distribution.Moment_P2.LC1=Moment_P2_Loadcase1;

    % Vertical shear force distributions for all load cases
    Load_distribution.Shear_P2.LC1=Shear_P2_Loadcase1;
    
    % Torque distribution for all load cases 
    Load_distribution.Torque.LC1=Torque_Loadcase1;
    

    
    %% Drag 
    
    % create an object array 
    All_FEM=flatlist(FEM_full);
    
    Wing_conn=All_FEM(ismember({All_FEM.Name}, 'Connector_Right'));
    Wing_inboard = All_FEM(ismember({All_FEM.Name}, 'A320Wing_right'));
    
    Wing_conn_ID=Wing_conn.AeroSets.E;
    Wing_inboard_ID=Wing_inboard.AeroSets(1).E;
    Wing_outboard_ID=Wing_inboard.AeroSets(2).E;
    
    
    % aeroforces at each segment
    res_aeroF = mni.result.f06(strcat('C:\Git\A321_sizing\hg_codes\results\extension10','/A321_36000ft_1g.f06')).read_aeroF;
    
    Index_conn=ismember([res_aeroF.PanelID(1,:)],Wing_conn_ID);
    Index_wing_inboard=ismember([res_aeroF.PanelID(1,:)],Wing_inboard_ID);
    Index_wing_outboard=ismember([res_aeroF.PanelID(1,:)],Wing_outboard_ID);
    
    Lift_conn=res_aeroF.aeroFz(Index_conn);
    Lift_wing_inboard=res_aeroF.aeroFz(Index_wing_inboard);
    Lift_wing_outboard=res_aeroF.aeroFz(Index_wing_outboard);
    
    
    % calculate panel width for each segment: conn + wing_bef_kink +
    % wing_aft_kink 
    
    num_panel_conn=numel(Wing_conn_ID)/10;
    num_panel_inboard=numel(Wing_inboard_ID)/10;
    num_panel_outboard=numel(Wing_outboard_ID)/10;
    

    kink_pos=Wingbox.YData(2);
    
    width_conn=2/num_panel_conn;
    
    width_wing1=kink_pos/num_panel_inboard;
    width_wing2=(Wingbox.YData(3)-Wingbox.YData(2))/num_panel_outboard;
    
    panel_width =[width_conn*ones(1,num_panel_conn),width_wing1*ones(1,num_panel_inboard),width_wing2*ones(1,num_panel_outboard)];

    
    % Y-Data
    % find corresponding Y positions
    Y_conn=0.5*width_conn:width_conn*0.999:2;
    
    Y_wing1=2+0.5*width_wing1:width_wing1*0.999:2+Wingbox.YData(2);
    
    Y_wing2=2+Wingbox.YData(2)+ 0.5*width_wing2:width_wing2*0.999:2+Wingbox.YData(3);
    
    Y=[Y_conn,Y_wing1,Y_wing2];
    
    % normalise lift by panel width to obtain lift per unit span
    lift_wing=[Lift_conn';Lift_wing_inboard';Lift_wing_outboard'];
    
    lift_wing_matrix=reshape(lift_wing, 10, numel(lift_wing)/10);
    
    % abs. value of the lift
    lift_wing_abs=sum(lift_wing_matrix);
    
    Cum_num=cumsum([num_panel_conn,num_panel_inboard,num_panel_outboard,num_panel_FWT]);
    
    lift_wing_matrix(:,1:Cum_num(1))=lift_wing_matrix(:,1:Cum_num(1))/width_conn;
    lift_wing_matrix(:,Cum_num(1)+1:Cum_num(2))=lift_wing_matrix(:,Cum_num(1)+1:Cum_num(2))/width_wing1;
    lift_wing_matrix(:,Cum_num(2)+1:Cum_num(3))=lift_wing_matrix(:,Cum_num(2)+1:Cum_num(3))/width_wing2;

    
    % lift per unit span
    lift_wing_var=sum(lift_wing_matrix);
    
    % whole wing span lift distribution
    Y_left=sort(-Y);
    lift_wing_var_left=flip(lift_wing_var);
    lift_wing_abs_left=flip(lift_wing_abs);
    panel_width_left=flip(panel_width);
    
    Y_all=[Y_left,Y]';
    Lift_all=[lift_wing_var_left,lift_wing_var];
    Lift_all_abs=[lift_wing_abs_left,lift_wing_abs];
    
    panel_width_all=[panel_width_left,panel_width];
    
    % calculate vorticity
    % Air conditions
    FlightPoint=awi.model.FlightPoint;
    FlightPoint.Mach=0.78;
    FlightPoint.AcVelocity=0.78*340;
    FlightPoint.Altitude = 36000;
    getFlightPointData(FlightPoint,'ISA');
    DynPressure = FlightPoint.DynPressure;
    
    Gamma_all=Lift_all/(FlightPoint.AirDensity*FlightPoint.AcVelocity);
    
    % calculate the derivative
    dGdy= gradient(Gamma_all(:)) ./ gradient(Y_all(:));
    
    % downwash
    wj=zeros(1,length(Y_all));
    alphai=zeros(1,length(Y_all));
    
    
    for i=1:length(Y_all)
        
        w=-(1/(4*pi))*dGdy* panel_width_all(i)./(Y_all(i)-Y_all);
        
        w=w( ~any( isnan( w ) | isinf( w ), 2 ),: );
        
        wj(i)=sum(w);
        
        alphai(i)=wj(i)/FlightPoint.AcVelocity;
        
    end
    
    Dragi_var=Lift_all_abs.*sin(alphai);
    Drag_force=sum([Dragi_var]);
    
    Total_area=(Wingbox_right.SurfaceArea+FWT_R.SurfaceArea+Connector_right.SurfaceArea)*2;
    
    Cdi=Drag_force/(DynPressure*Total_area);
    
    Cl=sum(lift_wing_abs)/(DynPressure*Total_area);
    
 
end
    




























