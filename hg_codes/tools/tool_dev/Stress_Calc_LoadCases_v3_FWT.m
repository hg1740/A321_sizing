
%% SOL 144 + SOL 146

%load cases:

% Case 1: 2.5 g 36000ft M =0.78
% Case 2: 2.5 g 3000ft M=0.48
% Case 3: g + 1MC 36000ft M = 0.78
% Case 4: g + 1MC 20000ft M = 0.6
% Case 5: g + 1MC 3000ft M = 0.48
% Case 6: g - 1MC 36000ft M = 0.78 
% Case 7: g - 1MC 20000ft M = 0.6
% Case 8: g - 1MC 3000ft M = 0.48

% -1MC results can be obtained by just take the mirror image of the results from +1MC about horizontal axis

function [Aircraft, FEM_full, Y_all, Load_distribution, Delta, Von_skn, Von_spr, sigmab_skn, taub_skn, sigma_pp, sigma_strg, sigma_crip, sigma_col]=Stress_Calc_LoadCases_v3_FWT(run_folder,x,Semi_span,Root_chord,LE_sweep,TE_sweep1,TE_sweep2,BeamLoc,Fuselage_length, Wing_position, Engine_position, Horizontal_tail_position, Vertical_tail_position, Secondary_mass,Fuel_mass,Fuselage_total_mass,fold_angle,flare_angle,fold_eta,hinge_stiffness,FWT_root_width,FWT_root_height,FWT_tip_width,FWT_tip_height, FWT_thickness)

    
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
    Connector_left.Name = 'Connector_Right';
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
    
    
    %% initialise properties 
    
    Wingbox_right.A   =  [1,1];
    Wingbox_right.A_eta=[0,1];
    
    Wingbox_right.I11 = [1,1];
    Wingbox_right.I11_eta=[0,1];
    
    Wingbox_right.I22 = [1,1];
    Wingbox_right.I22_eta = [0,1];
    
    Wingbox_right.J   = [1,1];
    Wingbox_right.J_eta= [0,1];
    
    %% FWT right definition
    
%     fold_angle  = -10;   %[deg],
%     flare_angle = 25;   %[deg],
%     fold_eta=0.75;
    
    
    FWT_R = insertWingFold(Wingbox_right, 'FlareAngle', flare_angle, 'FoldAngle', fold_angle,'EtaFold',fold_eta);
    FWT_R.HingeStiffness = [1e14 1e14 1e14 1e14 hinge_stiffness 1e14];
    
    FWT_R.NumAeroPanel=10;
    FWT_R.AeroPanelAR=2.5;
    
    FWT_R.NumBeamElem = 10;
    
    FWT_box_root=awi.model.BoxBeam;
    FWT_box_root.BoxType='SymmetricBox';
    FWT_box_root.Height=FWT_root_height;
    FWT_box_root.Width=FWT_root_width;
    
    FWT_box_root.CoverThickness=FWT_thickness;
    FWT_box_root.SparThickness=FWT_thickness;
    
    getGeometricProps(FWT_box_root)
    
    FWT_box_tip=awi.model.BoxBeam;
    FWT_box_tip.BoxType='SymmetricBox';
    FWT_box_tip.Height=FWT_tip_height;
    FWT_box_tip.Width=FWT_tip_width;
    
    FWT_box_tip.CoverThickness=FWT_thickness;
    FWT_box_tip.SparThickness=FWT_thickness;
    
    getGeometricProps(FWT_box_tip)
    
    FWT_eta=[0,1];
    
    FWT_R.A   =  [FWT_box_root.Abb, FWT_box_tip.Abb];
    FWT_R.A_eta=FWT_eta;
    
    FWT_R.I11   =  [FWT_box_root.Izz, FWT_box_tip.Izz];
    FWT_R.I11_eta=FWT_eta;
    
    FWT_R.I22   =  [FWT_box_root.Ixx, FWT_box_tip.Ixx];
    FWT_R.I22_eta=FWT_eta;
    
    FWT_R.J   =  [FWT_box_root.Jbb, FWT_box_tip.Jbb];
    FWT_R.J_eta=FWT_eta;
    
%     sc_num=numel(FWT_R.A);
% %     FWT_eta=[0,1];
%     
%     FWT_R.A   =  ones(1,sc_num).*FWT_box.Abb;
%     % FWT.A_eta=FWT_eta;
%     
%     FWT_R.I11 = ones(1,sc_num).*FWT_box.Izz;
%     % FWT.I11_eta=FWT_eta;
%     
%     FWT_R.I22 = ones(1,sc_num).*FWT_box.Ixx;
%     % FWT.I22_eta = FWT_eta;
%     
%     FWT_R.J   = ones(1,sc_num).*FWT_box.Jbb;
%     % FWT.J_eta= FWT_eta;
    
    
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
    
    
    % add point masses
    for i=1:1:2
        handle=strcat('PM_right','i');
        handle=awi.model.PointMass;
        handle.SOffset=0+i*0.2;
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
    
  
    %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=Wingbox_right.NumBeamElem+2;
    kink_eta=0.27*(1/fold_eta);
    
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
    
    
    %% initialise properties
    
    Wingbox_left.A   =  1;
    Wingbox_left.A_eta=[0,1];
    
    Wingbox_left.I11 = 1;
    Wingbox_left.I11_eta=[0,1];
    
    Wingbox_left.I22 = 1;
    Wingbox_left.I22_eta = [0,1];
    
    Wingbox_left.J   = 1;
    Wingbox_left.J_eta= [0,1];
    
    %% FWT left definition
    
    FWT_L = insertWingFold(Wingbox_left, 'FlareAngle', -flare_angle, 'FoldAngle', fold_angle,'EtaFold',fold_eta);
    FWT_L.HingeStiffness = [1e14 1e14 1e14 1e14 hinge_stiffness 1e14];
    
    FWT_L.NumAeroPanel=10;
    FWT_L.AeroPanelAR=2.5;
    
    FWT_L.NumBeamElem = 10;
    
    FWT_L.A   =  [FWT_box_root.Abb, FWT_box_tip.Abb];
    FWT_L.A_eta=FWT_eta;
    
    FWT_L.I11   =  [FWT_box_root.Izz, FWT_box_tip.Izz];
    FWT_L.I11_eta=FWT_eta;
    
    FWT_L.I22   =  [FWT_box_root.Ixx, FWT_box_tip.Ixx];
    FWT_L.I22_eta=FWT_eta;
    
    FWT_L.J   =  [FWT_box_root.Jbb, FWT_box_tip.Jbb];
    FWT_L.J_eta=FWT_eta;
   
%     FWT_L.A   =  ones(1,sc_num).*FWT_box.Abb;
%     % FWT.A_eta=FWT_eta;
%     
%     FWT_L.I11 = ones(1,sc_num).*FWT_box.Izz;
%     % FWT.I11_eta=FWT_eta;
%     
%     FWT_L.I22 = ones(1,sc_num).*FWT_box.Ixx;
%     % FWT.I22_eta = FWT_eta;
%     
%     FWT_L.J   = ones(1,sc_num).*FWT_box.Jbb;
%     % FWT.J_eta= FWT_eta;
    
    %Make the material for FWT
    
    FWT_L.Material_eta = [0, 1];
    FWT_L.Material     = [Mat_fwt, Mat_fwt];
        
    % add point masses
    for i=1:1:2
        handle=strcat('PM_right','i');
        handle=awi.model.PointMass;
        handle.SOffset=0+i*0.2;
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


    %% LoadCase 1: cruising altitude 2.5g pull up
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
    TrimLoadcase1.LoadFactor = 2.5;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase1.CsDeflection=flap_angle*pi/180;
    
    
    build(TrimLoadcase1)
    
    
    %% Loadcase 2: sea level 2.5g pull up
    TrimLoadcase2 = awi.model.LoadCase;

    acMass = 500;
    altitude          = 3000;
    mach_number       = 0.48;
    aircraft_velocity = 0.48*340;
    flap_angle=0;

    TrimLoadcase2.Name = 'A321_sea_level_2.5g';

    TrimLoadcase2.Altitude   = altitude;
    TrimLoadcase2.Mach       = mach_number;
    TrimLoadcase2.AcVelocity = aircraft_velocity;
    TrimLoadcase2.AcMass = acMass;

    TrimLoadcase2.PitchAngle=0;
    TrimLoadcase2.RollAngle =0;
    TrimLoadcase2.ID = 1030;
    TrimLoadcase2.LoadFactor = 2.5;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase2.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase2)
    
    %% Loadcase 3: sea level - g sea level dive
    TrimLoadcase3 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 3000;
    mach_number       = 0.48;
    aircraft_velocity = 0.48*340;
    flap_angle=0;

    TrimLoadcase3.Name = 'A321_sea_level_neg_g';

    TrimLoadcase3.Altitude   = altitude;
    TrimLoadcase3.Mach       = mach_number;
    TrimLoadcase3.AcVelocity = aircraft_velocity;
    TrimLoadcase3.AcMass = acMass;

    TrimLoadcase3.PitchAngle=0;
    TrimLoadcase3.RollAngle =0;
    TrimLoadcase3.ID = 1030;
    TrimLoadcase3.LoadFactor = -1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase3.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase3)
    
    %% steady state for gust analysis, 36000 ft, g
    TrimLoadcase4 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;

    TrimLoadcase4.Name = 'A321_36000ft_g';

    TrimLoadcase4.Altitude   = altitude;
    TrimLoadcase4.Mach       = mach_number;
    TrimLoadcase4.AcVelocity = aircraft_velocity;
    TrimLoadcase4.AcMass = acMass;

    TrimLoadcase4.PitchAngle=0;
    TrimLoadcase4.RollAngle =0;
    TrimLoadcase4.ID = 1030;
    TrimLoadcase4.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase4.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase4)
    
    %% steady state for gust analysis, 20000 ft, g
    TrimLoadcase5 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 20000;
    mach_number       = 0.6;
    aircraft_velocity = 0.6*340;
    flap_angle=0;
    
    TrimLoadcase5.Name = 'A321_20000ft_g';
    
    TrimLoadcase5.Altitude   = altitude;
    TrimLoadcase5.Mach       = mach_number;
    TrimLoadcase5.AcVelocity = aircraft_velocity;
    TrimLoadcase5.AcMass = acMass;
    
    TrimLoadcase5.PitchAngle=0;
    TrimLoadcase5.RollAngle =0;
    TrimLoadcase5.ID = 1030;
    TrimLoadcase5.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase5.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase5)
    
    %% steady state for gust analysis, 3000 ft, g
    TrimLoadcase6 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 3000;
    mach_number       = 0.48;
    aircraft_velocity = 0.48*340;
    flap_angle=0;
    
    TrimLoadcase6.Name = 'A321_3000ft_g';
    
    TrimLoadcase6.Altitude   = altitude;
    TrimLoadcase6.Mach       = mach_number;
    TrimLoadcase6.AcVelocity = aircraft_velocity;
    TrimLoadcase6.AcMass = acMass;
    
    TrimLoadcase6.PitchAngle=0;
    TrimLoadcase6.RollAngle =0;
    TrimLoadcase6.ID = 1030;
    TrimLoadcase6.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase6.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase6)

    %% Generate the FEM 
% 
%     run_folder = [
%         'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\full_A321_sizing_test3']; %[-], folder for exporting the NASTRAN model

    % Convert to a finite element model and draw it
    FEM_full = convertToFE(Aircraft);

    % %Export it to a file
    export(FEM_full, run_folder);

  %% NASTRAN method - RUN SOL 144

     NastranMethods1 = awi.methods.Nastran;
     NastranMethods1.AnalysisModel = FEM_full;
     MassCases=awi.model.MassCases.empty;
     ID0=200;    
     
    % sizing case 1
    trimFile1 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase1, MassCases,run_folder,'DatFilename','A321_cruise_2p5g');
    
    % sizing case 2
    trimFile2 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase2, MassCases,run_folder,'DatFilename','A321_sealevel_2p5g');
    
    % sizing case 3
    trimFile3 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase3, MassCases,run_folder,'DatFilename','A321_sealevel_neg_1g');
    
    % sizing case 4
    trimFile4 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase4, MassCases,run_folder,'DatFilename','A321_36000ft_1g');
    
    % sizing case 5
    trimFile5 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase5, MassCases,run_folder,'DatFilename','A321_20000ft_1g');
    
    % sizing case 6
    trimFile6 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase6, MassCases,run_folder,'DatFilename','A321_3000ft_1g');
    
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\A320_half_model_SOL144*.*')

   
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
    NastranMethods1.runNastran(trimFile2);
    NastranMethods1.runNastran(trimFile3);
    NastranMethods1.runNastran(trimFile4);
    NastranMethods1.runNastran(trimFile5);
    NastranMethods1.runNastran(trimFile6);
    
    
    %% Gust analysis Case 1 : crusing altitude + 1MC
    
    GustLoadcase1 = awi.model.LoadCase;
    %     GustLoadcase.LoadCaseType = 'Pratt Gust';
    GustLoadcase1.Altitude   = 36000;
    GustLoadcase1.AcVelocity = 0.78*340;
    GustLoadcase1.AcMass = 500;
    GustLoadcase1.Mach = 0.78;
    GustLoadcase1.GustLength = linspace(18,214,7);%[18:20:214];
    
    % Gust direction: positive or negative hit 
    GustLoadcase1.GustDirection=1;
    
    FlightPoint1=awi.model.FlightPoint;
    FlightPoint1.Mach=0.78;
    FlightPoint1.AcVelocity=FlightPoint1.Mach*340;
    FlightPoint1.Altitude = 36000;
    getFlightPointData(FlightPoint1,'ISA');
    
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_full;
    MassCases=awi.model.MassCases.empty;
    
    gustfile1=NastranMethods1.writeGustFile(Aircraft, GustLoadcase1, MassCases, FlightPoint1, run_folder,'DatFilename','gust_analysis_g_36000ft_pos_1MC');

    NastranMethods1.runNastran(gustfile1);
    
    %% Gust analysis Case 2 : 20000ft altitude,  + 1MC
    
    GustLoadcase2 = awi.model.LoadCase;
    %     GustLoadcase.LoadCaseType = 'Pratt Gust';
    GustLoadcase2.Altitude   = 20000;
    GustLoadcase2.AcVelocity = 0.6*340;
    GustLoadcase2.AcMass = 500;
    GustLoadcase2.Mach = 0.6;
    GustLoadcase2.GustLength = linspace(18,214,7);%[18:20:214];
    
    % Gust direction: positive or negative hit
    GustLoadcase2.GustDirection=1;
    
    FlightPoint2=awi.model.FlightPoint;
    FlightPoint2.Mach=0.6;
    FlightPoint2.AcVelocity=FlightPoint2.Mach*340;
    FlightPoint2.Altitude = 20000;
    getFlightPointData(FlightPoint2,'ISA');
    
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_full;
    MassCases=awi.model.MassCases.empty;
    
    gustfile2=NastranMethods1.writeGustFile(Aircraft, GustLoadcase2, MassCases, FlightPoint2, run_folder,'DatFilename','gust_analysis_g_20000ft_pos_1MC');
    
    NastranMethods1.runNastran(gustfile2);
    
    %% Gust analysis Case 3 : sea level,  + 1MC
    
    GustLoadcase3 = awi.model.LoadCase;
    %     GustLoadcase.LoadCaseType = 'Pratt Gust';
    GustLoadcase3.Altitude   = 3000;
    GustLoadcase3.AcVelocity = 0.48*340;
    GustLoadcase3.AcMass = 500;
    GustLoadcase3.Mach = 0.48;
    GustLoadcase3.GustLength = linspace(18,214,7);%[18:20:214];
    
    % Gust direction: positive or negative hit
    GustLoadcase3.GustDirection=1;
    
    FlightPoint3=awi.model.FlightPoint;
    FlightPoint3.Mach=0.48;
    FlightPoint3.AcVelocity=FlightPoint3.Mach*340;
    FlightPoint3.Altitude = 3000;
    getFlightPointData(FlightPoint3,'ISA');
    
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_full;
    MassCases=awi.model.MassCases.empty;
    
    gustfile3=NastranMethods1.writeGustFile(Aircraft, GustLoadcase3, MassCases, FlightPoint3, run_folder,'DatFilename','gust_analysis_g_3000ft_pos_1MC');
    
    NastranMethods1.runNastran(gustfile3);
   
    
    
    %% extract results
 
    Trim1_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_cruise_2p5g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim2_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_sealevel_2p5g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim3_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_sealevel_neg_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim4_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_36000ft_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim5_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_20000ft_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim6_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_3000ft_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    
    WingNodes=NastranMethods1.WingNode;
    FWTNodes=NastranMethods1.FWTNode;
    
    index1=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[WingNodes(1:end-1).GID]);    
    index2=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[FWTNodes(1:end).GID]);
    index3=ismember([Trim1_res.f06data.Bendingmoment.LGRID],[FWTNodes(end).GID]);
    
    % extract forces:
    
    % forces from trim 1: 2.5 g 36000ft
    Trim1.M_P1=[Trim1_res.f06data.Bendingmoment.UMPLN1(index1),Trim1_res.f06data.Bendingmoment.LMPLN1(index2), Trim1_res.f06data.Bendingmoment.LMPLN1(index3)]; % in-plane moment
    Trim1.M_P2=[Trim1_res.f06data.Bendingmoment.UMPLN2(index1),Trim1_res.f06data.Bendingmoment.LMPLN2(index2), Trim1_res.f06data.Bendingmoment.UMPLN2(index3)]; % out of plane moment
    Trim1.T=[Trim1_res.f06data.Bendingmoment.UTORQUE1(index1),Trim1_res.f06data.Bendingmoment.LTORQUE1(index2),Trim1_res.f06data.Bendingmoment.UTORQUE1(index3)];% torque

    Trim1.S_P1=[Trim1_res.f06data.Bendingmoment.USPLN1(index1),Trim1_res.f06data.Bendingmoment.LSPLN1(index2), Trim1_res.f06data.Bendingmoment.USPLN1(index3)]; % in plane shear
    Trim1.S_P2=[Trim1_res.f06data.Bendingmoment.USPLN2(index1),Trim1_res.f06data.Bendingmoment.LSPLN2(index2), Trim1_res.f06data.Bendingmoment.USPLN2(index3)]; % out of plane shear
    
    % forces from trim 2: 2.5g Sea level
    Trim2.M_P1=[Trim2_res.f06data.Bendingmoment.UMPLN1(index1),Trim2_res.f06data.Bendingmoment.LMPLN1(index2), Trim2_res.f06data.Bendingmoment.UMPLN1(index3)]; % in-plane moment
    Trim2.M_P2=[Trim2_res.f06data.Bendingmoment.UMPLN2(index1),Trim2_res.f06data.Bendingmoment.LMPLN2(index2), Trim2_res.f06data.Bendingmoment.UMPLN2(index3)]; % out of plane moment
    Trim2.T=[Trim2_res.f06data.Bendingmoment.UTORQUE1(index1),Trim2_res.f06data.Bendingmoment.LTORQUE1(index2), Trim2_res.f06data.Bendingmoment.UTORQUE1(index3)];% torque

    Trim2.S_P1=[Trim2_res.f06data.Bendingmoment.USPLN1(index1),Trim2_res.f06data.Bendingmoment.LSPLN1(index2), Trim2_res.f06data.Bendingmoment.USPLN1(index3)]; % in plane shear
    Trim2.S_P2=[Trim2_res.f06data.Bendingmoment.USPLN2(index1),Trim2_res.f06data.Bendingmoment.LSPLN2(index2),Trim2_res.f06data.Bendingmoment.USPLN2(index3)]; % out of plane shear
    
    % forces from trim 3: -g sea level
    Trim3.M_P1=[Trim3_res.f06data.Bendingmoment.UMPLN1(index1),Trim3_res.f06data.Bendingmoment.LMPLN1(index2),Trim3_res.f06data.Bendingmoment.UMPLN1(index3)]; % in-plane moment
    Trim3.M_P2=[Trim3_res.f06data.Bendingmoment.UMPLN2(index1),Trim3_res.f06data.Bendingmoment.LMPLN2(index2),Trim3_res.f06data.Bendingmoment.UMPLN2(index3)]; % out of plane moment
    Trim3.T=[Trim3_res.f06data.Bendingmoment.UTORQUE1(index1),Trim3_res.f06data.Bendingmoment.LTORQUE1(index2),Trim3_res.f06data.Bendingmoment.UTORQUE1(index3)];% torque
    
    Trim3.S_P1=[Trim3_res.f06data.Bendingmoment.USPLN1(index1),Trim3_res.f06data.Bendingmoment.LSPLN1(index2),Trim3_res.f06data.Bendingmoment.USPLN1(index3)]; % in plane shear
    Trim3.S_P2=[Trim3_res.f06data.Bendingmoment.USPLN2(index1),Trim3_res.f06data.Bendingmoment.LSPLN2(index2), Trim3_res.f06data.Bendingmoment.USPLN2(index3)]; % out of plane shear
    
    % forces from trim 4: g 36000ft
    Trim4.M_P1=[Trim4_res.f06data.Bendingmoment.UMPLN1(index1),Trim4_res.f06data.Bendingmoment.LMPLN1(index2), Trim4_res.f06data.Bendingmoment.UMPLN1(index3)]; % in-plane moment
    Trim4.M_P2=[Trim4_res.f06data.Bendingmoment.UMPLN2(index1),Trim4_res.f06data.Bendingmoment.LMPLN2(index2), Trim4_res.f06data.Bendingmoment.UMPLN2(index3)]; % out of plane moment
    Trim4.T=[Trim4_res.f06data.Bendingmoment.UTORQUE1(index1),Trim4_res.f06data.Bendingmoment.LTORQUE1(index2), Trim4_res.f06data.Bendingmoment.UTORQUE1(index3)];% torque
    
    Trim4.S_P1=[Trim4_res.f06data.Bendingmoment.USPLN1(index1),Trim4_res.f06data.Bendingmoment.LSPLN1(index2),Trim4_res.f06data.Bendingmoment.USPLN1(index3)]; % in plane shear
    Trim4.S_P2=[Trim4_res.f06data.Bendingmoment.USPLN2(index1),Trim4_res.f06data.Bendingmoment.LSPLN2(index2), Trim4_res.f06data.Bendingmoment.USPLN2(index3)]; % out of plane shear
    
    % forces from trim 5: g 20000ft
    Trim5.M_P1=[Trim5_res.f06data.Bendingmoment.UMPLN1(index1),Trim5_res.f06data.Bendingmoment.LMPLN1(index2),Trim5_res.f06data.Bendingmoment.UMPLN1(index3)]; % in-plane moment
    Trim5.M_P2=[Trim5_res.f06data.Bendingmoment.UMPLN2(index1),Trim5_res.f06data.Bendingmoment.LMPLN2(index2),Trim5_res.f06data.Bendingmoment.UMPLN2(index3)]; % out of plane moment
    Trim5.T=[Trim5_res.f06data.Bendingmoment.UTORQUE1(index1),Trim5_res.f06data.Bendingmoment.LTORQUE1(index2),Trim5_res.f06data.Bendingmoment.UTORQUE1(index3)];% torque
    
    Trim5.S_P1=[Trim5_res.f06data.Bendingmoment.USPLN1(index1),Trim5_res.f06data.Bendingmoment.LSPLN1(index2),Trim5_res.f06data.Bendingmoment.USPLN1(index3)]; % in plane shear
    Trim5.S_P2=[Trim5_res.f06data.Bendingmoment.USPLN2(index1),Trim5_res.f06data.Bendingmoment.LSPLN2(index2),Trim5_res.f06data.Bendingmoment.USPLN2(index3)]; % out of plane shear
    
    % forces from trim 6: g 3000ft
    Trim6.M_P1=[Trim6_res.f06data.Bendingmoment.UMPLN1(index1),Trim6_res.f06data.Bendingmoment.LMPLN1(index2), Trim6_res.f06data.Bendingmoment.UMPLN1(index3)]; % in-plane moment
    Trim6.M_P2=[Trim6_res.f06data.Bendingmoment.UMPLN2(index1),Trim6_res.f06data.Bendingmoment.LMPLN2(index2), Trim6_res.f06data.Bendingmoment.UMPLN2(index3)]; % out of plane moment
    Trim6.T=[Trim6_res.f06data.Bendingmoment.UTORQUE1(index1),Trim6_res.f06data.Bendingmoment.LTORQUE1(index2), Trim6_res.f06data.Bendingmoment.UMPLN2(index3)];% torque
    
    Trim6.S_P1=[Trim6_res.f06data.Bendingmoment.USPLN1(index1),Trim6_res.f06data.Bendingmoment.LSPLN1(index2),Trim6_res.f06data.Bendingmoment.USPLN1(index3)]; % in plane shear
    Trim6.S_P2=[Trim6_res.f06data.Bendingmoment.USPLN2(index1),Trim6_res.f06data.Bendingmoment.LSPLN2(index2),Trim6_res.f06data.Bendingmoment.USPLN2(index3)]; % out of plane shear
       
    % Gust responses
    WingNodes=NastranMethods1.WingNode;
    
    WingNodes_All=[WingNodes(1:end-1),FWTNodes];
    
    NumSteps=201;
    X=[WingNodes.X];
    X_FWT=[FWTNodes.X];
    
    Y_Data=X(2,:);
    Y_Data_FWT=X_FWT(2,:);
    Y_all=[Y_Data(1:end-1),Y_Data_FWT];

    [Max_Moment_LC1,Min_Moment_LC1, Max_Torque_LC1, Min_Torque_LC1, Max_Shear_LC1, Min_Shear_LC1]=Gust_peaks(WingNodes_All,GustLoadcase1,run_folder,'\gust_analysis_g_36000ft_pos_1MC.h5',NumSteps);
    [Max_Moment_LC2,Min_Moment_LC2, Max_Torque_LC2, Min_Torque_LC2, Max_Shear_LC2, Min_Shear_LC2]=Gust_peaks(WingNodes_All,GustLoadcase2,run_folder,'\gust_analysis_g_20000ft_pos_1MC.h5',NumSteps);
    [Max_Moment_LC3,Min_Moment_LC3, Max_Torque_LC3, Min_Torque_LC3, Max_Shear_LC3, Min_Shear_LC3]=Gust_peaks(WingNodes_All,GustLoadcase3,run_folder,'\gust_analysis_g_3000ft_pos_1MC.h5',NumSteps);
    
    
    %% Internal Load Summery 
    
    % Load case 1: 2.5 g pull up crusing altutude 
    Moment_P1_Loadcase1=abs(Trim1.M_P1); % internal load M1 for 2.5 g pull up
    Moment_P2_Loadcase1=abs(Trim1.M_P2); % internal load M2 for 2.5 g pull up
    Torque_Loadcase1=abs(Trim1.T);       % internal load T for 2.5 g pull up
    Shear_P1_Loadcase1=abs(Trim1.S_P1);
    Shear_P2_Loadcase1=abs(Trim1.S_P2);
    
    % Load case 2: 2.5 g pull up sea level 
    Moment_P1_Loadcase2=abs(Trim2.M_P1); % internal load M1 for 2.5 g pull up
    Moment_P2_Loadcase2=abs(Trim2.M_P2); % internal load M2 for 2.5 g pull up
    Torque_Loadcase2=abs(Trim2.T);       % internal load T for 2.5 g pull up
    Shear_P1_Loadcase2=abs(Trim2.S_P1);
    Shear_P2_Loadcase2=abs(Trim2.S_P2);
    
    % Load case 3: - g dive at sea level
    Moment_P1_Loadcase3=abs(Trim3.M_P1); % internal load M1 for 2.5 g pull up
    Moment_P2_Loadcase3=abs(Trim3.M_P2); % internal load M2 for 2.5 g pull up
    Torque_Loadcase3=abs(Trim3.T);       % internal load T for 2.5 g pull up
    Shear_P1_Loadcase3=abs(Trim3.S_P1);
    Shear_P2_Loadcase3=abs(Trim3.S_P2);
    
    % Load case 4: 36000ft g and + 1MC
    Moment_P1_Loadcase4=abs(Trim4.M_P1); 
    Moment_P2_Loadcase4=abs(Trim4.M_P2) + [Max_Moment_LC1; 0;]'; 
    Torque_Loadcase4=abs(Trim4.T) + [Max_Torque_LC1; 0]';       
    Shear_P1_Loadcase4=abs(Trim4.S_P1);
    Shear_P2_Loadcase4=abs(Trim4.S_P2) + [Max_Shear_LC1; 0]';
    
    % Load case 5: 20000ft g and + 1MC
    Moment_P1_Loadcase5=abs(Trim5.M_P1);
    Moment_P2_Loadcase5=abs(Trim5.M_P2) + [Max_Moment_LC2; 0]';
    Torque_Loadcase5=abs(Trim5.T) + [Max_Torque_LC2; 0]';      
    Shear_P1_Loadcase5=abs(Trim5.S_P1);
    Shear_P2_Loadcase5=abs(Trim5.S_P2) + [Max_Shear_LC2; 0]';
    
    % Load case 6: 3000ft g and + 1MC
    Moment_P1_Loadcase6=abs(Trim6.M_P1);
    Moment_P2_Loadcase6=abs(Trim6.M_P2) + [Max_Moment_LC3; 0]';
    Torque_Loadcase6=abs(Trim6.T) + [Max_Torque_LC3; 0]';
    Shear_P1_Loadcase6=abs(Trim6.S_P1);
    Shear_P2_Loadcase6=abs(Trim6.S_P2) + [Max_Shear_LC3; 0]';
    
    % Load case 7: 36000ft g and - 1MC
    Moment_P1_Loadcase7=abs(Trim4.M_P1);
    Moment_P2_Loadcase7=abs(Trim4.M_P2) + abs([Min_Moment_LC1; 0]');
    Torque_Loadcase7=abs(Trim4.T) + abs([Min_Torque_LC1; 0]');
    Shear_P1_Loadcase7=abs(Trim4.S_P1);
    Shear_P2_Loadcase7=abs(Trim4.S_P2) + abs([Min_Shear_LC1; 0]');
    
    % Load case 8: 20000ft g and - 1MC
    Moment_P1_Loadcase8=abs(Trim5.M_P1);
    Moment_P2_Loadcase8=abs(Trim5.M_P2) + abs([Min_Moment_LC2; 0]');
    Torque_Loadcase8=abs(Trim5.T) + abs([Min_Torque_LC2; 0]');
    Shear_P1_Loadcase8=abs(Trim5.S_P1);
    Shear_P2_Loadcase8=abs(Trim5.S_P2) +  abs([Min_Shear_LC2; 0]');
    
    % Load case 9: 3000ft g and - 1MC
    Moment_P1_Loadcase9=abs(Trim6.M_P1);
    Moment_P2_Loadcase9=abs(Trim6.M_P2) + abs([Min_Moment_LC3; 0]');
    Torque_Loadcase9=abs(Trim6.T) + abs([Min_Torque_LC3; 0]');
    Shear_P1_Loadcase9=abs(Trim6.S_P1);
    Shear_P2_Loadcase9=abs(Trim6.S_P2) + abs([Min_Shear_LC3; 0]');
    
    % Maximum forces used for sizing
    M_P1=max([Moment_P1_Loadcase1; Moment_P1_Loadcase2; Moment_P1_Loadcase3; Moment_P1_Loadcase4; Moment_P1_Loadcase5; Moment_P1_Loadcase6; Moment_P1_Loadcase7; Moment_P1_Loadcase8; Moment_P1_Loadcase9]);
    M_P2=max([Moment_P2_Loadcase1; Moment_P2_Loadcase2; Moment_P2_Loadcase3; Moment_P2_Loadcase4; Moment_P2_Loadcase5; Moment_P2_Loadcase6; Moment_P2_Loadcase7; Moment_P2_Loadcase8; Moment_P2_Loadcase9]);
    T=max([Torque_Loadcase1; Torque_Loadcase2; Torque_Loadcase3; Torque_Loadcase4; Torque_Loadcase5; Torque_Loadcase6; Torque_Loadcase7; Torque_Loadcase8; Torque_Loadcase9]);
    
    S_P1=max([Shear_P1_Loadcase1; Shear_P1_Loadcase2; Shear_P1_Loadcase3; Shear_P1_Loadcase4; Shear_P1_Loadcase5; Shear_P1_Loadcase6; Shear_P1_Loadcase7; Shear_P1_Loadcase8; Shear_P1_Loadcase9]);
    S_P2=max([Shear_P2_Loadcase1; Shear_P2_Loadcase2; Shear_P2_Loadcase3; Shear_P2_Loadcase4; Shear_P2_Loadcase5; Shear_P2_Loadcase6; Shear_P2_Loadcase7; Shear_P2_Loadcase8; Shear_P2_Loadcase9]);
    
    %% Result output
    
    % Out of plane bending moment distributions for all load cases
    Load_distribution.Moment_P2.LC1=Moment_P2_Loadcase1;
    Load_distribution.Moment_P2.LC2=Moment_P2_Loadcase2;
    Load_distribution.Moment_P2.LC3=Moment_P2_Loadcase3;
    Load_distribution.Moment_P2.LC4=Moment_P2_Loadcase4;
    Load_distribution.Moment_P2.LC5=Moment_P2_Loadcase5;
    Load_distribution.Moment_P2.LC6=Moment_P2_Loadcase6;
    Load_distribution.Moment_P2.LC7=Moment_P2_Loadcase7;
    Load_distribution.Moment_P2.LC8=Moment_P2_Loadcase8;
    Load_distribution.Moment_P2.LC9=Moment_P2_Loadcase9;
    
    % Vertical shear force distributions for all load cases
    Load_distribution.Shear_P2.LC1=Shear_P2_Loadcase1;
    Load_distribution.Shear_P2.LC2=Shear_P2_Loadcase2;
    Load_distribution.Shear_P2.LC3=Shear_P2_Loadcase3;
    Load_distribution.Shear_P2.LC4=Shear_P2_Loadcase4;
    Load_distribution.Shear_P2.LC5=Shear_P2_Loadcase5;
    Load_distribution.Shear_P2.LC6=Shear_P2_Loadcase6;
    Load_distribution.Shear_P2.LC7=Shear_P2_Loadcase7;
    Load_distribution.Shear_P2.LC8=Shear_P2_Loadcase8;
    Load_distribution.Shear_P2.LC9=Shear_P2_Loadcase9;
    
    % Torque distribution for all load cases 
    Load_distribution.Torque.LC1=Torque_Loadcase1;
    Load_distribution.Torque.LC2=Torque_Loadcase2;
    Load_distribution.Torque.LC3=Torque_Loadcase3;
    Load_distribution.Torque.LC4=Torque_Loadcase4;
    Load_distribution.Torque.LC5=Torque_Loadcase5;
    Load_distribution.Torque.LC6=Torque_Loadcase6;
    Load_distribution.Torque.LC7=Torque_Loadcase7;
    Load_distribution.Torque.LC8=Torque_Loadcase8;
    Load_distribution.Torque.LC9=Torque_Loadcase9;
    
    % Delta out of plane moment
    Delta.Moment_P2_Max.LC1=Max_Moment_LC1;
    Delta.Moment_P2_Max.LC2=Max_Moment_LC2;
    Delta.Moment_P2_Max.LC3=Max_Moment_LC3;
    
    Delta.Moment_P2_Min.LC1=Min_Moment_LC1;
    Delta.Moment_P2_Min.LC2=Min_Moment_LC2;
    Delta.Moment_P2_Min.LC3=Min_Moment_LC3;
    
    % Delta vertical shear
    Delta.Shear_Max.LC1=Max_Shear_LC1;
    Delta.Shear_Max.LC2=Max_Shear_LC2;
    Delta.Shear_Max.LC3=Max_Shear_LC3;
    
    Delta.Shear_Min.LC1=Min_Shear_LC1;
    Delta.Shear_Min.LC2=Min_Shear_LC2;
    Delta.Shear_Min.LC3=Min_Shear_LC3;
    
    % Delta torque
    Delta.Torque_Max.LC1=Max_Torque_LC1;
    Delta.Torque_Max.LC2=Max_Torque_LC2;
    Delta.Torque_Max.LC3=Max_Torque_LC3;
    
    Delta.Torque_Min.LC1=Min_Torque_LC1;
    Delta.Torque_Min.LC2=Min_Torque_LC2;
    Delta.Torque_Min.LC3=Min_Torque_LC3;
    
    
    

    %% calculate stresses
    
    E=70e9;

    sigma_skn=ones(1,NumSec);
    tau_skin=ones(1,NumSec);

    sigma_spr=ones(1,NumSec); 
    tau_spr=ones(1,NumSec);

    Von_skn=ones(1,NumSec);
    Von_spr=ones(1,NumSec);

    sigmab_skn=ones(1,NumSec);
    taub_skn=ones(1,NumSec);

    sigma_strg=ones(1,NumSec);
    sigma_crip=ones(1,NumSec);
    sigma_col=ones(1,NumSec);
    sigma_pp=ones(1,NumSec);

    for jj=1:NumSec

        % geometry info.
        ts=t_strg(jj);
        ds=d_strg(jj);

        hs=Bheight(jj);
        ws=Bwidth(jj);
        Ixx=Ixx_val(jj);
        Izz=Izz_val(jj);
        t1=thickness1(jj); % spar
        t2=thickness2(jj); %skin

        % load info
        m1=abs(M_P1(jj));
        m2=abs(M_P2(jj));
        tc=abs(T(jj));
        sp1=abs(S_P1(jj));
        sp2=abs(S_P2(jj));

        % Calc.
        sigma_skn(jj)=0.5*m2*hs/Ixx;
        tau_skin(jj)= tc/(2*hs*ws*t2)+ sp2/(2*hs*t1+2*ws*t2);
%         tau_skin(jj)=0;

        sigma_spr(jj)=0.5*m1*ws/Izz;
        % old shear
%         tau_spr(jj)=sp2/(2*hs*t1) + abs(tc)/(2*hs*ws*t1);
        Q=ws*t2*hs/2 + hs^2*t1/4;
       
        V=sp2/(2*hs*t1+2*ws*t2);
       
        tau_spr(jj)=V*Q/(2*Ixx*t1) + tc/(2*hs*ws*t1);

        % Von mises: skin and spar
        Von_skn(jj)=sqrt(sigma_skn(jj)^2);
        Von_spr(jj)=sqrt(sigma_spr(jj)^2+3*tau_spr(jj)^2);

        % stringers compressive stress

        sigma_strg(jj)=m2*(0.5*hs-0.5*ds)/Ixx;


        % critical buckling stress - for skin--------------------------
        sigmab_skn(jj)=4*pi^2*E*(t2/strg_n)^2/(12*(1-0.33^2));
        taub_skn(jj)=5.6*pi^2*E*(t2/strg_n)^2/(12*(1-0.33^2));      
        sigma_pp(jj) = 0.5*sigma_skn(jj)+sqrt((0.5*sigma_skn(jj))^2+tau_skin(jj)^2);

        % yield stress - skin and spar tensile failure
        sigma_Y=5.2E8;

        % striner crippling stress--------------------------------------
        A1=-0.7885; B1=0.6194; A2=-0.8046; B2=1.2117;

        Rc1_=B1*(sqrt(sigma_Y/E)*(ds/ts))^A1;
        Rc2_=B2*(sqrt(sigma_Y/E)*(ds/ts))^A2;

        sigmac1_=min(Rc1_,1.45)*sigma_Y;
        sigmac2_=min(Rc2_,1.45)*sigma_Y;

        sigma_crip(jj)=(ds*ts*sigmac1_*2 + ds*ts*sigmac2_)/(ts*ds*3);
        %-------------------------------------------------------------

        % column buckling

%         w_eff=1.7*t2*sqrt(E/abs(sigma_strg(jj)));
%         A_eff=w_eff*t2;
        L=0.65; %take as the rib pitch 0.65 m 
        A=3*ds*ts;

        I=(ts*ds^3/12)+(ds*ts^3/12 + ds*ts*(ds/2)^2)*2;
        rg=sqrt(I/A);
%         rg=sqrt(I/(ds*ts*3+A_eff));
%         K=L/(rg*1);
        
        K=L/(rg*2);

        Kcr=pi*sqrt(2*E/sigma_crip(1));

        if K < Kcr      
            sigma_col(jj)=sigma_crip(jj)*(1-sigma_crip(jj)*K^2/(4*pi^2*E));       
        elseif K > Kcr
            sigma_col(jj)=pi^2*E/K^2;
            
        end
    end
   
   
end
    




























