%% SOL 144
function [Von_skn, Von_spr, sigmab_skn, taub_skn, sigma_pp, sigma_strg, sigma_crip, sigma_col]=BeamStress_calc_v2(x)

    
%     thickness1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
%         x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
%         x(21),x(22)];
% 
%     thickness2=[x(23),x(24),x(25),x(26),x(27),x(28),x(29),x(30),x(31),x(32)...
%         x(33),x(34),x(35),x(36),x(37),x(38),x(39),x(40),x(41),x(42)...
%         x(43),x(44)];
% 
%     Astrg=[x(45),x(46),x(47),x(48),x(49),x(50),x(51),x(52),x(53),x(54)...
%         x(55),x(56),x(57),x(58),x(59),x(60),x(61),x(62),x(63),x(64)...
%         x(65),x(66)];

%     num=numel(x)/3;
% 
%     thickness1=x(1 : num);
%     thickness2=x(num+1 : num*2);
%     Astrg=x(num*2+1 : num*3);
    
   
     %% Mass configurations

    Payload_max=25300;
    Fuel_max=25000;

    MTOW=93500;
    OWE=48500;
    MWE=44057;
    Wing_mass=8906;

    Engine_mass=7362/2;
    Pylon=1239/2;
    Fuselage_mass=8986;
    
    %% Right root connector

    Connector_right = awi.model.LiftingSurface;
    Connector_right.Name = 'Connector_Right';
    Connector_right.Origin=[20,0,0]; %15

    %Use the Leading/Trailing edge sweep to define the planform
    Connector_right.ActiveSet = 'sSet';

    %Tail wing dimensions
    Connector_right.SpanVector  = 'Y';
    Connector_right.Span        = 2;
    Connector_right.LESweep     = [0, 0];
    Connector_right.LESweep_eta = [0, 1];
    Connector_right.TESweep     = [0,  0];
    Connector_right.TESweep_eta = [0,  1];
    Connector_right.RootChord   = 6;

    %Make sure the beam is at the midchord
    all_eta           = Connector_right.Eta_;
    Connector_right.BeamLoc     = repmat(0.5, size(all_eta));
    Connector_right.BeamLoc_eta = all_eta;
%     Connector_right.XOffset=35;
%     Tailwing_right.YOffset=1;

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
    E0  = 70e12; %[N/m^2], typical YM of aluminium
    nu0 = 0.333;
    rho0=2810;
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
    Connector_box_right.Height=0.2;
    Connector_box_right.Width=2.5;
    Connector_box_right.CoverThickness=0.001;
    Connector_box_right.SparThickness=0.002;
    getGeometricProps(Connector_box_right)
    
    Connector_right.BoxBeam = Connector_box_right;
    Connector_right.A   = Connector_box_right.Abb;
    Connector_right.I11 = Connector_box_right.Ixx;
    Connector_right.I22 = Connector_box_right.Izz;
    Connector_right.J   = Connector_box_right.Jbb;
    Connector_right.NSM = Connector_box_right.NSM;
    Connector_right.NSI = Connector_box_right.NSI;

    for i=1:1:3
        handle_connectorR=strcat('PM_tail_R','i');
        handle_connectorR=awi.model.PointMass;
        handle_connectorR.SOffset=-0.1+i*0.2;
        handle_connectorR.Mass=10;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle_connectorR.MassGroup='Group3';
        Connector_right.add(handle_connectorR);

    end
    
    % Aeropanel definition
    Connector_right.AeroPanelLength=0.5;

    build(Connector_right);
    
    
    %% Left root connector
    
    Connector_left = awi.model.LiftingSurface;
    Connector_left.Name = 'Connector_Right';
    Connector_left.Origin=[20,0,0]; %15
    
    %Use the Leading/Trailing edge sweep to define the planform
    Connector_left.ActiveSet = 'sSet';

    %Tail wing dimensions
    Connector_left.SpanVector  = 'Y';
    Connector_left.Span        = -2;
    Connector_left.LESweep     = [0, 0];
    Connector_left.LESweep_eta = [0, 1];
    Connector_left.TESweep     = [0,  0];
    Connector_left.TESweep_eta = [0,  1];
    Connector_left.RootChord   = 6;

    %Make sure the beam is at the midchord
    all_eta           = Connector_left.Eta_;
    Connector_left.BeamLoc     = repmat(0.5, size(all_eta));
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
    Connector_box_left.Height=0.2;
    Connector_box_left.Width=2.5;
    Connector_box_left.CoverThickness=0.001;
    Connector_box_left.SparThickness=0.002;
    getGeometricProps(Connector_box_left)
    
    Connector_left.BoxBeam = Connector_box_left;
    Connector_left.A   = Connector_box_left.Abb;
    Connector_left.I11 = Connector_box_left.Ixx;
    Connector_left.I22 = Connector_box_left.Izz;
    Connector_left.J   = Connector_box_left.Jbb;
    Connector_left.NSM = Connector_box_left.NSM;
    Connector_left.NSI = Connector_box_left.NSI;

    for i=1:1:3
        handle_connectorL=strcat('PM_tail_R','i');
        handle_connectorL=awi.model.PointMass;
        handle_connectorL.SOffset=-0.1+i*0.2;
        handle_connectorL.Mass=10;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle_connectorL.MassGroup='Group3';
        Connector_left.add(handle_connectorL);

    end
    
    % Aeropanel definition
    Connector_left.AeroPanelLength=0.5;

    build(Connector_left);
   
%     FEM_test=convertToFE(Connector_left);
%     draw(FEM_test)
    
  %% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[20,2,0]; %15
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';
    
    % Num of element 
    Wingbox_right.NumBeamElem = 23;

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = 16;   %34.1/2;
    Wingbox_right.LESweep     = [27, 27];
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [0, 16.59, 16.59];
    Wingbox_right.TESweep_eta = [0, 0.27, 1];
    Wingbox_right.RootChord   = 6;
    
    
    %Dihedral 
    Wingbox_right.Dihedral=[5,5];
    Wingbox_right.Dihedral_eta=[0,1];
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_right.Eta_;
    Wingbox_right.BeamLoc     = repmat(0.5, size(all_eta));

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
    
    
    %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=Wingbox_right.NumBeamElem+2;
    
    %% sizing variables ---------------------------------------------
%     x=[0.02*ones(1,NumSec),0.005*ones(1,NumSec),1e-7*ones(1,NumSec)];
    
    thickness1=x(1:NumSec);
    thickness2=x(NumSec+1:NumSec*2);
    Astrg=x(NumSec*2+1:NumSec*3);
        
    d_strg=sqrt(Astrg/0.36);
    t_strg=0.12*d_strg;
    %-----------------------------------------------------------------

    % etaS=linspace(0,Wingbox_right.Span,NumSec);

    % set width and height array 
    YData=Wingbox_right.YData;
    SparWidth=Wingbox_right.Chord*0.5;

    RootH=Wingbox_right.Chord(1)*0.15;
    MidH=Wingbox_right.Chord(2)*0.12;
    TipH=Wingbox_right.Chord(end)*0.11;


    % set up eta values
    elnum=Wingbox_right.NumBeamElem + 1;
    Num_seg1=ceil(elnum*0.27);
    Num_seg2=elnum - Num_seg1;
    
    Num_sec1=Num_seg1+1;
    Num_sec2=Num_seg2+1;
    
    eta1_=linspace(0,0.27, Num_sec1);
    eta2_=linspace(0.27,1,Num_sec2);
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
    Wingbox_right.NSM=NSM_val;
    Wingbox_right.NSM_eta= eta_;
    
    Wingbox_right.NSI=NSI_val;
    Wingbox_right.NSI_eta= eta_;
    
    
    
    % Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
%     Wingbox_right.NumAeroPanel=20;
    Wingbox_right.AeroPanelLength=0.4;
    
    build(Wingbox_right)
    
    
%     FEM_test=convertToFE(Wingbox_right);
%     draw(FEM_test)

    %% Mass definition
    
    % total wing mass
    
    [wing_mass,total_mass]=Mass_calc_v2(x);
    
    wingmass_eta=0.04:0.04:1;
%     wingmass_eta=0.2:0.2:1;
    Mwidth=interp1(RData/RData(end),SparWidth,wingmass_eta);
    
    mass_set=(total_mass-wing_mass)*(Mwidth)/sum(Mwidth);
      
%     m=total_mass/19;
    
    for i=1:1:25
        handle=strcat('PM_right','i');
        handle=awi.model.PointMass;
        handle.SOffset=0+i*0.04;
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
    
    %Engine location 
    Y_eng=Wingbox_right.PanelCoords.LE.Y(2);
    X_eng=(Wingbox_right.PanelCoords.LE.X(2)+Wingbox_right.PanelCoords.TE.X(2))/2;
    Z_eng=(Wingbox_right.PanelCoords.LE.Z(2)+Wingbox_right.PanelCoords.TE.Z(2))/2;
    
    Engine.Origin = [X_eng-3.5, Y_eng, Z_eng];
  
    Engine.Length = 3.5;
    
    %Make the material
    E1  = 76e9; %[N/m^2],set as a rigid body
    nu = 0.333;
    Mat1 = awi.model.Material;
    Mat1.E  = E1;
    Mat1.Nu = nu;
    Mat1.G  = E1 / (2 * (1 + nu));
    Mat1.Rho=2800;
    
    
    % use the strong material
    Engine.Material_eta = [0, 1];
    Engine.Material     = [Mat1, Mat1];

    Engine.A   = 0.04432;
    Engine.I11 = 0.002;
    Engine.I22 = 0.002;
    Engine.J   = 0.001636;
    
    %Aeropanel althoufh it is useless now
    Engine.AeroPanelLength=0.5;
    
    % add engine mass
    engine_mass=awi.model.PointMass;   
    engine_mass.SOffset=0.0;
    engine_mass.Mass=Engine_mass;
    Engine.add(engine_mass);
    
    % add pylon
    pylon_mass=awi.model.PointMass;   
    pylon_mass.SOffset=0.9;
    pylon_mass.Mass=Pylon;
    Engine.add(pylon_mass);

    build(Engine)
      
    Wingbox_right.add(Engine)
    
    
    
    %Control surfaces - flaps
    flap_R=awi.model.ControlSurface;
    flap_R.Eta=[0, 0.27];
    flap_R.xLE=[0.8,0.8];
    flap_R.xTE=[1,1];
    flap_R.Max_def=0.1;
    flap_R.Max_rate=0.1;
    flap_R.HingeLine='LE';
    flap_R.Label='FlapR';
    flap_R.FaceColor='m';
    
%     flap_R.NumAeroPanel=10;
    flap_R.AeroPanelLength=0.4;
    
    build(flap_R)
    Wingbox_right.add(flap_R);
    
    Wingbox_right.ModelControlSurf = 1;
    
    
    build(Wingbox_right);
    
%     FEM_test=convertToFE(Wingbox_right);
%     
%     draw(FEM_test)


%% Wingbox 2 - left and control surf.

    Wingbox_left = awi.model.LiftingSurface;
    Wingbox_left.Name = 'A320Wing_left';
    Wingbox_left.Origin=[20,-2,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_left.ActiveSet = 'sSet';

    %Wing dimensions
    Wingbox_left.SpanVector  = 'Y';
    Wingbox_left.Span        = -16;  
    Wingbox_left.LESweep     = [-27, -27];
    Wingbox_left.LESweep_eta = [0, 1];
    Wingbox_left.TESweep     = [0, -16.59, -16.59];
    Wingbox_left.TESweep_eta = [0, 0.27, 1];
    Wingbox_left.RootChord   = 6;   
    
    %Dihedral 
    Wingbox_left.Dihedral=[5,5];
    Wingbox_left.Dihedral_eta=[0,1];
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_left.Eta_;
    Wingbox_left.BeamLoc     = repmat(0.5, size(all_eta));
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

    %Make the material
    E  = 70e9; %[N/m^2], typical YM of aluminium
    nu = 0.333;
    rho=2810;
    Mat = awi.model.Material;
    Mat.E  = E;
    Mat.Nu = nu;
    Mat.G  = E / (2 * (1 + nu));
    Mat.Rho=rho;
    Wingbox_left.Material_eta = [0, 1];
    Wingbox_left.Material     = [Mat, Mat];
    
    build(Wingbox_left)
    
    
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

    % NSM and NSI
    Wingbox_left.NSM=NSM_val;
    Wingbox_left.NSM_eta= eta_;
    
    Wingbox_left.NSI=NSI_val;
    Wingbox_left.NSI_eta= eta_;
    
       
    % Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
%     Wingbox_left.NumAeroPanel=20;
    Wingbox_left.AeroPanelLength=0.4;
    
    build(Wingbox_left)
   

    %% Mass definition
    
    % total wing mass
    
%     [wing_mass,total_mass]=Mass_calc_v1(x);
%     
%     wingmass_eta=0.04:0.04:1;
% %     wingmass_eta=0.2:0.2:1;
%     Mwidth=interp1(RData/RData(end),SparWidth,wingmass_eta);
%     
%     mass_set=(total_mass-wing_mass)*(Mwidth)/sum(Mwidth);
      
%     m=total_mass/19;
    
    for i=1:1:25
        handle=strcat('PM_left','i');
        handle=awi.model.PointMass;
        handle.SOffset=0+i*0.04;
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
    
    %Engine location 
    Y_eng=Wingbox_left.PanelCoords.LE.Y(2);
    X_eng=(Wingbox_left.PanelCoords.LE.X(2)+Wingbox_left.PanelCoords.TE.X(2))/2;
    Z_eng=(Wingbox_left.PanelCoords.LE.Z(2)+Wingbox_left.PanelCoords.TE.Z(2))/2;
    
    Engine2.Origin = [X_eng-3.5, Y_eng, Z_eng];
    
  
    Engine2.Length = 3.5;
%     Engine.XOffset=16.471008-3.5;
%     Engine.YOffset=5.8170588;
%     Engine.ZOffset=-2;
    
    %Make the material
    E1  = 76e9; %[N/m^2],set as a rigid body
    nu = 0.333;
    Mat1 = awi.model.Material;
    Mat1.E  = E1;
    Mat1.Nu = nu;
    Mat1.G  = E1 / (2 * (1 + nu));
    Mat1.Rho=2800;
    
    
    % use the strong material
    Engine2.Material_eta = [0, 1];
    Engine2.Material     = [Mat1, Mat1];

    Engine2.A   = 0.04432;
    Engine2.I11 = 0.002;
    Engine2.I22 = 0.002;
    Engine2.J   = 0.001636;
    
    %Aeropanel althoufh it is useless now
    Engine2.AeroPanelLength=0.5;
    
    % add engine mass
    engine2_mass=awi.model.PointMass;   
    engine2_mass.SOffset=0.0;
    engine2_mass.Mass=Engine_mass;
    Engine2.add(engine2_mass);
    
    % add pylon
    pylon2_mass=awi.model.PointMass;   
    pylon2_mass.SOffset=0.9;
    pylon2_mass.Mass=Pylon;
    Engine2.add(pylon2_mass);

    build(Engine2)
      
    Wingbox_left.add(Engine2)
    
    %Control surfaces - flaps
    flap_L=awi.model.ControlSurface;
    flap_L.Eta=[0, 0.25];
    flap_L.xLE=[0.8,0.8];
    flap_L.xTE=[1,1];
    flap_L.Max_def=0.1;
    flap_L.Max_rate=0.1;
    flap_L.HingeLine='LE';
    flap_L.Label='FlapL';
    flap_L.FaceColor='m';
    
    flap_L.AeroPanelLength=0.4;
    
    build(flap_R)
    Wingbox_left.add(flap_L);
    
    Wingbox_left.ModelControlSurf = 1;
     
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
    Body.Length=45;
    % Body.XOffset=-15;

    
    %Make the material
    E1  = 76e19; %[N/m^2],set as a rigid body
    nu = 0.333;
    Mat1 = awi.model.Material;
    Mat1.E  = E1;
    Mat1.Nu = nu;
    Mat1.G  = E1 / (2 * (1 + nu));

    
    % use the strong material
    Body.Material=Mat1;
    Body.Material_eta = [0, 1];
    Body.Material     = [Mat1, Mat1];

    %define  panel size
    % Body.NumAeroPanel=5;
    Body.AeroPanelLength=0.5;

    % define cross sectional definition
    Bodybox=awi.model.BoxBeam;
    Bodybox.BoxType='SymmetricBox';
    Bodybox.Height=2;
    Bodybox.Width=2;
    Bodybox.CoverThickness=0.01;
    Bodybox.SparThickness=0.01;
    getGeometricProps(Bodybox)
    Body.BoxBeam = Bodybox;
    Body.A   = Bodybox.Abb;
    Body.I11 = Bodybox.Ixx;
    Body.I22 = Bodybox.Izz;
    Body.J   = Bodybox.Jbb;
%     Body.NSM = Bodybox.NSM;
    Body.NSI = Bodybox.NSI;

    % add point masses- fuselage 
    M_fs=MWE-Wing_mass;
    m_fs=M_fs/11;
    
    %payload
    M_p=Payload_max;
    m_p= M_p/11;
    
    for i=1:1:11
        handle=strcat('PM_body','i');
        handle=awi.model.PointMass;
        handle.SOffset=-0.1+i*0.1;
        handle.Mass=m_fs+m_p;
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
    Tailwing_right.XOffset = 42;
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

    % material properties
    Tailwing_right.Material_eta = [0, 1];
    Tailwing_right.Material     = [Mat1, Mat1];

    % Define box beam corss section
    tailbox_right=awi.model.BoxBeam;
    tailbox_right.BoxType='SymmetricBox';
    tailbox_right.Height=0.2;
    tailbox_right.Width=2.5;
    tailbox_right.CoverThickness=0.001;
    tailbox_right.SparThickness=0.002;
    getGeometricProps(tailbox_right)
    Tailwing_right.BoxBeam = tailbox_right;
    Tailwing_right.A   = tailbox_right.Abb;
    Tailwing_right.I11 = tailbox_right.Ixx;
    Tailwing_right.I22 = tailbox_right.Izz;
    Tailwing_right.J   = tailbox_right.Jbb;
%     Tailwing_right.NSM = tailbox_right.NSM;
%     Tailwing_right.NSI = tailbox_right.NSI;

    for i=1:1:3
        handle_tailR=strcat('PM_tail_R','i');
        handle_tailR=awi.model.PointMass;
        handle_tailR.SOffset=-0.1+i*0.2;
        handle_tailR.Mass=10;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle.MassGroup='Group3';
        Tailwing_right.add(handle_tailR);

    end
    

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
    Tailwing_left.XOffset=42;
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
    Tailwing_left.Material     = [Mat1, Mat1];

    % Define box beam corss section
    tailbox_left=awi.model.BoxBeam;
    tailbox_left.BoxType='SymmetricBox';
    tailbox_left.Height=0.2;
    tailbox_left.Width=2.5;
    tailbox_left.CoverThickness=0.001;
    tailbox_left.SparThickness=0.002;
    getGeometricProps(tailbox_left)
    Tailwing_left.BoxBeam = tailbox_left;
    Tailwing_left.A   = tailbox_left.Abb;
    Tailwing_left.I11 = tailbox_left.Ixx;
    Tailwing_left.I22 = tailbox_left.Izz;
    Tailwing_left.J   = tailbox_left.Jbb;
%     Tailwing_right.NSM = tailbox_right.NSM;
%     Tailwing_right.NSI = tailbox_right.NSI;

    for i=1:1:3
        handle_tailL=strcat('PM_tail_L','i');
        handle_tailL=awi.model.PointMass;
        handle_tailL.SOffset=-0.1+i*0.2;
        handle_tailL.Mass=10;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle_tailL.MassGroup='Group4';
        Tailwing_left.add(handle_tailL);

    end

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
    Verticalwing.XOffset=41;


    Verticalwing.Material_eta = [0, 1];
    Verticalwing.Material     = [Mat1, Mat1];

    % % Aeropanel definition
    Verticalwing.NumAeroPanel=8;

    % Define box beam corss section
    Verticalbox=awi.model.BoxBeam;
    Verticalbox.BoxType='SymmetricBox';
    Verticalbox.Height=0.2;
    Verticalbox.Width=2.5;
    Verticalbox.CoverThickness=0.001;
    Verticalbox.SparThickness=0.002;
    getGeometricProps(Verticalbox)
    Verticalwing.BoxBeam = Verticalbox;
    Verticalwing.A   = Verticalbox.Abb;
    Verticalwing.I11 = Verticalbox.Ixx;
    Verticalwing.I22 = Verticalbox.Izz;
    Verticalwing.J   = Verticalbox.Jbb;
%     Verticalwing.NSM = Verticalbox.NSM;
%     Verticalwing.NSI = Verticalbox.NSI;

    for i=1:1:2
        handle_vertical=strcat('PM_vertical_L','i');
        handle_vertical=awi.model.PointMass;
        handle_vertical.SOffset=-0.1+i*0.5;
        handle_vertical.Mass=10;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle_vertical.MassGroup='Group5';
        Verticalwing.add(handle_tailL);

    end

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
        Connector_right.SurfaceArea,  Connector_left.SurfaceArea]);
    
    Aircraft.RefSpan  = Wingbox_right.Span;
    Aircraft.RefChord = Wingbox_right.RootChord*0.697; %*0.697;
%     Aircraft.RefChord = Aircraft.RefArea/Aircraft.RefSpan; 


    %% Generate the loadcase object
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
    
%     TrimLoadcase1.LoadCaseType = 'Pratt Gust';  
%     
%     n=calculateGustLoadFactor(TrimLoadcase1,Aircraft);

    TrimLoadcase1.PitchAngle=0;
    TrimLoadcase1.RollAngle =0;
    TrimLoadcase1.ID = 1020;
    TrimLoadcase1.LoadFactor = 2.5;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase1.CsDeflection=flap_angle*pi/180;
    
    
    build(TrimLoadcase1)
    
    
    %% Trim loadcase 2
    TrimLoadcase2 = awi.model.LoadCase;

    acMass = 500;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;

    TrimLoadcase.Name = 'A321_cruise_1g';
    % TrimLoadcase.LoadCaseTypes = 'Static';(read only)
    % TrimLoadcase.CsDeflecTypes='fixed';(read only)
    TrimLoadcase2.Altitude   = altitude;
    TrimLoadcase2.Mach       = mach_number;
    TrimLoadcase2.AcVelocity = aircraft_velocity;
    TrimLoadcase2.AcMass = acMass;

    TrimLoadcase2.PitchAngle=0;
    TrimLoadcase2.RollAngle =0;
    TrimLoadcase2.ID = 1030;
    TrimLoadcase2.LoadFactor = -1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase2.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase2)

    %% Generate the FEM 

    run_folder = [
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\full_A321_sizing1']; %[-], folder for exporting the NASTRAN model


    % Convert to a finite element model and draw it
    FEM_half = convertToFE(Aircraft);

    % %Export it to a file
    export(FEM_half, run_folder);


  %% NASTRAN method - RUN SOL 144
  
%     % NastranExe: 'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe'
% 
%     NastranMethods1 = awi.methods.Nastran;
%     NastranMethods1.AnalysisModel = FEM_half;
%     MassCases=awi.model.MassCases.empty;
%     ID0=200;
%     % trimdata=NastranMethods1.getTrimData(FEM_half, Aircraft, TrimLoadcase, MassCases, ID0, NastranMethods1.RefNode);
    
    
     % Run the analysis- SOL 144 static trim analysis 
    
%      NastranExe: 'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe'
     
     NastranMethods1 = awi.methods.Nastran;
     NastranMethods1.AnalysisModel = FEM_half;
     MassCases=awi.model.MassCases.empty;
     ID0=200;    
     
    trimFile1 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase1, MassCases,run_folder,'DatFilename','A321_cruise_2p5g');
    
    trimFile2 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase2, MassCases,run_folder,'DatFilename','A321_cruise_1g');
    
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\A320_half_model_SOL144*.*')
   
    delete(strcat(run_folder, '\A321_cruise_2p5g*.xdb'));
    delete(strcat(run_folder, '\A321_cruise_2p5g*.h5'));
    delete(strcat(run_folder, '\A321_cruise_2p5g*.log'));
    delete(strcat(run_folder, '\A321_cruise_2p5g*.f06'));
    delete(strcat(run_folder, '\A321_cruise_2p5g*.f04'));
    
    delete(strcat(run_folder, '\A321_cruise_1g*.xdb.*'));
    delete(strcat(run_folder, '\A321_cruise_1g*.h5.*'));
    delete(strcat(run_folder, '\A321_cruise_1g*.log.*'));
    delete(strcat(run_folder, '\A321_cruise_1g*.f06.*'));
    delete(strcat(run_folder, '\A321_cruise_1g*.f04.*'));
      
    NastranMethods1.runNastran(trimFile1);
    NastranMethods1.runNastran(trimFile2);
    
    
    %% extract results

    Trim1_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_cruise_2p5g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim2_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_cruise_1g.f06'),'ReadF06',true,'ReadHDF5',false);

    % extract forces
    Node0=108;
    Node1=131;
    
    % forces from trim1
    Trim1.M_P1=[Trim1_res.f06data.Bendingmoment.UMPLN1(Node0:Node1),Trim1_res.f06data.Bendingmoment.LMPLN1(Node1)]; % in-plane moment
    Trim1.M_P2=[Trim1_res.f06data.Bendingmoment.UMPLN2(Node0:Node1),Trim1_res.f06data.Bendingmoment.LMPLN2(Node1)]; % out of plane moment
    Trim1.T=[Trim1_res.f06data.Bendingmoment.UTORQUE1(Node0:Node1),Trim1_res.f06data.Bendingmoment.LTORQUE1(Node1)];% torque

    Trim1.S_P1=[Trim1_res.f06data.Bendingmoment.USPLN1(Node0:Node1),Trim1_res.f06data.Bendingmoment.LSPLN1(Node1)]; % in plane shear
    Trim1.S_P2=[Trim1_res.f06data.Bendingmoment.USPLN2(Node0:Node1),Trim1_res.f06data.Bendingmoment.LSPLN2(Node1)]; % out of plane shear
    
    % forces from trim2
    Trim2.M_P1=[Trim2_res.f06data.Bendingmoment.UMPLN1(Node0:Node1),Trim2_res.f06data.Bendingmoment.LMPLN1(Node1)]; % in-plane moment
    Trim2.M_P2=[Trim2_res.f06data.Bendingmoment.UMPLN2(Node0:Node1),Trim2_res.f06data.Bendingmoment.LMPLN2(Node1)]; % out of plane moment
    Trim2.T=[Trim2_res.f06data.Bendingmoment.UTORQUE1(Node0:Node1),Trim2_res.f06data.Bendingmoment.LTORQUE1(Node1)];% torque

    Trim2.S_P1=[Trim2_res.f06data.Bendingmoment.USPLN1(Node0:Node1),Trim2_res.f06data.Bendingmoment.LSPLN1(Node1)]; % in plane shear
    Trim2.S_P2=[Trim2_res.f06data.Bendingmoment.USPLN2(Node0:Node1),Trim2_res.f06data.Bendingmoment.LSPLN2(Node1)]; % out of plane shear
       
    % Maximum forces used for sizing
    M_P1=max([abs(Trim1.M_P1); abs(Trim2.M_P1)]);
    M_P2=max([abs(Trim1.M_P2); abs(Trim2.M_P2)]);
    T=max([abs(Trim1.T); abs(Trim2.T)]);
    
    S_P1=max([abs(Trim1.S_P1); abs(Trim2.S_P1)]);
    S_P2=max([abs(Trim1.S_P2); abs(Trim2.S_P2)]); 
    
    
    % calculate stresses

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
        sigma_Y=5E8;

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
    




























