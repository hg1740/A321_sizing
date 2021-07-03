  
%% sizing parameters thickness1 = spar, thickness2 = skin   

% %     x=[0.0191135910000000,0.0193953730000000,0.0201137600000000,0.0213823650000000,0.0229999860000000,0.0255251240000000,0.0285624380000000,0.0372113470000000,0.0373207110000000,0.0372799900000000,0.0368848460000000,0.0378039760000000,0.0379687300000000,0.0379303110000000,0.0378464170000000,0.0368959070000000,0.0365928770000000,0.0358741000000000,0.0351157990000000,0.0325204110000000,0.0321366360000000,0.0317223020000000,0.0280889430000000,0.0226448800000000,0.0118821390000000,0.00707914700000000,0.00718361400000000,0.00730636500000000,0.00745189900000000,0.00763438400000000,0.00783787200000000,0.00810313600000000,0.00823100800000000,0.00803824900000000,0.00781619300000000,0.00760075000000000,0.00736034800000000,0.00709415700000000,0.00680991600000000,0.00650766900000000,0.00618441800000000,0.00597655500000000,0.00574729800000000,0.00545697300000000,0.00516124700000000,0.00473319700000000,0.00418339100000000,0.00343190400000000,0.00256919700000000,0.000896999000000000,6.20000000000000e-05,6.38000000000000e-05,6.58000000000000e-05,6.82000000000000e-05,7.14000000000000e-05,7.49000000000000e-05,7.99000000000000e-05,8.21000000000000e-05,7.78000000000000e-05,7.33000000000000e-05,6.92000000000000e-05,6.48000000000000e-05,6.03000000000000e-05,5.55000000000000e-05,5.06000000000000e-05,4.55000000000000e-05,4.24000000000000e-05,3.92000000000000e-05,3.52000000000000e-05,3.13000000000000e-05,2.63000000000000e-05,2.05000000000000e-05,1.40000000000000e-05,8.02000000000000e-06,9.15000000000000e-07];
% %     

    x=0.002*ones(1,75);

    %% Wing configurations for starboard wing
  
     Aspect_ratio=10.172; % Aspect ratio = 10.172 for A321 model
    
    Total_area=126;         % include two wing surface areas + floor size on the fuselage
    Fuselage_width=4;       % dimeter of the fuselage
  
    Wing_span = sqrt(Aspect_ratio*Total_area);
    BeamLoc = 0.4;          % choose a spar location: 0 --> 1
    Semi_span=(Wing_span-Fuselage_width)/2; % length of one wing: 16m for A321 model
    
    Root_chord =  Total_area/(1.064*Semi_span + 4);
    LE_sweep=27;            % deg
    
    Wing_area = (Total_area - Fuselage_width*Root_chord)/2;
    
    Mid_chord=0.63685*Root_chord;
    Tip_chord=0.2248*Root_chord;
    
    X0=Root_chord; 
    X1=0.27*Semi_span*tan(27*pi/180) + Mid_chord;
    X2=Semi_span*tan(27*pi/180) + Tip_chord;
    
    tan_TE_sweep1=(X1-X0)/(0.27*Semi_span);
    tan_TE_sweep2=(X2-X1)/(0.73*Semi_span);
    
    TE_sweep1=atan(tan_TE_sweep1)*180/pi; % deg
    TE_sweep2=atan(tan_TE_sweep2)*180/pi; % deg
      
    Taper_ratio=Tip_chord/Root_chord;
    
    Mean_cord_coefficient=(2/3)*(1+Taper_ratio+Taper_ratio^2)/(1+Taper_ratio);
    
    %% Mass configurations

    Payload_max=25000; % kg
    Fuel_fraction=0.723; % percentage of fuel in the tank

    Fuel_capacity=32940; % L

    MTOW=93500; % maximum take off mass
    OWE=48500;  % Operating empty mass
    MWE=44057;  % Manufacture's empty mass
    MZF=73000;  % Maxumum zero fuel mass

    Engine_mass=7362/2; % kg
    Pylon=1239/2; % kg
    Horizontal_tail=682; % kg
    Vertical_tail=522; % kg

    % For sized result, estimated wing mass = totol wing mass 2701.3 kg for
    % A 321 model
    
%     [Wingbox_mass, Secondary_mass, Wing_total_mass]=Wing_masses_v1(x,Height,Width,etaRL,MTOW,Surfacearea,LE_Sweep,Semi_span);
     Wing_total_mass=2747.3;
     Secondary_mass=835.2;

    % The tube shell is assumed to be made by 4mm thick Aluminium skin. 
    
    Fuselage_shell_mass=2*pi*2*0.004*2800*44.5; %kg
    
    Fuselage_structure_mass=(OWE - 2*Wing_total_mass - Horizontal_tail - Vertical_tail - Engine_mass*2 - Pylon*2 - Fuselage_shell_mass);
    
    % obtain lumped masses: fuel masses on the wing;
    % Fuselage_total_mass = the structure mass  +  payload 
    
    % the sizing is carroied out using max. payload with its corresponding
    % max fuel mass 
    
    [Fuel_mass,Fuselage_total_mass]= Lumped_masses(Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
    
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
    Connector_right.RootChord   = Root_chord;
    
    % Jig - twist
%     Connector_right.AoA=[3.5,3.5];
%     Connector_right.AoA_eta=[0,1];

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
    Connector_box_right.CoverThickness=0.008;
    Connector_box_right.SparThickness=0.008;
    getGeometricProps(Connector_box_right)
    
    Connector_right.BoxBeam = Connector_box_right;
    Connector_right.A   = Connector_box_right.Abb;
    Connector_right.I11 = Connector_box_right.Ixx;
    Connector_right.I22 = Connector_box_right.Izz;
    Connector_right.J   = Connector_box_right.Jbb;


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
    Connector_right.AeroPanelAR=2;

    build(Connector_right);
    
%     draw(Connector_right)
    
    
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
    Connector_left.RootChord   = Root_chord;
    
    % Jig - twist
%     Connector_left.AoA=[3.5,3.5];
%     Connector_left.AoA_eta=[0,1];

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
    Connector_box_left.CoverThickness=0.008;
    Connector_box_left.SparThickness=0.008;
    getGeometricProps(Connector_box_left)
    
    Connector_left.BoxBeam = Connector_box_left;
    Connector_left.A   = Connector_box_left.Abb; 
    Connector_left.I11 = Connector_box_left.Ixx;
    Connector_left.I22 = Connector_box_left.Izz;
    Connector_left.J   = Connector_box_left.Jbb;


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
    Connector_left.AeroPanelAR=2;
    
    build(Connector_left);
    

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
    rho_wing=2800; 
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
    Wingbox_right.NumAeroPanel=15;
    Wingbox_right.AeroPanelAR=2.5;
    %     Wingbox_right.AeroPanelLength=0.4;
    
    
    %% initialise properties
    
    
    Wingbox_right.A   =  1;
    Wingbox_right.A_eta=[0,1];
    
    Wingbox_right.I11 = 1;
    Wingbox_right.I11_eta=[0,1];
    
    Wingbox_right.I22 = 1;
    Wingbox_right.I22_eta = [0,1];
    
    Wingbox_right.J   = 1;
    Wingbox_right.J_eta= [0,1];
    
    %% FWT right definition
    
    fold_angle  = -10;   %[deg],
    flare_angle = 25;   %[deg],
    fold_eta=0.75;
    
    
    FWT_R = insertWingFold(Wingbox_right, 'FlareAngle', flare_angle, 'FoldAngle', fold_angle,'EtaFold',fold_eta);
    FWT_R.HingeStiffness = [1e14 1e14 1e14 1e14 1e-4 1e14];
    FWT_R.NumAeroPanel=15;
    
    FWT_R.NumBeamElem = 10;
    
    FWT_box=awi.model.BoxBeam;
    FWT_box.BoxType='SymmetricBox';
    FWT_box.Height=0.18;
    FWT_box.Width=0.18;
    FWT_box.CoverThickness=0.005;
    FWT_box.SparThickness=0.005;
    
    getGeometricProps(FWT_box)
    
    sc_num=numel(FWT_R.A);
    % FWT_eta=[0,1];
    
    FWT_R.A   =  ones(1,sc_num).*FWT_box.Abb;
    % FWT.A_eta=FWT_eta;
    
    FWT_R.I11 = ones(1,sc_num).*FWT_box.Izz;
    % FWT.I11_eta=FWT_eta;
    
    FWT_R.I22 = ones(1,sc_num).*FWT_box.Ixx;
    % FWT.I22_eta = FWT_eta;
    
    FWT_R.J   = ones(1,sc_num).*FWT_box.Jbb;
    % FWT.J_eta= FWT_eta;
    
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
       
%         NSM_val(ii)=boxname.NSM;
%         NSI_val(ii)=boxname.NSI;
        
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
    
    % total wing mass
    
%     [wing_mass,total_mass]=Mass_calc_v2(x);
    
%     wingmass_eta=0.04:0.04:1;
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
    
    %Engine location 
    Y_eng=Wingbox_right.PanelCoords.LE.Y(2);
    X_eng=(Wingbox_right.PanelCoords.LE.X(2)+Wingbox_right.PanelCoords.TE.X(2))/2 - (Wingbox_right.PanelCoords.TE.X(2)-Wingbox_right.PanelCoords.LE.X(2))*0.1;
    Z_eng=(Wingbox_right.PanelCoords.LE.Z(2)+Wingbox_right.PanelCoords.TE.Z(2))/2;
    
    Engine.Origin = [X_eng-3.5, Y_eng, Z_eng];
  
    Engine.Length = 3.5;
    
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
    Engine_Inertia = Engine_Mat.Rho*Engine_A*Engine_radius^2;
    
    Engine.A   = Engine_A*0.00001; % temp
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
    
    
    
    %Control surfaces - flaps
%     flap_R=awi.model.ControlSurface;
%     flap_R.Eta=[0, 0.24];
%     flap_R.xLE=[0.8,0.8];
%     flap_R.xTE=[1,1];
%     flap_R.Max_def=0.1;
%     flap_R.Max_rate=0.1;
%     flap_R.HingeLine='LE';
%     flap_R.Label='FlapR';
%     flap_R.FaceColor='r';
%     
% %     flap_R.NumAeroPanel=10;
%     flap_R.AeroPanelLength=0.4;
%     
%     build(flap_R)
%     Wingbox_right.add(flap_R);
%     
%     Wingbox_right.ModelControlSurf = 1;
    
    
    build(Wingbox_right);
    
%     draw(Wingbox_right)
%     FEM_test=convertToFE(Wingbox_right)
%     draw(FEM_test);
    


%% Wingbox 2 - left and control surf.

    Wingbox_left = awi.model.LiftingSurface;
    Wingbox_left.Name = 'A320Wing_left';
    Wingbox_left.Origin=[20,-2,0];
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

    %Make the material
%     E  = 70e9; %[N/m^2], typical YM of aluminium
%     nu = 0.333;
%     rho=2810*0.0001; %temp
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
    FWT_L.HingeStiffness = [1e14 1e14 1e14 1e14 1e-4 1e14];
    FWT_L.NumAeroPanel=15;
    
    FWT_L.NumBeamElem = 10;
   
    FWT_L.A   =  ones(1,sc_num).*FWT_box.Abb;
    % FWT.A_eta=FWT_eta;
    
    FWT_L.I11 = ones(1,sc_num).*FWT_box.Izz;
    % FWT.I11_eta=FWT_eta;
    
    FWT_L.I22 = ones(1,sc_num).*FWT_box.Ixx;
    % FWT.I22_eta = FWT_eta;
    
    FWT_L.J   = ones(1,sc_num).*FWT_box.Jbb;
    % FWT.J_eta= FWT_eta;
    
    
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

    % NSM and NSI: you will double the mass, nastran calculated mass based
    % on your input, no need to add extra mass components. 
    
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
    
    %Engine location 
    Y_eng=Wingbox_left.PanelCoords.LE.Y(2);
    X_eng=(Wingbox_left.PanelCoords.LE.X(2)+Wingbox_left.PanelCoords.TE.X(2))/2 - (Wingbox_left.PanelCoords.TE.X(2)-Wingbox_left.PanelCoords.LE.X(2))*0.1;
    Z_eng=(Wingbox_left.PanelCoords.LE.Z(2)+Wingbox_left.PanelCoords.TE.Z(2))/2;
    
    Engine2.Origin = [X_eng-3.5, Y_eng, Z_eng];
    
  
    Engine2.Length = 3.5;
%     Engine.XOffset=16.471008-3.5;
%     Engine.YOffset=5.8170588;
%     Engine.ZOffset=-2;
    
    
    
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
    
    %Control surfaces - flaps
%     flap_L=awi.model.ControlSurface;
%     flap_L.Eta=[0, 0.24];
%     flap_L.xLE=[0.8,0.8];
%     flap_L.xTE=[1,1];
%     flap_L.Max_def=0.1;
%     flap_L.Max_rate=0.1;
%     flap_L.HingeLine='LE';
%     flap_L.Label='FlapL';
%     flap_L.FaceColor='r';
%     
%     flap_L.AeroPanelLength=0.4;
%     
%     build(flap_R)
%     
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
    Body.Length=45;
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

    Body.Material=Body_Mat;
    Body.Material_eta = [0, 1];
    Body.Material     = [Body_Mat, Body_Mat];

    %define  panel size
    % Body.NumAeroPanel=5;
    Body.AeroPanelLength=0.5;

    % fuselage is assumed to be made by 4m diameter tube with wall
    % thickness of 4mm 
    
    Body_radius=2; 
    Body_thickness=0.004;
    CS_A=2*pi*Body_radius*Body_thickness;
    
    CS_I11=pi*Body_radius^3*Body_thickness;
    CS_I22=pi*Body_radius^3*Body_thickness;
    CS_J=2*pi*Body_radius^3*Body_thickness;
    CS_Inertia = Body_Mat.Rho*CS_A*Body_radius^2;
    CS_NSM = Body_Mat.Rho*CS_A;
    
    Body.A   = CS_A;
    Body.I11 = CS_I11;
    Body.I22 = CS_I22;
    Body.J   = CS_J;

    
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
    
%     Tailwing_right.NSM = 1000;
%     Tailwing_right.NSI = 100;

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
    myelevator_right.FaceColor='r';
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
    Tailwing_left.Dihedral=[5,5];
    Tailwing_left.Dihedral_eta=[0,1];

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
    Tailwing_left.A   = tailbox_left.Abb*0.0001; %temp
    Tailwing_left.I11 = tailbox_left.Ixx;
    Tailwing_left.I22 = tailbox_left.Izz;
    Tailwing_left.J   = tailbox_left.Jbb;

% 
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
    myelevator_left.FaceColor='r';
    
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
    Verticalwing.A   = Verticalbox.Abb; %temp
    Verticalwing.I11 = Verticalbox.Ixx;
    Verticalwing.I22 = Verticalbox.Izz;
    Verticalwing.J   = Verticalbox.Jbb;

% 
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



%     Body.add(Wingbox_right)
%     Body.add(Wingbox_left)
%     
%     Body.add(Tailwing_right)
%     Body.add(Tailwing_left)
%     Body.add(Verticalwing)
    


    %The analysis methods require an 'awi.model.Aircraft' object
    % This is because some information is only known at the aircraft level,
    % e.g. all-up mass, reference span, area, etc.
    % Aircraft = awi.model.Aircraft;
    % Aircraft.add(LS);

    Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea, Wingbox_left.SurfaceArea,...
        Connector_right.SurfaceArea,  Connector_left.SurfaceArea]);
    
    
    Aircraft.RefSpan  = Wingbox_right.Span*2+Connector_right.Span*2;
    Aircraft.RefChord = Wingbox_right.RootChord*Mean_cord_coefficient; %mean aerodynamic chord = 0.697 for A321 wing;

    
%     FEM_test2=convertToFE(Aircraft);
%     draw(Aircraft)
%     draw(FEM_test2)

    %% Trim loadcase 
    TrimLoadcase = awi.model.LoadCase;

    acMass = 500;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;

    TrimLoadcase.Name = 'A321_cruise_1g';
    TrimLoadcase.Altitude   = altitude;
    TrimLoadcase.Mach       = mach_number;
    TrimLoadcase.AcVelocity = aircraft_velocity;
    TrimLoadcase.AcMass = acMass;

    TrimLoadcase.PitchAngle=0;
    TrimLoadcase.RollAngle =0;
    TrimLoadcase.ID = 1030;
    TrimLoadcase.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase)

    %% Generate the FEM 

    run_folder = [
        'C:\Git\A321_sizing\hg_codes\results\test_temp']; %[-], folder for exporting the NASTRAN model

    % Convert to a finite element model and draw it
    FEM_full = convertToFE(Aircraft);

    % %Export it to a file
    export(FEM_full, run_folder);

  
     %% Run the analysis- SOL 144 static trim analysis 
    
%      NastranExe: 'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe'
     
     NastranMethods1 = awi.methods.Nastran;
     NastranMethods1.AnalysisModel = FEM_full;
     MassCases=awi.model.MassCases.empty;
     ID0=200;    
     
    trimFile = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase, MassCases,run_folder,'DatFilename','A321_range_calc');
   
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
      
    NastranMethods1.runNastran(trimFile);
    
    
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
    
    
   %% result extraction 
   
   % SOL 144
   
   Trim1_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_range_calc.f06'),'ReadF06',true,'ReadHDF5',false);
   
   WingNodes=NastranMethods1.WingNode;
   FWTNodes=NastranMethods1.FWTNode;
   
   X=[WingNodes.X];
   X_FWT=[FWTNodes.X];
   
   Y_Data=X(2,:);
   Y_Data_FWT=X_FWT(2,:);
   Y_all=[Y_Data(1:end-1),Y_Data_FWT];
   
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

   
   %SOL 146
   WingNodes_All=[WingNodes(1:end-1),FWTNodes];
   
   NumSteps=201;

   
   [Max_Moment_LC1,Min_Moment_LC1, Max_Torque_LC1, Min_Torque_LC1, Max_Shear_LC1, Min_Shear_LC1]=Gust_peaks(WingNodes_All,GustLoadcase1,run_folder,'\gust_analysis_g_36000ft_pos_1MC.h5',NumSteps);
   
   Max_Moment_LC1=Max_Moment_LC1';
   Min_Moment_LC1=Min_Moment_LC1';
   
   Max_Torque_LC1=Max_Torque_LC1';
   Min_Torque_LC1=Min_Torque_LC1';
   
   Max_Shear_LC1=Max_Shear_LC1';
   Min_Shear_LC1=Min_Shear_LC1';
   
   
 %% lift and drag calculateion   
 
%  % Induced drag and lift coefficient
% 
% Ajj=mni.result.op4(strcat(run_folder,'/ajj.op4')).read_matrix();
% FFaj=mni.result.op4(strcat(run_folder,'/ffaj.op4')).read_matrix();
% 
% res_aeroF = mni.result.f06(strcat(run_folder,'/A321_range_calc.f06')).read_aeroF;
% 
% % res_aeroF.aeroFz=res_aeroF.aeroFz;
% % 
% FlightPoint=awi.model.FlightPoint;
% FlightPoint.Mach=TrimLoadcase.Mach;
% FlightPoint.AcVelocity=FlightPoint.Mach*340;
% FlightPoint.Altitude = 36000;
% getFlightPointData(FlightPoint,'ISA');
% 
% DynPressure = FlightPoint.DynPressure;
% 
% WJ = Ajj*(FFaj./DynPressure);
% % surface_normal = model.fwt_normal_vector();
% % drag_mag = dot([1 0 0],surface_normal);
% 
% 
% % index for each lifting surfaces
% idx =407:958; % wing_right
% idx_c=1:48; %conn_right
% idx_t=97:184; %tail wing _right
% 
% idx1 =959:1510; % wing_left
% idx_c1=49:96; %conn_right
% idx_t1=188:278; %tail wing _left
% 
% % right side of AC
% lift_wing=res_aeroF.aeroFz(idx)';
% lift_conn=res_aeroF.aeroFz(idx_c)';
% lift_tail=res_aeroF.aeroFz(idx_t)';
% 
% % left side of AC
% lift_wing_L=res_aeroF.aeroFz(idx1)';
% lift_conn_L=res_aeroF.aeroFz(idx_c1)';
% lift_tail_L=res_aeroF.aeroFz(idx_t1)';
% 
% total_lift=sum(lift_wing)+sum(lift_conn)+sum(lift_tail)-sum(lift_wing_L)-sum(lift_conn_L)-sum(lift_tail_L);
% 
% total_area=(Wingbox_right.SurfaceArea+Connector_right.SurfaceArea+ Tailwing_right.SurfaceArea)*2;
% 
% % calculation of the lift coefficient: CL 
% Cl=total_lift/(DynPressure*total_area);
% 
% 
% % Induced drag - wing
% 
% % wing total on the wing
% lift_wing_left=[lift_conn_L;lift_wing_L];
% 
% lift_wing_left_matrix=reshape(lift_wing_left, 12, numel(lift_wing_left)/12);
% 
% % abs lift distribution on the span
% abs_lift_left=sum(lift_wing_left_matrix);
% 
% % normalise lift by panel width to obtain lift per unit span
% lift_wing_left_matrix(:,1:4)=lift_wing_left_matrix(:,1:4)/0.5;
% lift_wing_left_matrix(:,5:13)=lift_wing_left_matrix(:,5:13)/0.477;
% lift_wing_left_matrix(:,14:end)=lift_wing_left_matrix(:,14:end)/0.31371;
% 
% % lift per unit span
% lift_wing_left=sum(lift_wing_left_matrix);
% 
% %find Y-coords
% y0=0.25:0.5:1.75;
% y1=2.2385:0.477:6.0545;
% y2=6.449915:0.3137:17.74335;
% Y=[y0,y1,y2];
% 
% % Y VS lift
% Y_all=[-Y,Y];
% 
% %absolute lift
% abs_Lift_wing_all=[-abs_lift_left,-abs_lift_left];
% 
% %lift per unit span 
% Lift_wing_all=[-lift_wing_left,-lift_wing_left];
% 
% % circulation
% Gamma=Lift_wing_all/(FlightPoint.AirDensity*FlightPoint.AcVelocity);
% 
% % curve fitting
% Y_fit=linspace(-17.9,17.9,300);
% 
% P1=polyfit(Y_all,Gamma,10);
% % gamma_fit=polyval(P1,Y1_fit);
% P2 = polyder(P1);
% 
% wj=zeros(numel(Y_all),1);
% Fd=zeros(numel(Y_all),1);
% alphai=zeros(numel(Y_all),1);
% 
% for i=1:1:numel(Y_all)
%     
%     x0=Y_all(i);
%     fun_wj = @(x) (-1/(4*pi))*(P2(1)*x.^9+P2(2)*x.^8+P2(3)*x.^7+P2(4)*x.^6+ P2(5)*x.^5+ P2(6)*x.^4+ P2(7)*x.^3+ P2(8)*x.^2+ P2(9)*x +  P2(10))./(x0-x);
%     
%     wj_l = integral(fun_wj,-17.9, x0-0.00001);
%     wj_r = integral(fun_wj,x0+0.00001, 17.9);
%     
%     wj(i)= wj_l+wj_r;
%     alphai(i)=wj(i)/FlightPoint.AcVelocity;
%     Fd(i)=abs_Lift_wing_all(i)*sin(alphai(i));
%      
% end
% 
% % rm trivial data at the tip 
% %TODO: improve the way of integration/curve fiting
% 
% Fd(Fd>10)=0;
% Drag_force=sum(Fd);
% 
% Cdi=Drag_force/(DynPressure*126);
% 
% % figure 
% % plot(Y_all,alphai,'b.')
% % figure 
% % plot(Y_all,Lift_wing_all,'b.')
% % 
% % figure
% % plot(Y_all,abs_Lift_wing_all,'b.')
% 
% 
% 
% 
% % Parasite drag (zero-lift drag) Cd0
% 
% % wing 
% Re_c=1e6; % Ref. Reynold's number
% % Re=2.56e7; % Reynold's number
% mu=1.47e-5; % air viscosity 
% c=4.18; % mean aerodynamic chord
% rho=0.4; % air density 
% U=230; % free stream velocity
% Ma=0.78; % mach number
% x_tc=0.4; % relative chord position for the max. thickness 
% Q=1;    % Interface factor: 1 for wing 
% tc=Bheight./(2*Bwidth); % thickness to chord ratio along the span
% AreaR=2*(Bheight+2*Bwidth)/(2*Bwidth); % wetted area to reference area (planform area)
% Roughness=1.33e-5;
% 
% % Cut off Re number
% Xc=mu*Re_c/(rho*c*U); 
% Re_wingcut=38.21*(c/Roughness)^1.053;
% 
% % Reynold's number
% Re_wing_=rho*U*c/mu;
% Re_wing=min(Re_wing_,Re_wingcut);
% 
% % drag coefficient calculation
% C_lam=1.328/sqrt(Re_wing);
% C_turb=0.455/(log10(Re_wing)^2.58*(1+0.144*Ma^2)^0.65);
% Cf=Xc*C_lam+(1-Xc)*C_turb;
% 
% % Form factor to include the effect of pressure drag
% FF=(1+(0.6/x_tc)*tc+100*tc.^4).*(1.34*Ma^0.18*(cos(27*pi/180)^0.28));
% 
% Y_pt=etaS*Wingbox_right.Span;
% 
% Cd0=Cf.*FF.*Q.*AreaR;
% 
% % figure 
% % plot(Y_pt,Cd0,'b.')
% 
% drag_span=Cd0.*(2*Bwidth);
% 
% drag_force=area(drag_span,Y_pt);
% Wing_drag_force=trapz(Y_pt,drag_span);
% 
% wing_cd0=Wing_drag_force/Wingbox_right.SurfaceArea;
% 
% % fuselage 
% 
% Fuselage_length=45;
% Fuselage_diameter=4;
% FF_fuselage=1+60/(Fuselage_length/Fuselage_diameter)^3+Fuselage_length/(400*Fuselage_diameter);
% AreaR_fuselage=(2*pi*Fuselage_diameter*0.5*Fuselage_length)/(Fuselage_diameter*Fuselage_length);
% 
% % Cut off Re number
% Xc_fuse=mu*Re_c/(rho*Fuselage_length*U); 
% Re_fusecut=38.21*(Fuselage_length/Roughness)^1.053;
% 
% % Reynold's number
% Re_fuselage_=rho*U*Fuselage_length/mu;
% Re_fuselage=min(Re_fuselage_,Re_fusecut);
% 
% % drag coefficient calculation
% C_lam_fuselage=1.328/sqrt(Re_fuselage);
% C_turb_fuselage=0.455/(log10(Re_fuselage)^2.58*(1+0.144*Ma^2)^0.65);
% Cf_fuselage=Xc*C_lam_fuselage+(1-Xc_fuse)*C_turb_fuselage;
% 
% 
% Cd0f=Cf_fuselage.*FF_fuselage.*Q.*AreaR_fuselage;
% 
% % Engine 
% nacelles_length=3.5;
% nacelles_diameter=2;
% FF_nacelles=1+0.35/(nacelles_length/nacelles_diameter);
% Cf_nacelles=C_turb;
% AreaR_nacelles=(2*pi*nacelles_diameter*0.5*nacelles_length)/(nacelles_length*nacelles_diameter);
% 
% Q_en=1.3;
% Cd0n=Cf_nacelles.*FF_nacelles.*Q_en.*AreaR_nacelles;
% 
% % total zero-lift drag
% Cd0_total=(126*wing_cd0+180*Cd0f+14*Cd0n)/126;
% 
% % total_drag
% Cd_all=Cd0_total+abs(Cdi);




%% polar

% % A321
% CL=[0.6152
% 0.5987
% 0.5823
% 0.5658
% 0.5329
% 0.4835
% 0.4522
% 0.3864
% 0.3337
% 0.258
% 0.1622
% 0.08
% 0
% -0.08
% -0.1622
% -0.258
% -0.3337
% -0.3864
% -0.4522
% -0.4835
% -0.5329
% -0.5658
% -0.5823
% -0.5987
% -0.6152];
% 
% CD=[0.0314
% 0.0307
% 0.03
% 0.0293
% 0.0279
% 0.0261
% 0.025
% 0.023
% 0.0216
% 0.02
% 0.0184
% 0.0176
% 0.0174
% 0.0176
% 0.0184
% 0.02
% 0.0216
% 0.023
% 0.025
% 0.0261
% 0.0279
% 0.0293
% 0.03
% 0.0307
% 0.0314];
% 
% figure 
% 
% plot(CD,CL,'b-o','MarkerFaceColor','b','MarkerSize',4,'LineWidth',1)
% xlabel('Drag coefficient','Interpreter','latex','FontSize',12)
% ylabel('Lift coefficient','Interpreter','latex','FontSize',12)
% set(gcf,'Color','w')
% 
% figure 
% L_D=CL./CD;
% plot(CL,L_D,'b-o','MarkerFaceColor','b','MarkerSize',4,'LineWidth',1)
% xlabel('Lift coefficient','Interpreter','latex','FontSize',12)
% ylabel('Lift to drag ratio','Interpreter','latex','FontSize',12)
% set(gcf,'Color','w')

% Cdp


% load enso;
% f = fit(Y_all', Lift_wing_all', 'sin6');
% figure
% plot(f,Y_all', Lift_wing_all')
% figure
% plot(f)


%% elliptical test

% gamma0=100;
% b=30;
% 
% fun_e = @(x) gamma0*sqrt(1-(2.*x./b).^2);
% 
% x=-14.99:0.01:14.99;
% x_fit=-15:0.1:15;
% 
% gamma=fun_e(x);
% 
% P1=polyfit(x,gamma,12);
% gamma_fit=polyval(P1,x_fit);
% 
% % load enso;
% % f = fit(x', gamma', 'sin7');
% % 
% % a1 =       142.2;
% % b1 =       0.123;
% % c1 =       1.571;
% % a2 =       48.05;
% % b2 =       0.231;
% % c2 =      -1.571;
% % a3 =       13.89;
% % b3 =      0.4589;
% % c3 =       1.571;
% % a4 =       8.118 ;
% % b4 =      0.5147 ;
% % c4 =      -1.571 ;
% % a5 =      0.1311 ;
% % b5 =      0.9915 ;
% % c5 =       1.571 ;
% % a6 =     0.09931 ;
% % b6 =       1.227 ;
% % c6 =      -1.571  ;
% % a7 =      0.0904  ;
% % b7 =        1.51 ;
% % c7 =       1.571  ;
% % 
% % syms x
% % f = a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3) + a4*sin(b4*x+c4) + a5*sin(b5*x+c5) + a6*sin(b6*x+c6) + a7*sin(b7*x+c7);
% % 
% % ff=diff(f);
% % 
% figure 
% plot(x,gamma,'b.')
% hold on 
% plot(x_fit,gamma_fit,'b')
% % hold on 
% % plot(x_fit,fun_sin(x_fit))
% 
% 
% P2 = polyder(P1);
% 
% 
% wj=zeros(numel(x),1);
% 
% for i=1:1:numel(x)
%     
%     x0=x(i);
%     fun_wj = @(x) (-1/(4*pi))*(P2(1)*x.^11+P2(2)*x.^10+P2(3)*x.^9+P2(4)*x.^8+P2(5)*x.^7+P2(6)*x.^6+ P2(7)*x.^5+ P2(8)*x.^4+ P2(9)*x.^3+ P2(10)*x.^2+ P2(11)*x +  P2(12))./(x0-x);
%     
%     
%     wj_l = integral(fun_wj,-15, x0-0.0001);
%     wj_r = integral(fun_wj,x0+0.0001, 15);
%     
% 
%     wj(i)=wj_l+wj_r;
% %     alphai(i)=q(i)/V;
% %     Fd(i)=lift_var(i)*sin(alphai(i));
%      
% end
% 
% 
% figure 
% plot(x,wj,'b.')



 
     
     
     
     
    