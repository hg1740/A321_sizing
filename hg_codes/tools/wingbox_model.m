    %% sizing parameters thickness1 = spar, thickness2 = skin   

%     x=[0.0191135910000000,0.0193953730000000,0.0201137600000000,0.0213823650000000,0.0229999860000000,0.0255251240000000,0.0285624380000000,0.0372113470000000,0.0373207110000000,0.0372799900000000,0.0368848460000000,0.0378039760000000,0.0379687300000000,0.0379303110000000,0.0378464170000000,0.0368959070000000,0.0365928770000000,0.0358741000000000,0.0351157990000000,0.0325204110000000,0.0321366360000000,0.0317223020000000,0.0280889430000000,0.0226448800000000,0.0118821390000000,0.00707914700000000,0.00718361400000000,0.00730636500000000,0.00745189900000000,0.00763438400000000,0.00783787200000000,0.00810313600000000,0.00823100800000000,0.00803824900000000,0.00781619300000000,0.00760075000000000,0.00736034800000000,0.00709415700000000,0.00680991600000000,0.00650766900000000,0.00618441800000000,0.00597655500000000,0.00574729800000000,0.00545697300000000,0.00516124700000000,0.00473319700000000,0.00418339100000000,0.00343190400000000,0.00256919700000000,0.000896999000000000,6.20000000000000e-05,6.38000000000000e-05,6.58000000000000e-05,6.82000000000000e-05,7.14000000000000e-05,7.49000000000000e-05,7.99000000000000e-05,8.21000000000000e-05,7.78000000000000e-05,7.33000000000000e-05,6.92000000000000e-05,6.48000000000000e-05,6.03000000000000e-05,5.55000000000000e-05,5.06000000000000e-05,4.55000000000000e-05,4.24000000000000e-05,3.92000000000000e-05,3.52000000000000e-05,3.13000000000000e-05,2.63000000000000e-05,2.05000000000000e-05,1.40000000000000e-05,8.02000000000000e-06,9.15000000000000e-07];
%     
    x=0.001*ones(1,75);
    %% Wing configurations for starboard wing
  
     Aspect_ratio=18.172; % Aspect ratio = 10.172 for A321 model
    
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

%     build(Wingbox_right);
    
%     draw(Wingbox_right)
%     FEM_test=convertToFE(Wingbox_right)
%     draw(FEM_test);
%         
    
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
    
%     Engine=awi.model.BluffBody;
%     Engine.Name='Engine';
%     
%     % cylinder body
%     Engine.Radius=[1.4, 1.4, 1];
%     Engine.Eta =  [0, 0.6, 1];
%     
%     %Engine location 
%     Y_eng=Wingbox_right.PanelCoords.LE.Y(2);
%     X_eng=(Wingbox_right.PanelCoords.LE.X(2)+Wingbox_right.PanelCoords.TE.X(2))/2 - (Wingbox_right.PanelCoords.TE.X(2)-Wingbox_right.PanelCoords.LE.X(2))*0.1;
%     Z_eng=(Wingbox_right.PanelCoords.LE.Z(2)+Wingbox_right.PanelCoords.TE.Z(2))/2;
%     
%     Engine.Origin = [X_eng-3.5, Y_eng, Z_eng];
%   
%     Engine.Length = 3.5;
%     
%     %Make engine material
%     E1  = 76e9; %[N/m^2],set as a rigid body
%     nu = 0.333;
%     Engine_Mat = awi.model.Material;
%     Engine_Mat.E  = E1;
%     Engine_Mat.Nu = nu;
%     Engine_Mat.G  = E1 / (2 * (1 + nu));
%     Engine_Mat.Rho=1; % using lumped mass instead
%     
%     
%     % use the strong material
%     Engine.Material_eta = [0, 1];
%     Engine.Material     = [Engine_Mat, Engine_Mat];
%     
%     % Engine stiffness
%     Engine_radius=1;
%     Engine_thickness=0.015;
%     Engine_A=2*pi*Engine_radius*Engine_thickness; 
%     
%     Engine_I11=pi*Engine_radius^3*Engine_thickness;
%     Engine_I22=pi*Engine_radius^3*Engine_thickness;
%     Engine_J=2*pi*Engine_radius^3*Engine_thickness;
%     Engine_Inertia = Engine_Mat.Rho*Engine_A*Engine_radius^2;
%     
%     Engine.A   = Engine_A*0.00001; % temp
%     Engine.I11 = Engine_I11;
%     Engine.I22 = Engine_I22;
%     Engine.J   = Engine_J;
%     %     Body.NSM = Bodybox.NSM;
% %     Engine.NSI = Engine_Inertia;
% 
% %     Engine.A   = 0.04432;
% %     Engine.I11 = 0.002;
% %     Engine.I22 = 0.002;
% %     Engine.J   = 0.001636;
%     
%     %Aeropanel althoufh it is useless now
%     Engine.AeroPanelLength=0.5;
%     
%     % add engine mass
%     engine_mass=awi.model.PointMass;   
%     engine_mass.SOffset=0.1;
%     engine_mass.Mass=Engine_mass;
%     Engine.add(engine_mass);
%     
%     % add pylon
%     pylon_mass=awi.model.PointMass;   
%     pylon_mass.SOffset=0.9;
%     pylon_mass.Mass=Pylon;
%     Engine.add(pylon_mass);
% 
%     build(Engine)
%       
%     Wingbox_right.add(Engine)
    
    
    
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
    
    draw(Wingbox_right)
    FEM_test=convertToFE(Wingbox_right)
    draw(FEM_test);