%% SOL 144
function [Von_skn, Von_spr, sigmab_skn, taub_skn, sigma_pp, sigma_strg, sigma_crip, sigma_col]=BeamStress_calc_v0(x)

    
    thickness1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
        x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
        x(21),x(22),x(23),x(24),x(25)];

    thickness2=[x(26),x(27),x(28),x(29),x(30),x(31),x(32),x(33),x(34),x(35)...
        x(36),x(37),x(38),x(39),x(40),x(41),x(42),x(43),x(44),x(45)...
        x(46),x(47),x(48),x(49),x(50)];

    Astrg=[x(51),x(52),x(53),x(54),x(55),x(56),x(57),x(58),x(59),x(60)...
        x(61),x(62),x(63),x(64),x(65),x(66),x(67),x(68),x(69),x(70)...
        x(71),x(72),x(73),x(74),x(75)];
    
   
    %% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[15,0,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = 34.1/2;
    Wingbox_right.LESweep     = [25, 25];
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [0, 0, 18];
    Wingbox_right.TESweep_eta = [0, 0.3, 1];
    Wingbox_right.RootChord   = 6;

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
    E  = 70e9; %[N/m^2], typical YM of aluminium
    nu = 0.333;
    Mat = awi.model.Material;
    Mat.E  = E;
    Mat.Nu = nu;
    Mat.G  = E / (2 * (1 + nu));
    Wingbox_right.Material_eta = [0, 1];
    Wingbox_right.Material     = [Mat, Mat];
    
    build(Wingbox_right)
    
    
    %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=25;
    d_strg=sqrt(Astrg/0.36);
    t_strg=0.12*d_strg;

    % etaS=linspace(0,Wingbox_right.Span,NumSec);

    % set width and height array 
    YData=Wingbox_right.YData;
    SparWidth=Wingbox_right.Chord*0.5;

    RootH=Wingbox_right.Chord(1)*0.13;
    TipH=Wingbox_right.Chord(end)*0.09;


    % set up eta values
    eta1_=linspace(0,0.3,8);
    eta2_=linspace(0.3,1,18);
    etaS=[eta1_(1:end-1),eta2_(1:end)];

    RData=Wingbox_right.RData;
    eta_R=RData/RData(end);
    eta_Y=YData/YData(end);
    etaRS=interp1(eta_Y,eta_R,etaS);

    Bwidth=interp1(RData/RData(end),SparWidth,etaRS);
    Bheight=interp1([0,1],0.79*[RootH,TipH],etaRS);

    % stringer pitch 
    strg_n=0.24;

    %intialise data array
    A_val=zeros(1,numel(thickness1));
    Ixx_val=zeros(1,numel(thickness1));
    Izz_val=zeros(1,numel(thickness1));
    J_val=zeros(1,numel(thickness1));
    

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

    end


    eta_=etaRS;
    Wingbox_right.A   =  A_val;
    Wingbox_right.A_eta=eta_;

    Wingbox_right.I11 = Ixx_val;
    Wingbox_right.I11_eta=eta_;

    Wingbox_right.I22 = Izz_val;
    Wingbox_right.I22_eta = eta_;

    Wingbox_right.J   = J_val;
    Wingbox_right.J_eta= eta_;

  
    % Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
    Wingbox_right.NumAeroPanel=20;
   

    %% Mass definition
    
    % total wing mass
    
    [wing_mass,total_mass]=Mass_calc(x);
    
    m=total_mass/12;
    
    for i=1:1:12
        handle=strcat('PM_right','i');
        handle=awi.model.PointMass;
        handle.SOffset=0+i*0.08;
        handle.Mass=m;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle.MassGroup='Group1';
        Wingbox_right.add(handle);

    end
    
    % attachments
    engine_mass=awi.model.PointMass;   
    engine_mass.SOffset=0.32;
    engine_mass.Mass=3681+1239/2;
    Wingbox_right.add(engine_mass);
    
    Landing_Gear=awi.model.PointMass;   
    Landing_Gear.SOffset=0.23;
    Landing_Gear.Mass=2491/2;
    Wingbox_right.add(Landing_Gear);


    build(Wingbox_right);


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
    Body.Length=37.57;
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
%     Body.NSI = Bodybox.NSI;

    % add point masses- fuselage 
    M_fs=30000;
    m_fs=M_fs/11;
    
    %payload
    M_p=13608;
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

    %Make sure the beam is at the midchord
    all_eta           = Tailwing_right.Eta_;
    Tailwing_right.BeamLoc     = repmat(0.5, size(all_eta));
    Tailwing_right.BeamLoc_eta = all_eta;
    Tailwing_right.XOffset=35;

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

    % Aeropanel definition
    Tailwing_right.AeroPanelLength=0.3;

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
    build(myelevator_right)
    Tailwing_right.add(myelevator_right);

    Tailwing_right.ModelControlSurf = 1;


    build(Tailwing_right);


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
    Verticalwing.XOffset=35;


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

    build(Verticalwing);


    %% Build aircraft model
    Aircraft = awi.model.Aircraft;

    Aircraft.add(Body);
    Body.add(Wingbox_right)
    Body.add(Tailwing_right)
    Body.add(Verticalwing)


    %The analysis methods require an 'awi.model.Aircraft' object
    % This is because some information is only known at the aircraft level,
    % e.g. all-up mass, reference span, area, etc.
    % Aircraft = awi.model.Aircraft;
    % Aircraft.add(LS);

    Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea]);
    Aircraft.RefSpan  = Wingbox_right.Span;
    Aircraft.RefChord = Wingbox_right.RootChord;


    %% Generate the loadcase object
    TrimLoadcase = awi.model.LoadCase;

    acMass = 500;
    altitude          = 36000;
    mach_number       = 0.7;
    aircraft_velocity = 250;

    TrimLoadcase.Name = 'A320_half_model_SOL144';
    % TrimLoadcase.LoadCaseTypes = 'Static';(read only)
    % TrimLoadcase.CsDeflecTypes='fixed';(read only)
    TrimLoadcase.Altitude   = altitude;
    TrimLoadcase.Mach       = mach_number;
    TrimLoadcase.AcVelocity = aircraft_velocity;
    TrimLoadcase.AcMass = acMass;

    TrimLoadcase.PitchAngle=0;
    TrimLoadcase.RollAngle =0;
    TrimLoadcase.ID = 1020;
    TrimLoadcase.LoadFactor = 2.5;
    build(TrimLoadcase)

    %% Generate the FEM - SOL 103 Modal analysis

    run_folder = [
        'D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12']; %[-], folder for exporting the NASTRAN model


    % Convert to a finite element model and draw it
    FEM_half = convertToFE(Aircraft);
%     draw(FEM_half);

    % %Export it to a file
    export(FEM_half, run_folder);


    %% Run the analysis- SOL 144 static trim analysis

    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_half;
    MassCases=awi.model.MassCases.empty;
    ID0=200;
    % trimdata=NastranMethods1.getTrimData(FEM_half, Aircraft, TrimLoadcase, MassCases, ID0, NastranMethods1.RefNode);
    trimFile = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase, MassCases,run_folder);
    
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\A320_half_model_SOL144*.*')
   
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.xdb');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.h5');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.log');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.f06');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.f04');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.h5.*');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.log.*');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.f06.*');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144*.f04.*');
      
    NastranMethods1.runNastran(trimFile);
    
    
    %% extract results

    F061=NastranMethods1.extractNastranResults('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144.f06','ReadF06',true,'ReadHDF5',false);

    % extract forces
    M_P1=[F061.f06data.Bendingmoment.UMPLN1(69:92),F061.f06data.Bendingmoment.LMPLN1(92)]; % in-plane moment
    M_P2=[F061.f06data.Bendingmoment.UMPLN2(69:92),F061.f06data.Bendingmoment.LMPLN2(92)]; % out of plane moment
    T=[F061.f06data.Bendingmoment.UTORQUE1(69:92),F061.f06data.Bendingmoment.LTORQUE1(92)];% torque

    S_P1=[F061.f06data.Bendingmoment.USPLN1(69:92),F061.f06data.Bendingmoment.LSPLN1(92)]; % in plane shear 
    S_P2=[F061.f06data.Bendingmoment.USPLN2(69:92),F061.f06data.Bendingmoment.LSPLN2(92)]; % out of plane shear

    % calculate stresses


    sigma_skn=ones(1,25);
    tau_skin=ones(1,25);

    sigma_spr=ones(1,25); 
    tau_spr=ones(1,25);

    Von_skn=ones(1,25);
    Von_spr=ones(1,25);

    sigmab_skn=ones(1,25);
    taub_skn=ones(1,25);

    sigma_strg=ones(1,25);
    sigma_crip=ones(1,25);
    sigma_col=ones(1,25);
    sigma_pp=ones(1,25);

    for jj=1:25

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
    
%     k=4;
% 
%     sigma_b=k*(pi)^2*E*(mybox_right.CoverThickness/mybox_right.Width)^2/(12*(1-0.33^2));
    
    
    
    
% end

% Trim_res=NastranMethods1.trim(Aircraft, TrimLoadcase);


%% Run flutter - SOL 145 flutter analysis
% 
% FlightPoint=awi.model.FlightPoint;
% 
% FlightPoint.Mach=0.5;
% FlightPoint.Altitude = 36000;
% 
% getFlightPointData(FlightPoint)
% 
% flutterFile = NastranMethods1.writeFlutterFile(Aircraft, FlightPoint, MassCases, '123456', run_folder)
% 
% NastranMethods1.runNastran(flutterFile);


%% SOL 146 - Guest

% t11=x(1);t12=x(2);t13=x(3);t14=x(4);t15=x(5);
%     t16=x(6);t17=x(7);t18=x(8);t19=x(9);t110=x(10);
%     t111=x(11);t112=x(12);t113=x(13);t114=x(14);t115=x(15);
%     t116=x(16);t117=x(17);t118=x(18);t119=x(19);t120=x(20);
%     t121=x(21);t122=x(22);t123=x(23);t124=x(24);t125=x(25);
% 
%     
%     t21=x(26);t22=x(27);t23=x(28);t24=x(29);t25=x(30);
%     t26=x(31);t27=x(32);t28=x(33);t29=x(34);t210=x(35);
%     t211=x(36);t212=x(37);t213=x(38);t214=x(39);t215=x(40);
%     t216=x(41);t217=x(42);t218=x(43);t219=x(44);t220=x(45);
%     t221=x(46);t222=x(47);t223=x(48);t224=x(49);t225=x(50);
%     
%     
%     Ag1=x(51);Ag2=x(52);Ag3=x(53);Ag4=x(54);Ag5=x(55);
%     Ag6=x(56);Ag7=x(57);Ag8=x(58);Ag9=x(59);Ag10=x(60);
%     Ag11=x(61);Ag12=x(62);Ag13=x(63);Ag14=x(64);Ag15=x(65);
%     Ag16=x(66);Ag17=x(67);Ag18=x(68);Ag19=x(69);Ag20=x(70);
%     Ag21=x(71);Ag22=x(72);Ag23=x(73);Ag24=x(74);Ag25=x(75);
%     
%         
%     thickness1=[t11,t12,t13,t14,t15...
%         t16,t17,t18,t19,t110...
%         t111,t112,t113,t114,t115...
%         t116,t117,t118,t119,t120...
%         t121,t122,t123,t124,t125]; 
%     
%     thickness2=[t21,t22,t23,t24,t25...
%         t26,t27,t28,t29,t210...
%         t211,t212,t213,t214,t215...
%         t216,t217,t218,t219,t220...
%         t221,t222,t223,t224,t225];
%  
%    
%     Astrg=[Ag1,Ag2,Ag3,Ag4,Ag5,Ag6,Ag7,Ag8,Ag9,Ag10...
%         Ag11,Ag12,Ag13,Ag14,Ag15,Ag16,Ag17,Ag18,Ag19,Ag20...
%         Ag21,Ag22,Ag23,Ag24,Ag25];



























