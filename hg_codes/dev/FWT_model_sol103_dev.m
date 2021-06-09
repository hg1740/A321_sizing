  
%% sizing parameters thickness1 = spar, thickness2 = skin

    thickness1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
        x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
        x(21),x(22)];

    thickness2=[x(23),x(24),x(25),x(26),x(27),x(28),x(29),x(30),x(31),x(32)...
        x(33),x(34),x(35),x(36),x(37),x(38),x(39),x(40),x(41),x(42)...
        x(43),x(44)];

    Astrg=[x(45),x(46),x(47),x(48),x(49),x(50),x(51),x(52),x(53),x(54)...
        x(55),x(56),x(57),x(58),x(59),x(60),x(61),x(62),x(63),x(64)...
        x(65),x(66)];
    
   
  %% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[15,0,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = 15/0.75;   %34.1/2;
    Wingbox_right.LESweep     = [27, 27]*0.75;
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [0, 16.59, 16.59]*0.75;
    Wingbox_right.TESweep_eta = [0, 0.27, 1];
    Wingbox_right.RootChord   = 6;
    
    % testing non-sweep    
%     Wingbox_right.SpanVector  = 'Y';
%     Wingbox_right.Span        = 15;   %34.1/2;
%     Wingbox_right.LESweep     = [20, 20];
%     Wingbox_right.LESweep_eta = [0, 1];
%     Wingbox_right.TESweep     = [20, 20, 20];
%     Wingbox_right.TESweep_eta = [0, 0.27, 1];
%     Wingbox_right.RootChord   = 4;
% 
%     % testing single sweep    
%     Wingbox_right.SpanVector  = 'Y';
%     Wingbox_right.Span        = 15;   %34.1/2;
%     Wingbox_right.LESweep     = [0, 0];
%     Wingbox_right.LESweep_eta = [0, 1];
%     Wingbox_right.TESweep     = [0, 0, 0];
%     Wingbox_right.TESweep_eta = [0, 0.27, 1];
%     Wingbox_right.RootChord   = 5;
    
    %Dihedral 
    Wingbox_right.Dihedral=[0,0];
    Wingbox_right.Dihedral_eta=[0,1];
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_right.Eta_;
    Wingbox_right.BeamLoc     = repmat(0.5, size(all_eta));
%     Wingbox_right.BeamLoc     = [0.34,0.4,0.4];
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
    rho=2810;
    Mat = awi.model.Material;
    Mat.E  = E;
    Mat.Nu = nu;
    Mat.G  = E / (2 * (1 + nu));
    Mat.Rho=rho;
    Wingbox_right.Material_eta = [0, 1];
    Wingbox_right.Material     = [Mat, Mat];
    
    build(Wingbox_right)
    
    
%% FWT definition 
% fold_angle  = -10;   %[deg],
% flare_angle = -10;   %[deg],
% FWT = insertWingFold(Wingbox_right, 'FlareAngle', flare_angle, 'FoldAngle', fold_angle,'EtaFold',0.75);
% FWT.HingeStiffness = [0 0 0 1e-4 0 0];
% FWT.NumAeroPanel=20;
% % 
% % 
% % 
% % 
% % FWT_box=awi.model.BoxBeam;
% % FWT_box.BoxType='SymmetricBox';
% % FWT_box.Height=0.5;
% % FWT_box.Width=0.5;
% % FWT_box.CoverThickness=2;
% % FWT_box.SparThickness=2;
% % 
% % getGeometricProps(boxname)
% % 
% % FWT_eta=[0,1];
% % 
% % FWT.A   =  FWT_box.Abb;
% % % FWT.A_eta=FWT_eta;
% % 
% % FWT.I11 = Izz_val;
% % % FWT.I11_eta=FWT_eta;
% % 
% % FWT.I22 = Ixx_val;
% % % FWT.I22_eta = FWT_eta;
% % 
% % FWT.J   = J_val;
% % % FWT.J_eta= FWT_eta;
% % 
% % % NSM and NSI
% % FWT.NSM=NSM_val;
% % % FWT.NSM_eta= FWT_eta;
% % 
% % FWT.NSI=NSI_val;
% % % FWT.NSI_eta= FWT_eta;
% %         
% %         
% % % add point masses
% % for i=1:1:10
% %     handle=strcat('PM_right','i');
% %     handle=awi.model.PointMass;
% %      handle.SOffset=0+i*0.1;
% %      handle.Mass=0.5;
% %      handle.Inertia11 = 0.1;
% %      handle.Inertia22 =  0.1;
% %      handle.Inertia33 =  0.1;
% %      handle.Inertia12 =  0.1;
% %      handle.Inertia23 =  0.1;
% %      handle.Inertia13 =  0.1;
% %      handle.MassGroup='Group1'; 
% %      FWT.add(handle);
% %     
% % end

    
    
    %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=21;
    d_strg=sqrt(Astrg/0.36);
    t_strg=0.12*d_strg;

    % etaS=linspace(0,Wingbox_right.Span,NumSec);

    % set width and height array 
    YData=Wingbox_right.YData;
    SparWidth=Wingbox_right.Chord*0.5;

    RootH=Wingbox_right.Chord(1)*0.15;
    MidH=Wingbox_right.Chord(2)*0.12;
    TipH=Wingbox_right.Chord(end)*0.11;


    % set up eta values
    eta1_=linspace(0,0.27,8);
    eta2_=linspace(0.27,1,14);
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
    
    % offset
    
%     Wingbox_right.SCy=
    
    
    
    % Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
    Wingbox_right.NumAeroPanel=20;
    
    build(Wingbox_right)
    
    
    
 %% FWT definition 
fold_angle  = -10;   %[deg],
flare_angle = -15;   %[deg],
FWT = insertWingFold(Wingbox_right, 'FlareAngle', flare_angle, 'FoldAngle', fold_angle,'EtaFold',0.75);
FWT.HingeStiffness = [1e14 1e14 1e14 1e14 1e14 1e14];
FWT.NumAeroPanel=20;

FWT_box=awi.model.BoxBeam;
FWT_box.BoxType='SymmetricBox';
FWT_box.Height=0.5;
FWT_box.Width=0.5;
FWT_box.CoverThickness=0.002;
FWT_box.SparThickness=0.002;

getGeometricProps(FWT_box)

sc_num=numel(FWT.A);
% FWT_eta=[0,1];

FWT.A   =  ones(1,sc_num).*FWT_box.Abb;
% FWT.A_eta=FWT_eta;

FWT.I11 = ones(1,sc_num).*FWT_box.Izz;
% FWT.I11_eta=FWT_eta;

FWT.I22 = ones(1,sc_num).*FWT_box.Ixx;
% FWT.I22_eta = FWT_eta;

FWT.J   = ones(1,sc_num).*FWT_box.Jbb;
% FWT.J_eta= FWT_eta;

% NSM and NSI
FWT.NSM = 1.*ones(1,sc_num).*FWT_box.NSI;
% FWT.NSM_eta= FWT_eta;

FWT.NSI = 1*ones(1,sc_num).*FWT_box.NSM;

%Make the material for FWT

E_fwt  = 70e9; %[N/m^2], typical YM of aluminium
nu_fwt = 0.333;
rho_fwt=2810;
Mat_fwt = awi.model.Material;
Mat_fwt.E  = E_fwt;
Mat_fwt.Nu = nu_fwt;
Mat_fwt.G  = E_fwt / (2 * (1 + nu_fwt));
Mat_fwt.Rho=rho_fwt;
FWT.Material_eta = [0, 1];
FWT.Material     = [Mat_fwt, Mat_fwt];



% FWT.NSI_eta= FWT_eta;
%         
% %         
% % % add point masses
% % for i=1:1:10
% %     handle=strcat('PM_right','i');
% %     handle=awi.model.PointMass;
% %      handle.SOffset=0+i*0.1;
% %      handle.Mass=0.5;
% %      handle.Inertia11 = 0.1;
% %      handle.Inertia22 =  0.1;
% %      handle.Inertia33 =  0.1;
% %      handle.Inertia12 =  0.1;
% %      handle.Inertia23 =  0.1;
% %      handle.Inertia13 =  0.1;
% %      handle.MassGroup='Group1'; 
% %      FWT.add(handle);
% %     
% % end   
    

    %% Mass definition
    
    % total wing mass
    
    [wing_mass,total_mass]=Mass_calc_v1(x);
    
    wingmass_eta=0.04:0.04:0.96;
%     wingmass_eta=0.2:0.2:1;
    Mwidth=interp1(RData/RData(end),SparWidth,wingmass_eta);
    
    mass_set=(total_mass-wing_mass)*(Mwidth)/sum(Mwidth);
      
%     m=total_mass/19;
    
    for i=1:1:24
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
    
%     % Landing Gear
%     Landing_Gear=awi.model.PointMass;   
%     Landing_Gear.SOffset=0.23;
%     Landing_Gear.Mass=2491/2;
%     Wingbox_right.add(Landing_Gear);
%     
%     % engine
%     engine=awi.model.PointMass;   
%     engine.SOffset=0.25;
% %     engine.XOffset=-1.5;
% %     engine.YOffset=4;
% %     engine.ZOffset=-1.8;
%     engine.Mass=3681+1239/2;
%     Wingbox_right.add(engine);
    
 
    %% attachments - engine
    
    Engine=awi.model.BluffBody;
    Engine.Name='Engine';
    
    % cylinder body
    Engine.Radius=[1.4, 1.4, 1];
    Engine.Eta =  [0, 0.6, 1];
    
%     Engine.Origin = [16.238147-3.5, 4.05 , -1.8];
    Engine.Origin = [16.238147-3.5, 5.4 , 0.35432909];
  
    Engine.Length = 3.5;
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
    engine_mass.SOffset=0.4;
    engine_mass.Mass=3681;
    Engine.add(engine_mass);
    
    % add pylon
    pylon_mass=awi.model.PointMass;   
    pylon_mass.SOffset=0.9;
    pylon_mass.Mass=1239/2;
    Engine.add(pylon_mass);

    build(Engine)
      
    Wingbox_right.add(Engine)
    build(Wingbox_right);
    
    draw(Wingbox_right)
    

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
    
%     FEM_body=convertToFE(Body)


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

    Aircraft.add(Wingbox_right);
    
%     Body.add(Engine)
    
%     Body.add(Wingbox_right)
%     Body.add(Tailwing_right)
%     Body.add(Verticalwing)
    


    %The analysis methods require an 'awi.model.Aircraft' object
    % This is because some information is only known at the aircraft level,
    % e.g. all-up mass, reference span, area, etc.
    % Aircraft = awi.model.Aircraft;
    % Aircraft.add(LS);

    Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea]);
    Aircraft.RefSpan  = Wingbox_right.Span;
%     Aircraft.RefChord = Wingbox_right.RootChord*0.697;
    Aircraft.RefChord =  Aircraft.RefArea/Aircraft.RefSpan;


    %% Generate the loadcase object
    TrimLoadcase = awi.model.LoadCase;

    acMass = 500;
    altitude          = 10000;
    mach_number       = 0.78;
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
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\FWT_SOL103_test1']; %[-], folder for exporting the NASTRAN model


    % Convert to a finite element model and draw it
    FEM_half = convertToFE(Aircraft);

    % %Export it to a file
    export(FEM_half, run_folder);


    %% Run the analysis- SOL 103 static trim analysis

      % insert SPC
    fid = fopen('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\FWT_SOL103_test1\NastranHeaderFile.dat');
    
    cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true );
    cac=cac{1};
    fclose( fid )
    
    mid_line0=22;
    mid_line=length(cac)-1; 
    fid = fopen( 'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\FWT_SOL103_test1\FWT_103_T1.dat', 'w' );
    
    % write existing head file
    for jj = 1 : mid_line0
        fprintf(fid,'%s\n',cac{jj})       
    end
    
    fprintf(fid,'SPC = 1 \r\n')
    
    for ii =  mid_line0 : mid_line 
        fprintf(fid,'%s\n',cac{ii})       
    end
    
    
    %write spc part
    line='$.1.....2.......3.......4.......5.......6.......7.......8.......9.......10......\r\n';
    fprintf(fid,line);
    
    spc_format='%-8s%-8i%-8i%-8i\r\n';
    fprintf(fid,spc_format,'SPC1',1,123456,1008);
    
    % write the end 
    fprintf(fid,'%s\n',cac{end})
    
    