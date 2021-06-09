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
    Wingbox_right.Origin=[15,0,0]; %15
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';
%      Wingbox_right.BeamLength = 0.5;
     Wingbox_right.NumBeamElem = 25;

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = 15;   %34.1/2;
    Wingbox_right.LESweep     = [27, 27];
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [0, 16.59, 16.59];
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
    Wingbox_right.Dihedral=[5,5];
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
    
    
    %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=22;
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
    eta1_=linspace(0,0.27,7);
    eta2_=linspace(0.27,1,16);
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
   

    %% Mass definition
    
    % total wing mass
    
    [wing_mass,total_mass]=Mass_calc_v1(x);
    
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
    
%     Engine.Origin = [16.238147-3.5, 4.05 , -1.8];
    Engine.Origin = [16.238147-3.5, 4.05 , 0.35432909];
  
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
    
    FEM_test=convertToFE(Wingbox_right);
    
    draw(FEM_test)
    
    run_folder = [
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\wing_test']; %[-], folder for exporting the NASTRAN model


    % %Export it to a file
    export(FEM_test, run_folder);