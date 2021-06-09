  
%% sizing parameters thickness1 = spar, thickness2 = skin



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

    thickness1=x0(1:22);
    thickness2=x0(23:44); 
    Astrg=x0(45:66);
   
  %% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[15,0,0]; %15
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';

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
    
    [wing_mass,total_mass]=Mass_calc_v1(x0);
    
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
    flap_R.Eta=[0, 0.24];
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
    Wingbox_left.Origin=[15,0,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_left.ActiveSet = 'sSet';

    %Wing dimensions
    Wingbox_left.SpanVector  = 'Y';
    Wingbox_left.Span        = -15;   %34.1/2;
    Wingbox_left.LESweep     = [-27, -27];
    Wingbox_left.LESweep_eta = [0, 1];
    Wingbox_left.TESweep     = [0, -16.59, -16.59];
    Wingbox_left.TESweep_eta = [0, 0.27, 1];
    Wingbox_left.RootChord   = 6;
    
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

%     NumSec=22;
%     d_strg=sqrt(Astrg/0.36);
%     t_strg=0.12*d_strg;
% 
%     % etaS=linspace(0,Wingbox_right.Span,NumSec);
% 
%     % set width and height array 
%     YData=Wingbox_right.YData;
%     SparWidth=Wingbox_right.Chord*0.5;
% 
%     RootH=Wingbox_right.Chord(1)*0.15;
%     MidH=Wingbox_right.Chord(2)*0.12;
%     TipH=Wingbox_right.Chord(end)*0.11;
% 
% 
%     % set up eta values
%     eta1_=linspace(0,0.27,7);
%     eta2_=linspace(0.27,1,16);
%     etaS=[eta1_(1:end-1),eta2_(1:end)];
% 
%     RData=Wingbox_right.RData;
%     eta_R=RData/RData(end);
%     eta_Y=YData/YData(end);
%     etaRS=interp1(eta_Y,eta_R,etaS);
% 
%     Bwidth=interp1(RData/RData(end),SparWidth,etaRS);
%     Bheight=interp1(RData/RData(end),0.79*[RootH,MidH,TipH],etaRS);
% %      Bheight=interp1(RData/RData(end),0.79*[RootH,RootH,RootH],etaRS);
% 
%     % stringer pitch 
%     strg_n=0.24;
% 
%     %intialise data array
%     A_val=zeros(1,NumSec);
%     Ixx_val=zeros(1,NumSec);
%     Izz_val=zeros(1,NumSec);
%     J_val=zeros(1,NumSec);
%     
%     %NSM
%     NSM_val=zeros(1,NumSec);
%     NSI_val=zeros(1,NumSec);
%     
%     %offset from shear center
%     SCy_val=zeros(1,NumSec);
%     SCz_val=zeros(1,NumSec);
%     NAy_val=zeros(1,NumSec);
%     NAz_val=zeros(1,NumSec);
%     CMy_val=zeros(1,NumSec);
%     CMz_val=zeros(1,NumSec);
%     
%     %offset from CoG
%     xOff=interp1(YData/YData(end),Wingbox_right.XData,etaS);
%     xOff_1=xOff(1:end-1);
%     xOff_2=xOff(2:end);
%     xOff_val=[0,xOff_2-xOff_1];
%     
% 
%     for ii=1:NumSec
% 
%         boxname=strcat('Box',string(ii));
%         boxname=awi.model.BoxBeam;
%         boxname.BoxType='SymmetricBox';
%         boxname.Height=Bheight(ii);
%         boxname.Width=Bwidth(ii);
%         boxname.CoverThickness=thickness2(ii);
%         boxname.SparThickness=thickness1(ii);
% 
%         NumStrg=floor(Bwidth(ii)/strg_n);
% 
%         ts=t_strg(ii);
%         ds=d_strg(ii);
%         hs=Bheight(ii);
%         ws=Bwidth(ii);
% 
%         Istrg_xx=((ts*ds^3/12)+(ds*ts^3/12 + ds*ts*(ds/2)^2)*2 + 3*ds*ts*(hs/2-ds/2)^2)*NumStrg*2;
%         Istrg_zz_=(ds*ts^3/12 + (ts*ds^3/12 + ts*ds*(ds/2)^2)*2);
% 
%         if mod(NumStrg,2)==0
%             offset=0.12:strg_n:ws/2;
%             Istrg_zz=(NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;
% 
%         elseif mod(NumStrg,2)==1
%             offset=0:strg_n:ws/2;
%             Istrg_zz=(NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;
% 
%         end
% 
%         getGeometricProps(boxname)
%         A_val(ii)=boxname.Abb+0;
%         Ixx_val(ii)=boxname.Ixx+Istrg_xx;
%         Izz_val(ii)=boxname.Izz+Istrg_zz;
%         J_val(ii)=boxname.Jbb;
%         
%         % NSM
%        
%         NSM_val(ii)=boxname.NSM;
%         NSI_val(ii)=boxname.NSI;
%         
%         % offset
%         SCy_val(ii)=boxname.xSC;
%         SCz_val(ii)=boxname.zSC;
%         NAy_val(ii)=boxname.xNA;
%         NAz_val(ii)=boxname.zNA;
%         CMy_val(ii)=boxname.xCM;
%         CMz_val(ii)=boxname.zCM;
% 
%     end


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
    
    % offset
    
%     Wingbox_right.SCy=
    
    
    
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
    
 
    %% attachments 2  - engine_left
    
    Engine2=awi.model.BluffBody;
    Engine2.Name='Engine_left';
    
    % cylinder body
    Engine2.Radius=[1.4, 1.4, 1];
    Engine2.Eta =  [0, 0.6, 1];
    
%     Engine.Origin = [16.238147-3.5, 4.05 , -1.8];
    Engine2.Origin = [16.238147-3.5, -4.05 , 0.35432909];
  
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
    engine2_mass.SOffset=0.4;
    engine2_mass.Mass=3681;
    Engine2.add(engine2_mass);
    
    % add pylon
    pylon2_mass=awi.model.PointMass;   
    pylon2_mass.SOffset=0.9;
    pylon2_mass.Mass=1239/2;
    Engine2.add(pylon2_mass);

    build(Engine2)
      
    Wingbox_left.add(Engine2)
    
    %Control surfaces - flaps
    flap_L=awi.model.ControlSurface;
    flap_L.Eta=[0, 0.24];
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
%     
%     draw(Wingbox_right)
    

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
    Tailwing_left.XOffset=35;
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
    
    
    Body.add(Wingbox_right)
    Body.add(Wingbox_left)
    Body.add(Tailwing_right)
    Body.add(Tailwing_left)
    Body.add(Verticalwing)
    


    %The analysis methods require an 'awi.model.Aircraft' object
    % This is because some information is only known at the aircraft level,
    % e.g. all-up mass, reference span, area, etc.
    % Aircraft = awi.model.Aircraft;
    % Aircraft.add(LS);

    Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea+Wingbox_left.SurfaceArea]);
    Aircraft.RefSpan  = Wingbox_right.Span;
    Aircraft.RefChord = Wingbox_right.RootChord*0.697; %*0.697;
%     Aircraft.RefChord = Aircraft.RefArea/Aircraft.RefSpan; 


    %% Generate the loadcase object
    TrimLoadcase = awi.model.LoadCase;

    acMass = 500;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;

    TrimLoadcase.Name = 'A320_half_model_SOL144_1';
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
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase.CsDeflection=flap_angle*pi/180;
    
    
    build(TrimLoadcase)
    
    
    %% Trim loadcase 2
    TrimLoadcase2 = awi.model.LoadCase;

    acMass = 500;
    altitude          = 16000;
    mach_number       = 0.6;
    aircraft_velocity = 0.78*340;
    flap_angle=20;

    TrimLoadcase.Name = 'A320_half_model_SOL144_2';
    % TrimLoadcase.LoadCaseTypes = 'Static';(read only)
    % TrimLoadcase.CsDeflecTypes='fixed';(read only)
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

    %% Generate the FEM 

    run_folder = [
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\gust_debug']; %[-], folder for exporting the NASTRAN model

%     run_folder = [
%         'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\NSM_test']; %[-], folder for exporting the NASTRAN model


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
     
     
     
    trimFile = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase, MassCases,run_folder);
    
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\A320_half_model_SOL144*.*')
   
    delete(strcat(run_folder, '\A320_half_model_SOL144*.xdb'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.h5'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.log'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.f06'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.f04'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.xdb.*'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.h5.*'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.log.*'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.f06.*'));
    delete(strcat(run_folder, '\A320_half_model_SOL144*.f04.*'));
      
    NastranMethods1.runNastran(trimFile);
    
    draw(Aircraft)
% 
%     % result plotting for SOL 144
%     
%     F061=NastranMethods1.extractNastranResults(strcat(run_folder,'\A320_half_model_SOL144.f06'),'ReadF06',true,'ReadHDF5',false);
% 
%     % extract forces
%     M_P1=[F061.f06data.Bendingmoment.UMPLN1(69:89),F061.f06data.Bendingmoment.LMPLN1(89)]; % in-plane moment
%     M_P2=[F061.f06data.Bendingmoment.UMPLN2(69:89),F061.f06data.Bendingmoment.LMPLN2(89)]; % out of plane moment
%     T=[F061.f06data.Bendingmoment.UTORQUE1(69:89),F061.f06data.Bendingmoment.LTORQUE1(89)];% torque
% 
%     S_P1=[F061.f06data.Bendingmoment.USPLN1(69:89),F061.f06data.Bendingmoment.LSPLN1(89)]; % in plane shear
%     S_P2=[F061.f06data.Bendingmoment.USPLN2(69:89),F061.f06data.Bendingmoment.LSPLN2(89)]; % out of plane shear
% 
% 
%     Grid_coord = h5read(strcat(run_folder,'\A320_half_model_SOL144.h5'),'/NASTRAN/INPUT/NODE/GRID');
%     Displacement = h5read(strcat(run_folder,'\A320_half_model_SOL144.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');
%     
%     Y=Grid_coord.X(2,346:367);
%     Displacement_Z=Displacement.Z(346:367);
% 
%     figure % bending moment
%     plot(Y,M_P2,'b-s','MarkerFaceColor','b')
% 
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Bending moment (Nm)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
% 
%     figure %shear force
%     plot(Y, S_P2,'b-s','MarkerFaceColor','b')
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Vertical shear force (N)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
% 
%     figure % torque
%     plot(Y, T,'b-s','MarkerFaceColor','b')
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Torque (Nm)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
%     
%     figure 
%     plot(Y,Displacement_Z,'b-s','MarkerFaceColor','b')
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Deflection (m)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')


%% Run Sol 103


    
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_half;
    MassCases=awi.model.MassCases.empty;
    RefGrid=NastranMethods1.RefNode;

    fid = fopen(strcat(run_folder,'\NastranHeaderFile.dat'),'r');
    
    cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true );
    cac=cac{1};
    fclose( fid );
    
    mid_line0=22;
    mid_line=length(cac)-1; 
    fid = fopen( strcat(run_folder,'\full_model_103.dat'), 'w' );
    
    % write existing head file
    for jj = 1 : mid_line0
        fprintf(fid,'%s\n',cac{jj})       
    end
    
%     fprintf(fid,'SUPORT = %i\r\n',201)
    fprintf(fid,'SPC = %i\r\n',202)
    
    for ii =  mid_line0 : mid_line 
        fprintf(fid,'%s\n',cac{ii})       
    end
    
    
    %write spc part
    line='$.1.....2.......3.......4.......5.......6.......7.......8.......9.......10......\r\n';
    fprintf(fid,line);
    ID0=250;
    dof = 35;
    spc = 246;
    SPC_id = ID0;
    
    TrimData.RefGrid = RefGrid;
    TrimData.SPC_id  = SPC_id;
    TrimData.SPCdof  = spc;
    TrimData.SUPdof  = dof;
                        
%     spc_format='%-8s%-8i%-8i%-8i\r\n';
%     fprintf(fid,spc_format,'SPC1',1,123456,1005);
    
    fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPCADD', TrimData.SPC_id+2, TrimData.SPC_id, TrimData.SPC_id+1);
    %   - SPC
    rows=floor((numel(TrimData.RefGrid)-3)/5);
    remd=rem((numel(TrimData.RefGrid)-3),5);
    fprintf(fid, '%-8s%-8i%-8i%-8i%-8i%-8i\r\n', 'SPC1', TrimData.SPC_id, TrimData.SPCdof, TrimData.RefGrid(1:3).GID);
    fprintf(fid, '        %-8i%-8i%-8i%-8i%-8i\r\n', TrimData.RefGrid(4:4+rows*5-1).GID);
    format_last=strcat(repmat('%-8i',1,remd+1),'\r\n');
    fprintf(fid, format_last, '', TrimData.RefGrid(4+rows*5:end).GID);
    
    %charles added SPC 2
    % -SPC 2 at COG
    fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', TrimData.SPC_id+1, 1246, TrimData.RefGrid(34).GID);
    fprintf(fid, '%-8s%-8i%-8i\r\n', 'SUPORT', TrimData.RefGrid(34).ID, TrimData.SUPdof);
    
    % write the end 
    fprintf(fid,'%s\n',cac{end});


    %% Run the analysis - SOL 145 flutter 
    
        % NastranExe: 'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe'

%     NastranMethods1 = awi.methods.Nastran;
%     NastranMethods1.AnalysisModel = FEM_half;
%     MassCases=awi.model.MassCases.empty;
%     
%     FlightPoint=awi.model.FlightPoint;
%     
%     FlightPoint.Mach=0.78;
% %     FlightPoint.AcvELOCITY=50;
%     FlightPoint.Altitude = 36000;
%     
%     getFlightPointData(FlightPoint)
%     
%     flutterFile = NastranMethods1.writeFlutterFile(Aircraft, FlightPoint, MassCases, '1246', run_folder, 'RequestModeshapes',true,'FlutterMethod','pk');
%     
%     delete(strcat(run_folder,'\flutter_analysis*.xdb'));
%     delete(strcat(run_folder,'\flutter_analysis*.h5'));
%     delete(strcat(run_folder,'\flutter_analysis*.log'));
%     delete(strcat(run_folder,'\flutter_analysis*.f06'));
%     delete(strcat(run_folder,'\flutter_analysis*.f04'));
%     
%     NastranMethods1.runNastran(flutterFile);
% 
% 
% % flutter results Vg Vf;
% 
% flutter_data = h5read(strcat(run_folder,'\flutter_analysis.h5'),'/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');
% 
% modes=[1:35];
% 
% figure
% 
% for i=1:6
%     Modes_pt=modes(i);
%     [index1,~]=find(flutter_data.POINT==Modes_pt);
%     
%     velocity=flutter_data.VELOCITY(index1);
%     frequency=flutter_data.FREQUENCY(index1);
%     
%     
%     plot(velocity,frequency,'s-','LineWidth',1)
%        
%     set(gcf,'Color','w')
%     xlabel('Velocity','Interpreter','latex','FontSize',12)
%     ylabel('Frequency','Interpreter','latex','FontSize',12)
%     hold on
%     axis([260 400 0 25])
%     
% end
% 
% legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','Interpreter','latex','FontSize',10)
%  
%  
% figure
% for i=1:6
%     Modes_pt=modes(i);
%     [index1,~]=find(flutter_data.POINT==Modes_pt);
%     
%     velocity=flutter_data.VELOCITY(index1);
%     
%     damping=flutter_data.DAMPING(index1);
%     
%     plot(velocity,damping,'s-','LineWidth',1)
%     set(gcf,'Color','w')
%     xlabel('Velocity','Interpreter','latex','FontSize',12)
%     ylabel('Damping','Interpreter','latex','FontSize',12)
%     hold on
%     
% end    
% legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','Interpreter','latex','FontSize',10)    
    
   

%% Run gust analysis SOL 146

GustLoadcase = awi.model.LoadCase;

GustLoadcase.Altitude   = 36000;
GustLoadcase.AcVelocity = 0.78*340;
GustLoadcase.AcMass = 500;
GustLoadcase.Mach = 0.78;
GustLoadcase.GustLength = [20:20:200];

FlightPoint=awi.model.FlightPoint;
FlightPoint.Mach=0.78;
% FlightPoint.AcVelocity=FlightPoint.Mach*340;
FlightPoint.Altitude = 36000;
getFlightPointData(FlightPoint,'ISA');

NastranMethods1 = awi.methods.Nastran;
NastranMethods1.AnalysisModel = FEM_half;
MassCases=awi.model.MassCases.empty;


gustfile=NastranMethods1.writeGustFile(Aircraft, GustLoadcase, MassCases, FlightPoint, run_folder,'DatFilename','gust_analysis');
NastranMethods1.runNastran(gustfile);

%% Result 
    
    % Bending moment 
    Gust_res=Gust_extract(strcat(run_folder,'\gust_analysis.f06')); 
    
    NumGust=numel(GustLoadcase.GustLength);
    t=linspace(0,2,201);
    
    M_root=zeros(NumGust,201);

    peaks=zeros(NumGust,2);
    
    for j=1:NumGust
        
        for i=1:201
%             f_name=fieldnames(F);
            
            M_root(j,i)=Gust_res.Bendingmoment(i+(j-1)*201).UMPLN2(108); % out of plane moment 69 at root
            
            
        end
        
        peaks(j,:)=[max(M_root(j,:)),min(M_root(j,:))];
        
    end
    
figure 
plot(t,M_root(4,:),'b-')
% hold on 
% plot(t,M_root(2,:),'r-')
% hold on 
% plot(t,M_root(3,:),'r-')
% hold on 
% plot(t,M_root(4,:),'r-')
% hold on 
% plot(t,M_root(5,:),'k.')
% hold on 
% plot(t,M_root(6,:),'k.')
% hold on 
% plot(t,M_root(7,:),'k.')
% hold on 
% plot(t,M_root(8,:),'k.')

xlabel('Time (s)','FontSize',12,'Interpreter','latex')
ylabel('Root Bending moment (Nm)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

%  peak plot 
figure 
plot(GustLoadcase.GustLength,peaks(:,1),'b.')
hold on 
plot(GustLoadcase.GustLength,peaks(:,2),'r.')

xlabel('Gust length (m)','FontSize',12,'Interpreter','latex')
ylabel('Root Bending moment (Nm)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

%  z-motion of cog and wing tip

Gust_data = h5read(strcat(run_folder,'\gust_analysis.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');
% Z_disp_data=Gust_data.Z;

[ind1,~]=find(Gust_data.ID==1366);
[ind2,~]=find(Gust_data.ID==1387);

temp1=Gust_data.Z([ind1]);
temp2=Gust_data.Z([ind2]);

Z_disp_cog=reshape(temp1,[201,NumGust]);
Z_disp_tip=reshape(temp2,[201,NumGust]);



figure 
plot(t,Z_disp_cog(:,1))
hold on 
plot(t,Z_disp_cog(:,2))
hold on 
plot(t,Z_disp_cog(:,3))
hold on 
plot(t,Z_disp_cog(:,4))
hold on 
plot(t,Z_disp_cog(:,5))
hold on 
plot(t,Z_disp_cog(:,6))
hold on 
plot(t,Z_disp_cog(:,7))
hold on 
plot(t,Z_disp_cog(:,8))


figure 
plot(t,Z_disp_tip(:,1))
hold on 
plot(t,Z_disp_tip(:,2))
hold on 
plot(t,Z_disp_tip(:,3))
hold on 
plot(t,Z_disp_tip(:,4))
hold on 
plot(t,Z_disp_tip(:,5))
hold on 
plot(t,Z_disp_tip(:,6))
hold on 
plot(t,Z_disp_tip(:,7))
hold on 
plot(t,Z_disp_tip(:,8))
hold on 
plot(t,Z_disp_tip(:,9))
hold on 
plot(t,Z_disp_tip(:,10))

figure 
plot(t,Z_disp_tip(:,1)-Z_disp_cog(:,1))
hold on 
plot(t,Z_disp_tip(:,2)-Z_disp_cog(:,2))
hold on 
plot(t,Z_disp_tip(:,3)-Z_disp_cog(:,3))
hold on 
plot(t,Z_disp_tip(:,4)-Z_disp_cog(:,4))
hold on 
plot(t,Z_disp_tip(:,5)-Z_disp_cog(:,5))
hold on 
plot(t,Z_disp_tip(:,6)-Z_disp_cog(:,6))
hold on 
plot(t,Z_disp_tip(:,7)-Z_disp_cog(:,7))




%% Run the analysis - SOL 146 gust (old code)

   
%     GustLoadcase = awi.model.LoadCase;
% %     GustLoadcase.LoadCaseType = 'Pratt Gust';
%     GustLoadcase.Altitude   = 3000;
%     GustLoadcase.AcVelocity = 0.70*340;
%     GustLoadcase.AcMass = 70000;
%     GustLoadcase.Mach = 0.7;
%     GustLoadcase.GustLength = [50,200];
% %     obj=GustLoadcase;
% 
% % [a,b]=calculateGustLoadFactor(GustLoadcase,Aircraft);
% 
% 
% %     acMass = 500;
% %     altitude          = 36000; % crusing altitude feet
% %     mach_number       = 0.5;
% %     aircraft_velocity = mach_number*340;
% % 
% %     GustLoadcase.Name = 'A320_half_model_SOL146';
% %     GustLoadcase.Altitude   = altitude;
% %     GustLoadcase.Mach       = mach_number;
% %     GustLoadcase.AcVelocity = aircraft_velocity;
% %     GustLoadcase.AcMass = acMass;
% %     FP   = getFlightPointData(GustLoadcase);
%     
%     
%     FlightPoint=awi.model.FlightPoint;  
%     FlightPoint.Mach=0.78;
%     FlightPoint.AcVelocity=FlightPoint.Mach*340;
%     FlightPoint.Altitude = 36000;
%     FP=getFlightPointData(FlightPoint,'ISA');
%     
%     % gust parameters
% % for H=20:10:100    
%     % Gust gradient ranged from 9 - 107 m
%     H=100;
%     L=2*H;
%     Kg=1; % gust alleviation factor
%     Ng=64; %number of points 
%     
%     tmax=L/FP.AcVelocity;
%     ts=linspace(0,tmax,Ng);
%           
% %     CL_alfa=2*pi;
% %     GustLoadcase.calculateGustLoadFactor(Aircraft, CL_alfa);
%     
%     %Calculate gust velocity
%     Uref = [ ...
%         0    , 15000, 60000 ; ...
%         17.07, 13.41, 6.36 ];
%     Ug = interp1(Uref(1, :), Uref(2, :),FP.Altitude);  %[m/s], EAS
%     U_ref = Ug ./ sqrt(FP.DensityRatio);  %[m/s], TAS
%     
%     Uds=U_ref*Kg*(H/107)^(1/6);
%     
%     gust_1mc=Uds*(1-cos(2*pi*FP.AcVelocity*ts/L))*0.5;
%     
% %     gust_1mc_w=-Uds*(1-cos(2*pi*FP.AcVelocity*(ts-ts(end))/L))/2;
%     
%     figure 
%     
% %     plot(ts+ts(end),gust_1mc_w,'b.')
% %     
% %     hold on 
%     
%     plot(ts*FP.AcVelocity,gust_1mc,'b.')
%     hold on
%     set(gcf,'Color','w')
%     xlabel('Guest length (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Gust Velocity (m/s)','FontSize',12,'Interpreter','latex')
%     
%     
% % end
% 
%     
%     
%     bulk_data=zeros(1,numel(ts)*2);
%     
%     for i=0:numel(ts)-1
%         
%         index1=1+2*i;
%         index2=2+2*i;
%         bulk_data(index1)=ts(i+1);
%         bulk_data(index2)=gust_1mc(i+1)/max(gust_1mc);
%         
%     end
%     
%     %% data entry summery 
%     % enter TLOAD data
%     gust_data.TLid=100;
%     
%     % TABLED1 ID
%     gust_data.TABLED1.id=1;
%     
%     % enter DLOAD data
%     gust_data.DLid=10;
%     gust_data.DLtab=1;
%     gust_data.DLS=1.0; % scale factor 1
%     gust_data.DLSi=1.0; % scale factor 2
%     
%     % enter DAREA data
%     gust_data.DAREA.sid=200;
%     gust_data.DAREA.p1=NastranMethods1.RefNode(34).GID; % choose the wing root Grid ID
%     gust_data.DAREA.c1=1;
%     gust_data.DAREA.a1=1;
%     
%     % enter GUST data
%     gust_data.GUSTid=300;
%     gust_Data.Wg=Uds/FP.AcVelocity;
%     gust_data.X0=-100; %initial position of the ac
%     
%     % FREQUQ
%     gust_data.FREQ.id=600;
%     gust_data.FREQ.F1=0;
%     gust_data.FREQ.DF=0.05;
%     gust_data.FREQ.NDF=600;
%     
%     % TSTEP
%     gust_data.TSTEP.id=700;
%     gust_data.TSTEP.N=200;
%     gust_data.TSTEP.DT=0.01;
%     gust_data.TSTEP.NO1=1;
%     
%     
%     % TABDAMP1
%     
%     gust_data.TABDMP1.TID = 1000;
%     gust_data.TABDMP1.x   = [0   , 1000];
%     gust_data.TABDMP1.y   = [0.02, 0.02];
%     
%     dat = [gust_data.TABDMP1.x ; gust_data.TABDMP1.y];
%     str = cellstr(num2str(dat(:), '%#-8.2g'));
%     str = [str ; {'ENDT'}];
%     if numel(str) > 8
%         error('Update code for writing structural damping table');
%     end
%     
%     % - EIGR
%     gust_data.EIGRL.SID = 20;
%     gust_data.EIGRL.V0  = 0;
%     gust_data.EIGRL.V1  = [];
%     gust_data.EIGRL.ND  = 30;
%     
%     
%     %---------------------------------------------------------------
%       
% %     plot(ts,gust_1mc,'b.')
% 
%      NastranMethods1 = awi.methods.Nastran;
%      NastranMethods1.AnalysisModel = FEM_half;
%      MassCases=awi.model.MassCases.empty;
%      ID0=200;   
% % 
% %     NastranMethods1.writeGustFile(Aircraft, GustLoadcase, MassCases, FlightPoint, run_folder)
% 
% % write Gust file 
% 
%    fid = fopen( strcat(run_folder,'\guest_analysis_corrected.dat'), 'w' );
%    
% % write headlines
% 
%  %Executive Control
%  
%  awi.fe.FEBaseClass.writeHeading(fid, 'E X E C U T I V E  C O N T R O L');
%  fprintf(fid, 'SOL %i\r\n', 146);
%  fprintf(fid, 'ECHOOFF         $ SUPPRESSES THE ECHO OF EXECUTIVE CONTROL\r\n');
%  fprintf(fid, 'CEND\r\n');
%  
%  %Case Control
% awi.fe.FEBaseClass.writeHeading(fid, 'C A S E  C O N T R O L');
% awi.fe.FEBaseClass.writeHeading(fid, 'O U T P U T  O P T I O N S');
%  fprintf(fid, 'LINE = 99999   $ SPECIFIES THE NUMBER OF LINES PER PRINTED PAGE\r\n');
%  fprintf(fid, 'ECHO = NONE    $ SUPPRESSES THE ECHO OF BULK DATA\r\n');
%  fprintf(fid, 'METHOD = %i\r\n', gust_data.EIGRL.SID );
%  %
% awi.fe.FEBaseClass.writeHeading(fid, 'O U T P U T  Q U A N T I T I E S');
%  fprintf(fid, 'DISP(SORT1)  = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
%  fprintf(fid, 'FORCE(SORT1) = ALL $ OUTPUT ALL ELEMENT FORCES\r\n');
%  fprintf(fid, 'STRESS = ALL $ OUTPUT ALL ELEMENT FORCES\r\n');
%  fprintf(fid, 'AEROF       = ALL $ OUTPUT ALL AERO FORCES\r\n');
%  fprintf(fid, 'APRESSURE   = ALL $ OUTPUT ALL AERO PRESSURES\r\n');
%  fprintf(fid, 'MONITOR  = ALL\r\n');
%  %
%  awi.fe.FEBaseClass.writeHeading(fid, 'G L O B A L  C A R D S')
%  
%  fprintf(fid, 'SPC  = %i\r\n', 1);
%  
%  awi.fe.FEBaseClass.writeSubHeading(fid, 'S U B C A S E S');
%  fprintf(fid, 'SUBCASE %i\r\n', 1);
%  fprintf(fid, 'LABEL = %s\r\n','Gust response');
%  fprintf(fid, 'SDAMP  = %i\r\n', gust_data.TABDMP1.TID);
%  fprintf(fid, 'FREQ   = %i\r\n', gust_data.FREQ.id);
%  fprintf(fid, 'TSTEP  = %i\r\n', gust_data.TSTEP.id);
%  fprintf(fid, 'GUST  = %i\r\n', gust_data.GUSTid);
% %  fprintf(fid, 'DLOAD  = %i\r\n', gust_data.DLid);
% 
%  NastranMethods1.defaultBulkStatement(fid);
%  fprintf(fid, 'PARAM,MACH,%8.2f\r\n',FP.Mach);
%  fprintf(fid, 'PARAM,Q,%8.2f\r\n',FP.DynPressure);
%  fprintf(fid, 'PARAM,GUSTAERO,-1\r\n');
%  
%  awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
%  
%  %   - SPC symmetric BC
%  fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPCADD', 1, 10, 20);
%  
%  SCP_nodes=NastranMethods1.RefNode;
%  
%  % constraint along fuselage  -----------------------------------------------------
%  
%  rows=floor((numel( SCP_nodes)-3)/5);
%  remd=rem((numel( SCP_nodes)-3),5);
%  fprintf(fid, '%-8s%-8i%-8i%-8i%-8i%-8i\r\n', 'SPC1', 10, 246,  SCP_nodes(1:3).GID);
%  fprintf(fid, '        %-8i%-8i%-8i%-8i%-8i\r\n', SCP_nodes(4:4+rows*5-1).GID);
%  format_last=strcat(repmat('%-8i',1,remd+1),'\r\n');
%  fprintf(fid, format_last, '', SCP_nodes(4+rows*5:end).GID);
% 
% % constraint at cog  -----------------------------------------------------
% 
%  fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', 20, 1,  NastranMethods1.RefNode(34).GID);
%  
%  
%  
% %  fprintf(fid, '%-8s%-8i%-8s\r\n', 'SUPORT', NastranMethods1.RefNode(34).GID, '35');
%  
%  awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
%     
%  % frequency 
%  fprintf(fid, '%-8s%-8i%-8.3f%-8.3f%-8i\r\n', 'FREQ1',gust_data.FREQ.id,gust_data.FREQ.F1,gust_data.FREQ.DF,gust_data.FREQ.NDF);
%  
%   % TSTEP
%  fprintf(fid, '%-8s%-8i%-8i%-8.2f%-8i\r\n', 'TSTEP',gust_data.TSTEP.id,gust_data.TSTEP.N,gust_data.TSTEP.DT,gust_data.TSTEP.NO1);
%  
%  % TABDMP
%  fprintf(fid, ['%-8s%-8i%-8s\r\n%-8s', repmat('%-8s', [1, numel(str)]), '\r\n'], ...
%      'TABDMP1', gust_data.TABDMP1.TID, 'CRIT', blanks(8), str{:});
%  
%  % DAREA  
%  fprintf(fid, '%-8s%-8i%-8i%-8i%-8.2f\r\n', 'DAREA',gust_data.DAREA.sid,gust_data.DAREA.p1,gust_data.DAREA.c1,gust_data.DAREA.a1);
%            
%  
%  awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
%   
%  % - TLOAD1
%  fprintf(fid, ['%-8s%-8i%-8i',blanks(16),'%-8i\r\n'], 'TLOAD1',gust_data.TLid,gust_data.DAREA.sid,gust_data.TABLED1.id);
%  
%   % - DLOAD
% %   fprintf(fid, '%-8s%-8i%-8.3f%-8.3f%-8i\r\n', 'DLOAD',gust_data.DLid,gust_data.DLS,gust_data.DLSi,gust_data.TLid);
%    
%  % - Gust card
%  fprintf(fid, '%-8s%-8i%-8i%-8.3f%-8.1f%-8.1f\r\n', 'GUST',gust_data.GUSTid,gust_data.TLid,gust_Data.Wg,gust_data.X0,FP.AcVelocity);
%  
%   % - TABLED1
%   fprintf(fid, '%-8s%-8i\r\n', 'TABLED1',gust_data.TABLED1.id);
%   
%   fprintf(fid, '        %-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f\r\n',bulk_data);
%   
%   fprintf(fid, [blanks(8),'%-8s\r\n'],'ENDT');
%   
%   awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
%  
%   
%   % - EIGR 
% %   fprintf(fid, '%-8s%-8i%-#8.3g%-8s%-8i\r\n', 'EIGRL', gust_data.EIGRL.SID, 0, blanks(8),gust_data.EIGRL.ND);
%  fprintf(fid, '%-8s%-8i%-8s%-#8.3g%-16s%-8i\r\n', 'EIGR', gust_data.EIGRL.SID, 'MGIV', 0, blanks(16),gust_data.EIGRL.ND);
%  
%   % - AERO
%   NastranMethods1.writeUnsteadyAeroEntry(fid, Aircraft, FP(1));
%   
%   
%   % - MKAERO1
%   gust_data.MKAERO.M = unique([FP.Mach]);
%   gust_data.MKAERO.K = [0.001, 0.005, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6, 2, 2.5, 3, 3.5];
%   awi.methods.Nastran.writeMKAERO1(fid, gust_data.MKAERO);
%   
%   
%   % Including files
%   [~, includeFiles] = export(NastranMethods1.AnalysisModel, run_folder, 'WriteHeaderFile', false);
%   awi.fe.FEBaseClass.writeSubHeading(fid, 'I N C L U D E  F I L E S');
%   awi.fe.FEBaseClass.writeIncludeStatement(fid, includeFiles);
%   
%   %End of file
%   fprintf(fid, 'ENDDATA\r\n');
%   
%   %Close the file
%   fclose(fid);
%             
%   NastranMethods1.runNastran(strcat(run_folder,'\guest_analysis_corrected.dat'));          
%    
% % NastranMethods1 = awi.methods.Nastran;
% % NastranMethods1.AnalysisModel = FEM_half;
% 
% 
% 
% 
% % F.H10=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH10.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H20=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH20.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H30=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH30.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H40=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH40.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H50=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH50.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H60=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH60.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H70=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH70.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H80=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH80.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H90=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH90.f06'),'ReadF06',true,'ReadHDF5',false);
% % F.H100=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysisH100.f06'),'ReadF06',true,'ReadHDF5',false);
% 
% F06_nospc=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysis_no_spc1.f06'),'ReadF06',true,'ReadHDF5',false);
% F06_sg_spc=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysis_test.f06'),'ReadF06',true,'ReadHDF5',false);
% F06_spc=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysish100.f06'),'ReadF06',true,'ReadHDF5',false);
% F06_fixed_cog=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysis_fixed_cog.f06'),'ReadF06',true,'ReadHDF5',false);
% F06_spc_no_suport=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysis_spc_no_suport.f06'),'ReadF06',true,'ReadHDF5',false);
% F06_BC_config1=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysis_new_BC.f06'),'ReadF06',true,'ReadHDF5',false); % fuselage 246 cog 1 no suport
% F06_BC_config1_changed_ref_c=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysis_new_BC_changed_refchord.f06'),'ReadF06',true,'ReadHDF5',false); % fuselage 246 cog 1 no suport
% F06_no_cs=NastranMethods1.extractNastranResults(strcat(run_folder,'\guest_analysis_no_cs.f06'),'ReadF06',true,'ReadHDF5',false); % fuselage 246 cog 1 no suport
% F06_corrected=NastranMethods1.extractNastranResults(strcat(run_folder,'\Loadcase.f06'),'ReadF06',true,'ReadHDF5',false); % fuselage 246 cog 1 no suport
% Force= Gust_extract(strcat(run_folder,'\Loadcase.f06'));
% t=linspace(0,2,201);
% M_root1=zeros(1,201);
% M_root2=zeros(1,201);
% M_root3=zeros(1,201);
% M_root4=zeros(1,201);
% M_root5=zeros(1,201);
% M_root6=zeros(1,201);
% M_root7=zeros(1,201);
% M_root8=zeros(1,201);
% M_root9=zeros(1,201);
% M_root10=zeros(1,201);
% M_root11=zeros(1,201);
% S_root9=zeros(1,201);
% % peaks=zeros(10,2);
% % 
% % for j=1:10
% %     
% %     for i=1:201
% %         f_name=fieldnames(F);
% %          
% %         M_root(j,i)=F.(f_name{j}).f06data.Bendingmoment(i).UMPLN2(70); % out of plane moment 69 at root
% %        
% %     
% %     end
% %     
% %     peaks(j,:)=[max(M_root(j,:)),min(M_root(j,:))];
% %     
% % end
% 
% % for i=1:201
% %     M_rand(i)=F10.f06data.Bendingmoment(i).UMPLN2(75); % out of plane moment 69 at root
% %     
% % end
% 
%     for i=1:201
%       
%          
% %         M_root1(i)=F06_nospc.f06data.Bendingmoment(i).UMPLN2(71); % out of plane moment 69 at root
% %         M_root2(i)=F06_sg_spc.f06data.Bendingmoment(i).UMPLN2(71); % out of plane moment 69 at root
% %         M_root3(i)=F06_spc.f06data.Bendingmoment(i).UMPLN2(71); % out of plane moment 69 at root
% %         M_root4(i)=F06_fixed_cog.f06data.Bendingmoment(i).UMPLN2(71); % out of plane moment 69 at root
% %         M_root5(i)=F06_spc_no_suport.f06data.Bendingmoment(i).UMPLN2(71); % out of plane moment 69 at root
% %         M_root6(i)=F06_BC_config1.f06data.Bendingmoment(i).UMPLN2(71); % out of plane moment 69 at root
% %         M_root7(i)=F06_BC_config1_changed_ref_c.f06data.Bendingmoment(i).UMPLN2(71); % out of plane moment 69 at root
% %         
% % %         M_root8(i)=F06_BC_config1_changed_ref_c.f06data.Bendingmoment(i).UMPLN2(90); % out of plane moment 69 at root
% %         M_root9(i)=F06_no_cs.f06data.Bendingmoment(i).UMPLN2(69); % out of plane moment 69 at root
% 
% %         M_root10(i)=F06_corrected.f06data.Bendingmoment(i+201).UMPLN2(69); % out of plane moment 69 at root
% %         S_root9(i)=F06_no_cs.f06data.Bendingmoment(i).USPLN2(70); % out of plane moment 69 at root
%         M_root11(i)=F06_corrected.f06data.Bendingmoment(i).UMPLN2(69);
%         
%         
%         
%         
%     
%     end
% 
%     
% 
% figure 
% % plot(t,M_root1,'bs')
% % hold on 
% % plot(t,M_root2,'r.')
% % hold on 
% % plot(t,M_root3,'k.')
% % hold on 
% % plot(t,M_root4,'m.')
% % hold on 
% % plot(t,M_root5,'g.')
% % hold on 
% % plot(t,M_root6,'gs')
% % % hold on 
% % plot(t,M_root7,'b.')
% % hold on 
% % plot(t,M_root8,'rs')
% % hold on 
% plot(t,M_root11,'r.')
% % hold on 
% % plot(t,M_root11,'b.')
% % plot(t,S_root9,'b.')
% 
% % plot(t,M_root(2,:),'r.')
% % hold on 
% % plot(t,M_root(3,:),'k.')
% % hold on 
% % plot(t,M_root(4,:),'m.')
% % hold on 
% % plot(t,M_root(5,:),'g.-')
% % hold on 
% % plot(t,M_root(6,:),'b-')
% % hold on 
% % plot(t,M_root(7,:),'r-')
% % hold on 
% % plot(t,M_root(8,:),'k-')
% % hold on 
% % plot(t,M_root(9,:),'m-')
% % hold on 
% % plot(t,M_root(10,:),'g-')
% 
% xlabel('Time (s)','FontSize',12,'Interpreter','latex')
% ylabel('Root Bending moment (Nm)','FontSize',12,'Interpreter','latex')
% set(gcf,'color','w')
% 
% legend('No SPC','SPC on CoG','SPC on Fuselage','Fixed CoG','SPC on cog no suport','BC config1','BC config1_changed_c')
% % 
% % figure 
% % Hx=[10:10:100]*2;
% % 
% % plot(Hx',peaks(:,1),'rs-')
% % hold on 
% % plot(Hx',peaks(:,2),'bs-')
% % xlabel('Gust length (m)','FontSize',12,'Interpreter','latex')
% % ylabel('Root Bending moment (Nm)','FontSize',12,'Interpreter','latex')
% % set(gcf,'color','w')
% % axis([0 200 -2e7 2e7])
% %  % - EPOINT
%  
%  
%  
% %  % Result extract
% %  F061=NastranMethods1.extractNastranResults('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\gust_test\guest_analysis.f06','ReadF06',true,'ReadHDF5',false);
% %  
% %  % extract forces
% %  M_P1=[F061.f06data.Bendingmoment.UMPLN1(69:89),F061.f06data.Bendingmoment.LMPLN1(89)]; % in-plane moment
% %  M_P2=[F061.f06data.Bendingmoment.UMPLN2(69:89),F061.f06data.Bendingmoment.LMPLN2(89)]; % out of plane moment
% %  T=[F061.f06data.Bendingmoment.UTORQUE1(69:89),F061.f06data.Bendingmoment.LTORQUE1(89)];% torque
% %  
% %  S_P1=[F061.f06data.Bendingmoment.USPLN1(69:89),F061.f06data.Bendingmoment.LSPLN1(89)]; % in plane shear
% %  S_P2=[F061.f06data.Bendingmoment.USPLN2(69:89),F061.f06data.Bendingmoment.LSPLN2(89)]; % out of plane shear
% %  
% %  
% %  Grid_coord = h5read('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\NSM_test\A320_half_model_SOL144.h5','/NASTRAN/INPUT/NODE/GRID');
% %  Displacement = h5read('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\NSM_test\A320_half_model_SOL144.h5','/NASTRAN/RESULT/NODAL/DISPLACEMENT');
% %  
% %  Y=Grid_coord.X(2,346:367);
% %  Displacement_Z=Displacement.Z(346:367);
% %  
% %  figure % bending moment
% %  plot(Y,M_P2,'b-s','MarkerFaceColor','b')
% %  
% %  xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
% %  ylabel('Bending moment (Nm)','FontSize',12,'Interpreter','latex')
% %  set(gcf,'color','w')
% %  
% %  figure %shear force
% %  plot(Y, S_P2,'b-s','MarkerFaceColor','b')
% %  xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
% %  ylabel('Vertical shear force (N)','FontSize',12,'Interpreter','latex')
% %  set(gcf,'color','w')
% %  
% %  figure % torque
% %  plot(Y, T,'b-s','MarkerFaceColor','b')
% %  xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
% %  ylabel('Torque (Nm)','FontSize',12,'Interpreter','latex')
% %  set(gcf,'color','w')
%  
%  
%  
%  
%  %% test code
%  
%     H=[50,100];
%     N=numel(H);
%     Ng=64; %number of points
%     bulk_data=zeros(N,Ng*2);
%     
%     for n=1:N
%         L=2*H(n);
%         Kg=1; % gust alleviation factor
%        
%         
%         tmax=L/FP.AcVelocity;
%         ts=linspace(0,tmax,Ng);
%         
%         
%         %Calculate gust velocity
%         Uref = [ ...
%             0    , 15000, 60000 ; ...
%             17.07, 13.41, 6.36 ];
%         Ug = interp1(Uref(1, :), Uref(2, :),FP.Altitude);  %[m/s], EAS
%         U_ref = Ug ./ sqrt(FP.DensityRatio);  %[m/s], TAS
%         
%         Uds=U_ref*Kg*(H(n)/107)^(1/6);
%         
%         gust_1mc=Uds*(1-cos(2*pi*FP.AcVelocity*ts/L))*0.5;
%         
%         plot(ts*FP.AcVelocity,gust_1mc,'b.')
%         hold on
%         set(gcf,'Color','w')
%         xlabel('Guest length (m)','FontSize',12,'Interpreter','latex')
%         ylabel('Gust Velocity (m/s)','FontSize',12,'Interpreter','latex')
%         
% %         bulk_data=zeros(n,numel(ts)*2);
%         
%         for i=0:numel(ts)-1
%             
%             index1=1+2*i;
%             index2=2+2*i;
%             bulk_data(n,index1)=ts(i+1);
%             bulk_data(n,index2)=gust_1mc(i+1)/max(gust_1mc);
%             
%         end
%         
%     end
%  
%  
%   figure
%         
%         %     plot(ts+ts(end),gust_1mc_w,'b.')
%         %
%         %     hold on
%         
%         plot(ts*FP.AcVelocity,gust_1mc,'b.')
%         hold on
%         set(gcf,'Color','w')
%         xlabel('Guest length (m)','FontSize',12,'Interpreter','latex')
%         ylabel('Gust Velocity (m/s)','FontSize',12,'Interpreter','latex')
%         
%         
%         % end 
%  
%  
%  
gust = importdata('wind_sample1.txt');
 
 x1=gust(:,1);
 y1=gust(:,2);
 
 x2=gust(:,3);
 y2=gust(:,4);
 
 x3=gust(:,5);
 y3=gust(:,6);
 
 x4=gust(:,7);
 y4=gust(:,8);
 
 figure 
 
 plot(x1,y1,'b.')
 hold on
 plot(x2,y2,'b.')
  hold on
 plot(x3,y3,'b.')
   hold on
 plot(x4,y4,'b.')
%  
%  
% %  pad = blanks(2);
% %  for iS = 1 : numel(TrimData.Trim)
% %      fprintf(fid, 'SUBCASE %i\r\n', iS);
% %      fprintf(fid, '%sTRIM = %i\r\n', pad, TrimData.Trim(iS).SID);
% %  end
% %    
%    
%   %Bulk Data
% %             obj.defaultBulkStatement(fid);
% %             %Write the cards
% %             awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
% %             %   - SPC
% %             fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', TrimData.SPC_id, TrimData.SPCdof, TrimData.RefGrid.ID);
% %             %   - SUPORT
% %             fprintf(fid, '%-8s%-8i%-8i\r\n', 'SUPORT', TrimData.RefGrid.ID, TrimData.SUPdof);
% %             %   - AEROS
% %             obj.writeSteadyAeroEntry(fid, Aircraft);
% %             %   - AESTAT
% %             data = [num2cell(TrimData.AESTAT_id) ; TrimData.AESTAT];
% %             str  = sprintf('AESTAT,%i,%s-', data{:});
% %             data = strsplit(str, '-');
% %             %             data = strcat(repmat({['AESTAT  ', blanks(8)]}, [numel(AESTAT), 1]), AESTAT');
% %             fprintf(fid, '%s\r\n', data{:});
% %             %   - TRIM
% %             obj.writeTrimEntries(fid, TrimData.Trim, LoadCases);
% %             fprintf(fid, 'PARAM,AUNITS,0.1020\r\n');
% %             
% %             %Additional bulk data files
% %             if ~iscell(includeFiles)
% %                 includeFiles = {includeFiles};
% %             end
% %             awi.fe.FEBaseClass.writeSubHeading(fid, 'I N C L U D E  F I L E S');
% %             awi.fe.FEBaseClass.writeIncludeStatement(fid, includeFiles);
% %             
% %             %End of file
% %             fprintf(fid, 'ENDDATA\r\n');
% %             
% %             %Close the file
% %             fclose(fid);      
     
     
     
     
    