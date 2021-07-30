 

%% Import sized wing model

Wing_Model=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A140\AR10\Res_AR10_Eta_60_Model.mat');

Param=Wing_Model.Param;

% update aeropanel size
Param.Wing.AeroPanel_AR=2;

% %% Model check
% 
% if isfield(Param,'FWT')
%     
%     [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v1(Param);
%     
%     
% else
%     
%     [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
%     
% end
% 
% draw(FEM_full);


%% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[Param.Layout.Wing_Position,0,0]; %15
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
%     Wingbox_right.Dihedral=[0,0];
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

        FWT_R = insertWingFold(Wingbox_right, 'FlareAngle', 12.5, 'FoldAngle', Param.FWT.Fold_angle,'EtaFold',Param.FWT.Fold_eta);
%         FWT_R.HingeStiffness = [1e14 1e14 1e14 1e14 Param.FWT.Hinge_Stiffness 1e14];
        FWT_R.HingeStiffness = [0 0 0 0 1 0];
        
        
        %Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
        FrontSpar_FWT = awi.model.Spar;
        FrontSpar_FWT.XLoc = [0.15, 0.15];
        FrontSpar_FWT.Eta  = [0   , 1];
        RearSpar_FWT = awi.model.Spar;
        RearSpar_FWT.XLoc = [0.65, 0.65];
        RearSpar_FWT.Eta  = [0   , 1];
        
        FWT_R.add([FrontSpar_FWT, RearSpar_FWT]);

        
        FWT_R.AeroPanelLength=[];
        FWT_R.NumAeroPanel=Param.Wing.AeroPanel_Number;
        FWT_R.AeroPanelAR=Param.Wing.AeroPanel_AR;
 
        FWT_eta_= 0:0.1:1;
        FWT_Bheight=interp1([0,1],[Param.FWT.Root_Height,Param.FWT.Tip_Height],FWT_eta_);
        FWT_Bwidth=interp1([0,1],0.5*[Param.FWT.Root_Chord,Param.FWT.Tip_Chord],FWT_eta_);
        
        FWT_Spar=Param.FWT.Thickness(1:11);
        FWT_Skin=Param.FWT.Thickness(12:22);
        FWT_Strg=Param.FWT.Thickness(23:33);

%         FWT_Spar=0.02*ones(1:11);
%         FWT_Skin=0.003*ones(1:11);
%         FWT_Strg=0.18e-4*ones(1:11);

        
        
        fwt_d_strg=sqrt(FWT_Strg/0.36);
        fwt_t_strg=0.12*fwt_d_strg;
        

        % stringer pitch
        strg_n=Param.Wing.StringerPitch;
        
        %intialise data array
        FWT_A_val=zeros(1,11);
        FWT_Ixx_val=zeros(1,11);
        FWT_Izz_val=zeros(1,11);
        FWT_J_val=zeros(1,11);


    for ii=1:11

        boxname=strcat('Box',string(ii));
        boxname=awi.model.BoxBeam;
        boxname.BoxType='SymmetricBox';
        boxname.Height=FWT_Bheight(ii);
        boxname.Width=FWT_Bwidth(ii);
        boxname.CoverThickness=FWT_Skin(ii);
        boxname.SparThickness=FWT_Spar(ii);

        FWT_NumStrg=floor(FWT_Bwidth(ii)/strg_n);

        ts=fwt_t_strg(ii);
        ds=fwt_d_strg(ii);
        hs=FWT_Bheight(ii);
        ws=FWT_Bwidth(ii);

        Istrg_xx=((ts*ds^3/12)+(ds*ts^3/12 + ds*ts*(ds/2)^2)*2 + 3*ds*ts*(hs/2-ds/2)^2)*FWT_NumStrg*2;
        Istrg_zz_=(ds*ts^3/12 + (ts*ds^3/12 + ts*ds*(ds/2)^2)*2);

        if mod(FWT_NumStrg,2)==0
            offset=0.12:strg_n:ws/2;
            Istrg_zz=(FWT_NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;

        elseif mod(FWT_NumStrg,2)==1
            offset=0:strg_n:ws/2;
            Istrg_zz=(FWT_NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;

        end

        getGeometricProps(boxname)
        
        FWT_A_val(ii)=boxname.Abb + 3*ts*ds*FWT_NumStrg*2;
        FWT_Ixx_val(ii)=boxname.Ixx+Istrg_xx;
        FWT_Izz_val(ii)=boxname.Izz+Istrg_zz;
        FWT_J_val(ii)=boxname.Jbb;
        
    end
    
    
    % FWT structural properties
    
    FWT_R.A   =  FWT_A_val;
    FWT_R.A_eta=FWT_eta_;
    
    FWT_R.I11 = FWT_Izz_val;
    FWT_R.I11_eta=FWT_eta_;
    
    FWT_R.I22 = FWT_Ixx_val;
    FWT_R.I22_eta = FWT_eta_;
    
    FWT_R.J   = FWT_J_val;
    FWT_R.J_eta= FWT_eta_;
    
    
    % FWT Jig Twist
%     FWT_R.Twist = Param.FWT.Jig_Twist;
%     FWT_R.Twist_eta = Param.FWT.Jig_Eta;

    FWT_R.Twist = [0,0];
    FWT_R.Twist_eta = [0,1];
    
    

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
    for i=1:1:3
        handle=strcat('PM_right','i');
        handle=awi.model.PointMass;
        handle.SOffset=0.2+i*0.2;
        handle.Mass=0.1;
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
    
    %%sizing variables ---------------------------------------------
    
    thickness1=Param.Wing.Thickness(1:25);
    thickness2=Param.Wing.Thickness(26:50);
    Astrg=Param.Wing.Thickness(51:75);
        
    d_strg=sqrt(Astrg/0.36);
    t_strg=0.12*d_strg;
    
    % -------------------------------------------------------------

    % etaS=linspace(0,Wingbox_right.Span,NumSec);

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
%      Bheight=interp1(RData/RData(end),0.79*[RootH,RootH,RootH],etaRS);

    % stringer pitch 
    strg_n=Param.Wing.StringerPitch;

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
        A_val(ii)=boxname.Abb + 3*ts*ds*NumStrg*2;
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
    
    Y_eta=etaRS*RData(end);
    
    Wingbox_right.A   =  A_val;
    Wingbox_right.A_eta=eta_;

    Wingbox_right.I11 = Izz_val;
    Wingbox_right.I11_eta=eta_;

    Wingbox_right.I22 = Ixx_val;
    Wingbox_right.I22_eta = eta_;

    Wingbox_right.J   = J_val;
    Wingbox_right.J_eta= eta_;


    % Jig shape
    
%     Wingbox_right.Twist = Param.Wing.Jig_Twist;
%     Wingbox_right.Twist_eta = Param.Wing.Jig_Eta;

    Wingbox_right.Twist = [0,0];
    Wingbox_right.Twist_eta = [0,1];

    
    build(Wingbox_right)
    
    if isfield(Param,'FWT')
        
        Box_dimensions.Inboard.Height = [Bheight,FWT_Bheight(2:end)];
        Box_dimensions.Inboard.Width = [Bwidth,FWT_Bwidth(2:end)];
        Box_dimensions.Inboard.Stringer_length = [d_strg, fwt_d_strg(2:end)];
        Box_dimensions.Inboard.Stringer_thickness = [t_strg,fwt_t_strg(2:end)];
        
        Box_CrossSec.Izz=[Izz_val,FWT_Izz_val(2:end)];
        Box_CrossSec.Ixx=[Ixx_val,FWT_Ixx_val(2:end)];
        Box_CrossSec.A=[A_val,FWT_A_val(2:end)];
        
    else
        
        Box_dimensions.Inboard.Height = Bheight;
        Box_dimensions.Inboard.Width = Bwidth;
        Box_dimensions.Inboard.Stringer_length = d_strg;
        Box_dimensions.Inboard.Stringer_thickness = t_strg;
        
        Box_CrossSec.Izz=Izz_val;
        Box_CrossSec.Ixx=Ixx_val;
        Box_CrossSec.A=A_val;
        
    end
    
    %% Mass definition
    

    wingmass_eta=eta_;
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
    
    
    
    Aircraft = awi.model.Aircraft;
    
    Aircraft.add(Wingbox_right);
    
    if isfield(Param,'FWT')
        
        Aircraft.RefArea = Wingbox_right.SurfaceArea + FWT_R.SurfaceArea;
        Aircraft.RefSpan  = Wingbox_right.Span + FWT_R.Span;
        
    else
        Aircraft.RefArea = Wingbox_right.SurfaceArea;
        Aircraft.RefSpan  = Wingbox_right.Span;

    end
    
    
    
    
    Aircraft.RefChord = Wingbox_right.RootChord; %mean aerodynamic chord = 0.697 for A321 wing;
    
    FEM_full=convertToFE(Aircraft);
    

    
    % run folder
    run_folder='C:\Git\A321_sizing\hg_codes\results\wing_flutter';
    
    % %Export it to a file
    export(FEM_full, run_folder);
    
    draw(FEM_full)
    
    %% Flutter loadcase
    
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_full;
    MassCases=awi.model.MassCases.empty;
    
    FlightPoint=awi.model.FlightPoint;
    
    FlightPoint.Mach=0.1;
    FlightPoint.AcVelocity=50;
    FlightPoint.Altitude = 3000;
    
    getFlightPointData(FlightPoint)
    
    flutterFile = NastranMethods1.writeFlutterFile(Aircraft, FlightPoint, MassCases, '123456', run_folder, 'RequestModeshapes',true,'FlutterMethod','pk');
    
    delete(strcat(run_folder, '\flutter_analysis.xdb'));
    delete(strcat(run_folder, '\flutter_analysis.h5'));
    delete(strcat(run_folder, '\flutter_analysis.log'));
    delete(strcat(run_folder, '\flutter_analysis.f06'));
    delete(strcat(run_folder, '\flutter_analysis.f04'));
    delete(strcat(run_folder, '\flutter_analysis.op4'));
    
    delete(strcat(run_folder, '\flutter_analysis.xdb.*'));
    delete(strcat(run_folder, '\flutter_analysis.h5.*'));
    delete(strcat(run_folder, '\flutter_analysis.log.*'));
    delete(strcat(run_folder, '\flutter_analysis.f06.*'));
    delete(strcat(run_folder, '\flutter_analysis.f04.*'));
    
    NastranMethods1.runNastran(flutterFile);
    
    
    
    
    %% flutter results Vg Vf;
    
    
    flutter_data = h5read(strcat(run_folder,'\flutter_analysis.h5'),'/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');
    
    modes=[1:35];
    
    colors=[0, 0, 1; ...
        0.6350, 0.0780, 0.1840;...
        0, 0.5, 0;...
        1 0 0;
        0.9290, 0.6940, 0.1250;...
        0 0 0;...
        0.8500, 0.3250, 0.0980
        0 1 0;...
        1 0 1];
    
    
    figure
    
    for i=[1:9]
        Modes_pt=modes(i);
        [index1,~]=find(flutter_data.POINT==Modes_pt);
        
        velocity=flutter_data.VELOCITY(index1);
        frequency=flutter_data.FREQUENCY(index1);
        
        
        plot(velocity,frequency,'-','Color',[colors(i,:)],'LineWidth',1.5)
        
        set(gcf,'Color','w')
        xlabel('Velocity (m/s)','Interpreter','latex','FontSize',12)
        ylabel('Frequency (Hz)','Interpreter','latex','FontSize',12)
        hold on
        axis([150 600 0 30])
        %
    end
    grid on
    grid minor
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10,'Location','northwest')
    
    figure
    for i=[1:9]
        Modes_pt=modes(i);
        [index1,~]=find(flutter_data.POINT==Modes_pt);
        
        velocity=flutter_data.VELOCITY(index1);
        
        damping=flutter_data.DAMPING(index1);
        
        plot(velocity,damping,'-','Color',[colors(i,:)],'LineWidth',1.5)
        set(gcf,'Color','w')
        xlabel('Velocity (m/s)','Interpreter','latex','FontSize',12)
        ylabel('Damping','Interpreter','latex','FontSize',12)
        hold on
        axis([150 500 -1.5 0.5])
        
    end
    
    hold on
    plot([150,600],[0,0],'k--')
grid on
grid minor

    legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10,'Location','southwest')
    
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)


    
    %% Run the analysis- SOL 103
    % insert SPC
    fid = fopen(strcat(run_folder,'\NastranHeaderFile.dat'));
    
    cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true );
    cac=cac{1};
    fclose( fid )
    
    mid_line0=22;
    mid_line1=33;
    
    end_line=length(cac)-1;
    fid = fopen( strcat(run_folder,'\Modal_analysis.dat'), 'w' );
    
    % write existing head file
    for jj = 1 : mid_line0
        fprintf(fid,'%s\n',cac{jj})
    end
    
    fprintf(fid,'SPC = 1 \r\n')
    
    for ii =  mid_line0 : mid_line1
        fprintf(fid,'%s\n',cac{ii})
    end
    
    % correct eigen value extraction methods
    fprintf(fid, '%-8s%-8i%-8s%-#8.3g%-16s%-8i\r\n', 'EIGR', 117, 'MGIV', 0, blanks(16),30);
    
    
    for ii =  mid_line1+2 : end_line
        fprintf(fid,'%s\n',cac{ii})
    end
    
    
    %write spc part
    line='$.1.....2.......3.......4.......5.......6.......7.......8.......9.......10......\r\n';
    fprintf(fid,line);
    
    spc_format='%-8s%-8i%-8i%-8i\r\n';
    fprintf(fid,spc_format,'SPC1',1,123456,1008);
    
    % write the end
    fprintf(fid,'%s\n',cac{end})
    
    Sol103_file=strcat(run_folder,'\Modal_analysis.dat');
    
    
    delete(strcat(run_folder, '\Modal_analysis.xdb'));
    delete(strcat(run_folder, '\Modal_analysis.h5'));
    delete(strcat(run_folder, '\Modal_analysis.log'));
    delete(strcat(run_folder, '\Modal_analysis.f06'));
    delete(strcat(run_folder, '\Modal_analysis.f04'));
    delete(strcat(run_folder, '\Modal_analysis.op2'));
    
    delete(strcat(run_folder, '\Modal_analysis.xdb.*'));
    delete(strcat(run_folder, '\Modal_analysis.h5.*'));
    delete(strcat(run_folder, '\Modal_analysis.log.*'));
    delete(strcat(run_folder, '\Modal_analysis.f06.*'));
    delete(strcat(run_folder, '\Modal_analysis.f04.*'));
    delete(strcat(run_folder, '\Modal_analysis.op2.*'));
    
    NastranMethods1.runNastran(Sol103_file);
    
    
%     %% Result plot
%     % load the model
%     model = mni.import_matran(fullfile(run_folder,'modal_analysis.dat'));
%     model.draw
%      
%     % get modal data
%     res_modeshape = mni.result.f06(fullfile(run_folder,'modal_analysis.f06')).read_modeshapes;
%     res_freq = mni.result.f06(fullfile(run_folder,'modal_analysis.f06')).read_modes;
%     %% apply deformation result
%     %pick which mode to plot (1->6)
%     modeshape_num = 6;
%     
%     % apply the modal deformation
%     [~,i] = ismember(model.GRID.GID,res_modeshape.GID(modeshape_num,:));
%     model.GRID.Deformation = [res_modeshape.T1(modeshape_num,i);...
%         res_modeshape.T2(modeshape_num,i);res_modeshape.T3(modeshape_num,i)];
%     
%     % animate the plot to show the mode shape. 'scale' is the scale of the
%     % displacement. 'Frequency' is the frequency.
%     %
%     % PRESS SPACE TO END THE ANIMATION
%     %
%     model.update()
%     model.animate('Frequency',2,'Scale',10)
% %     
    
    
    
    
    
    
    
    
%% Damping 

% damping=[0
% 0.02
% 0.04
% 0.06
% 0.08
% 0.1
% 0.12
% 0.14
% 0.16];
%     
% 
% Flutter=[381.5
% 392
% 400
% 406.5
% 412
% 418
% 423.5
% 428.7
% 435];
% 
% figure 
% 
% plot(damping,Flutter,'bs-','LineWidth',1)
% xlabel('Damping','interpreter','latex','FontSize',12)
% ylabel('Flutter speed (m/s)','interpreter','latex','FontSize',12)
% 
% set(gcf,'Color','w')
% set(gca,'TickLabelInterpreter','latex')


%% FWT size
    
%     eta=[0,10,20,30,40];
%     
%     Flutter_Speed=[315,310,345,380,490];
%     
%     figure 
%     plot(eta, Flutter_Speed,'rs-','LineWidth',1,'MarkerFaceColor','r')
%     xlabel('Folding wingtip size ($\%$)','interpreter','latex','FontSize',12)
%     ylabel('Flutter speed (m/s)','interpreter','latex','FontSize',12)
%     
%     set(gcf,'Color','w')
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'FontSize',12)
    %% result 101
% 
%     % load the model
%     model = mni.import_matran(fullfile(run_folder,'SOL101.dat'));
%     model.draw
%     % extract the displacment data
%     res_disp = mni.result.f06(fullfile(run_folder,'SOL101.f06')).read_disp;
% 
%     % apply deformation results
%     [~,i] = ismember(model.GRID.GID,res_disp.GP);
%     model.GRID.Deformation = [res_disp.dX(:,i);res_disp.dY(:,i);res_disp.dZ(:,i)];
% 
%     % update the figure (Scale is the normalised scale of displacements 
%     % to plot)
%     model.update('Scale',0.000001)
    

    
    
    
   %% plot thickness 
   
%    spar1=Param.Wing.Thickness(1:25);
%    
%    spar_fwt=Param.FWT.Thickness(1:11);
%    
%    spar=[  spar1,  spar_fwt(2:end)];
%    
%    figure 
%    plot(spar,'bs')


%% FFAST Model result


%     run_folder='D:\pdf_nastran\pdf_nastran\FFAST_Test_FWT\FFAST Test';
% 
%     flutter_data = h5read(strcat(run_folder,'\145.h5'),'/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');
% 
%     modes=[1:35];
% 
%     figure
%     
%     for i=1:20
%         Modes_pt=modes(i);
%         [index1,~]=find(flutter_data.POINT==Modes_pt);
%         
%         velocity=flutter_data.VELOCITY(index1);
%         frequency=flutter_data.FREQUENCY(index1);
%         
%         
%         plot(velocity,frequency,'-','LineWidth',1)
%         
%         set(gcf,'Color','w')
%         xlabel('Velocity','Interpreter','latex','FontSize',12)
%         ylabel('Frequency','Interpreter','latex','FontSize',12)
%         hold on
%         %     axis([260 400 0 25])
%         %
%     end
%     
%     legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10)
%     
%     figure
%     for i=1:20
%         Modes_pt=modes(i);
%         [index1,~]=find(flutter_data.POINT==Modes_pt);
%         
%         velocity=flutter_data.VELOCITY(index1);
%         
%         damping=flutter_data.DAMPING(index1);
%         
%         plot(velocity,damping,'-','LineWidth',1)
%         set(gcf,'Color','w')
%         xlabel('Velocity','Interpreter','latex','FontSize',12)
%         ylabel('Damping','Interpreter','latex','FontSize',12)
%         hold on
%         
%     end
%     
%     hold on
%     plot([0,300],[0,0],'k--')
%     legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10)
% 



%% Flutter results got by the FFAST model script 


%     run_folder1='C:\Git\A321_sizing\hg_codes\results\wing_flutter';
% 
%     flutter_data = h5read(strcat(run_folder1,'\flutter_test.h5'),'/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');
% 
%     modes=[1:35];
% 
%     figure
%     
%     for i=1:20
%         Modes_pt=modes(i);
%         [index1,~]=find(flutter_data.POINT==Modes_pt);
%         
%         velocity=flutter_data.VELOCITY(index1);
%         frequency=flutter_data.FREQUENCY(index1);
%         
%         
%         plot(velocity,frequency,'-','LineWidth',1)
%         
%         set(gcf,'Color','w')
%         xlabel('Velocity','Interpreter','latex','FontSize',12)
%         ylabel('Frequency','Interpreter','latex','FontSize',12)
%         hold on
%         %     axis([260 400 0 25])
%         %
%     end
%     
%     legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10)
%     
%     figure
%     for i=1:20
%         Modes_pt=modes(i);
%         [index1,~]=find(flutter_data.POINT==Modes_pt);
%         
%         velocity=flutter_data.VELOCITY(index1);
%         
%         damping=flutter_data.DAMPING(index1);
%         
%         plot(velocity,damping,'-','LineWidth',1)
%         set(gcf,'Color','w')
%         xlabel('Velocity','Interpreter','latex','FontSize',12)
%         ylabel('Damping','Interpreter','latex','FontSize',12)
%         hold on
%         
%     end
%     
%     hold on
%     plot([0,300],[0,0],'k--')
%     legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10)



%% Flutter results Fintan example


%     run_folder2='C:\Git\A321_sizing\hg_codes\results\wing_flutter\FWT_example';
% 
%     flutter_data = h5read(strcat(run_folder2,'\flutter_analysis.h5'),'/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');
% 
%     modes=[1:35];
% 
%     figure
%     
%     for i=1:4
%         Modes_pt=modes(i);
%         [index1,~]=find(flutter_data.POINT==Modes_pt);
%         
%         velocity=flutter_data.VELOCITY(index1);
%         frequency=flutter_data.FREQUENCY(index1);
%         
%         
%         plot(velocity,frequency,'-','LineWidth',1)
%         
%         set(gcf,'Color','w')
%         xlabel('Velocity','Interpreter','latex','FontSize',12)
%         ylabel('Frequency','Interpreter','latex','FontSize',12)
%         hold on
%         %     axis([260 400 0 25])
%         %
%     end
%     
%     legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10)
%     
%     figure
%     for i=1:10
%         Modes_pt=modes(i);
%         [index1,~]=find(flutter_data.POINT==Modes_pt);
%         
%         velocity=flutter_data.VELOCITY(index1);
%         
%         damping=flutter_data.DAMPING(index1);
%         
%         plot(velocity,damping,'-','LineWidth',1)
%         set(gcf,'Color','w')
%         xlabel('Velocity','Interpreter','latex','FontSize',12)
%         ylabel('Damping','Interpreter','latex','FontSize',12)
%         hold on
%         
%     end
%     
%     hold on
%     plot([0,300],[0,0],'k--')
%     legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10)
    
    
    
    %% FFAST MODEL
    
%     run_folder3='D:\pdf_nastran\pdf_nastran\FFAST_Test_FWT\FFAST Test\Source\Wing_only';
%     
%     flutter_data = h5read(strcat(run_folder3,'\flutter_analysis.h5'),'/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');
%     
%     modes=[1:35];
% 
%     figure
%     
%     for i=1:25
%         Modes_pt=modes(i);
%         [index1,~]=find(flutter_data.POINT==Modes_pt);
%         
%         velocity=flutter_data.VELOCITY(index1);
%         frequency=flutter_data.FREQUENCY(index1);
%         
%         
%         plot(velocity,frequency,'-','LineWidth',1)
%         
%         set(gcf,'Color','w')
%         xlabel('Velocity','Interpreter','latex','FontSize',12)
%         ylabel('Frequency','Interpreter','latex','FontSize',12)
%         hold on
%         %     axis([260 400 0 25])
%         %
%     end
%     
%     legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10)
%     
%     figure
%     for i=1:25
%         Modes_pt=modes(i);
%         [index1,~]=find(flutter_data.POINT==Modes_pt);
%         
%         velocity=flutter_data.VELOCITY(index1);
%         
%         damping=flutter_data.DAMPING(index1);
%         
%         plot(velocity,damping,'-','LineWidth',1)
%         set(gcf,'Color','w')
%         xlabel('Velocity','Interpreter','latex','FontSize',12)
%         ylabel('Damping','Interpreter','latex','FontSize',12)
%         hold on
%         
%     end
%     
%     hold on
%     plot([0,300],[0,0],'k--')
%     legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','13th Mode','Interpreter','latex','FontSize',10)


    
    %% flutter speed contour
    
%     x = linspace(-2,2);
% y = linspace(0,4);
% [X,Y] = meshgrid(x,y);
    
AR=[10,13,16,19,22];

ETA=[0,10,20,30,40];

[AR_mesh,ETA_mesh] = meshgrid(AR,ETA);

Flutter_Speed_matrix=[428,343,320,270,261;...
    425,351,312,292,296;...
    500,397,345,332,309;...
    500,500,380,358,369;500,500,490,397,367];

Flutter_Speed_matrix140=[393,350,305,274,255;...
    431,365,320,285,258;...
    460,396,362,332,307;...
    590,486,431,372,354;...
    600,528,469,408,378];

Flutter_Speed_matrix160=[406,355,291,284,265;...
    411,370,315,287,263;...
    469,381,321,312,309;...
    600,439,396,366,353;...
    750,512,412,422,345];


figure 

contourf(AR_mesh,ETA_mesh,Flutter_Speed_matrix,'ShowText','on')

colormap(jet(602))


xlabel('Aspect ratio','interpreter','latex')

ylabel('Folding wingtip size, $\eta$ ($\%$)','interpreter','latex')
set(gcf,'Color','w')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)



figure 

Flutter_Speed_matrix140(Flutter_Speed_matrix140>=500)=500;
contourf(AR_mesh,ETA_mesh,Flutter_Speed_matrix140,[200,280,300,320,340,360,380,400,420,440,460,480,500],'ShowText','on')

colormap(jet(602))


xlabel('Aspect ratio','interpreter','latex')

ylabel('Folding wingtip size, $\eta$ ($\%$)','interpreter','latex')
set(gcf,'Color','w')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)



figure 

Flutter_Speed_matrix160(Flutter_Speed_matrix160>=500)=500;
contourf(AR_mesh,ETA_mesh,Flutter_Speed_matrix160,'ShowText','on')

colormap(jet(602))


xlabel('Aspect ratio','interpreter','latex')

ylabel('Folding wingtip size, $\eta$ ($\%$)','interpreter','latex')


c=colorbar;
c.Location='southoutside';
c.Label.String = 'Flutter speed (m/s)';
c.Label.FontSize = 20;
c.FontSize=16;
c.TickLabelInterpreter='latex';
set(gca,'FontSize',14)   
set(gcf,'Color','w')
set(gca,'TickLabelInterpreter','latex')    