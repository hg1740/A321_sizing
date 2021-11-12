
Param=eval('A321_v1');

%% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[Param.Layout.Wing_Position,2,0]; %15
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
        
        FWT_R = insertWingFold(Wingbox_right, 'FlareAngle', Param.FWT.Flare_angle, 'FoldAngle', Param.FWT.Fold_angle,'EtaFold',Param.FWT.Fold_eta);
        FWT_R.HingeStiffness = [1e14 1e14 1e14 1e14 Param.FWT.Hinge_Stiffness 1e14];
        
        FWT_R.AeroPanelLength=[];
        FWT_R.NumAeroPanel=Param.Wing.AeroPanel_Number;
        FWT_R.AeroPanelAR=Param.Wing.AeroPanel_AR;
        
        FWT_R.NumBeamElem = 10;
        
        FWT_eta_= 0:0.1:1;
        FWT_Bheight=interp1([0,1],[Param.FWT.Root_Height,Param.FWT.Tip_Height],FWT_eta_);
        FWT_Bwidth=interp1([0,1],0.5*[Param.FWT.Root_Chord,Param.FWT.Tip_Chord],FWT_eta_);
        
        %%sizing variables ---------------------------------------------
            
        FWT_Spar_cap=Param.FWT.Thickness(1:11);
        FWT_Spar_web=Param.FWT.Thickness(12:22);
        FWT_Skin=Param.FWT.Thickness(23:33);
        FWT_Strg=Param.FWT.Thickness(34:44);
        
        fwt_d_strg=sqrt(FWT_Strg/0.36);
        fwt_t_strg=0.12*fwt_d_strg;
        
        % computing beam properties ------------------------------------
        
        [FWT_A_val,FWT_Iyy_val, FWT_Izz_val, FWT_J_val]=Beam_proeprties(FWT_Bwidth, FWT_Bheight, FWT_Spar_cap, FWT_Spar_web, FWT_Skin, FWT_Strg);
        
        % --------------------------------------------------------------
        
        % FWT structural properties
        
        FWT_R.A   =  FWT_A_val.cross_section;
        FWT_R.A_eta=FWT_eta_;
        
        FWT_R.I11 = FWT_Izz_val;
        FWT_R.I11_eta=FWT_eta_;
        
        FWT_R.I22 = FWT_Iyy_val;
        FWT_R.I22_eta = FWT_eta_;
        
        FWT_R.J   = FWT_J_val;
        FWT_R.J_eta= FWT_eta_;
        
        
        % FWT Jig Twist
        FWT_R.Twist = Param.FWT.Jig_Twist;
        FWT_R.Twist_eta = Param.FWT.Jig_Eta;
        
        
        
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
        for i=1:1:2
            handle=strcat('PM_right','i');
            handle=awi.model.PointMass;
            handle.SOffset=0.2+i*0.2;
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
    
    Thcikness_spar_cap=Param.Wing.Thickness(1:25);
    Thcikness_spar_web=Param.Wing.Thickness(26:50);
    Thcikness_skin=Param.Wing.Thickness(51:75);
    Area_strg=Param.Wing.Thickness(76:100);
    
    d_strg=sqrt(Area_strg/0.36);
    t_strg=0.12*d_strg;

    % -------------------------------------------------------------

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
    
    % computing beam properties -----------------------------------------
    
    [A_val,Iyy_val, Izz_val, J_val]=Beam_proeprties(Bwidth, Bheight, Thcikness_spar_cap, Thcikness_spar_web, Thcikness_skin, Area_strg);
    
    % -------------------------------------------------------------------

    eta_=etaRS;
    Y_eta=etaRS*RData(end);
    
    Wingbox_right.A   =  A_val.cross_section;
    Wingbox_right.A_eta=eta_;

    Wingbox_right.I11 = Izz_val;
    Wingbox_right.I11_eta=eta_;

    Wingbox_right.I22 = Iyy_val;
    Wingbox_right.I22_eta = eta_;

    Wingbox_right.J   = J_val;
    Wingbox_right.J_eta= eta_;


    % Jig shape
    
    Wingbox_right.Twist = Param.Wing.Jig_Twist;
    Wingbox_right.Twist_eta = Param.Wing.Jig_Eta;
    
    build(Wingbox_right)
    
    if isfield(Param,'FWT')
        
        Box_dimensions.Inboard.Height = [Bheight,FWT_Bheight(2:end)];
        Box_dimensions.Inboard.Width = [Bwidth,FWT_Bwidth(2:end)];
        Box_dimensions.Inboard.Stringer_length = [d_strg, fwt_d_strg(2:end)];
        Box_dimensions.Inboard.Stringer_thickness = [t_strg,fwt_t_strg(2:end)];
        
        Box_CrossSec.Izz=[Izz_val,FWT_Izz_val(2:end)];
        Box_CrossSec.Ixx=[Iyy_val,FWT_Iyy_val(2:end)];
        Box_CrossSec.A=[A_val.cross_section,FWT_A_val.cross_section(2:end)];
        
    else
        
        Box_dimensions.Inboard.Height = Bheight;
        Box_dimensions.Inboard.Width = Bwidth;
        Box_dimensions.Inboard.Stringer_length = d_strg;
        Box_dimensions.Inboard.Stringer_thickness = t_strg;
        
        Box_CrossSec.Izz=Izz_val;
        Box_CrossSec.Ixx=Iyy_val;
        Box_CrossSec.A=A_val.cross_section;
        
    end
    
    %% Mass definition
    

    wingmass_eta=etaRS;
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
    
    
    FEM_wing=convertToFE(Wingbox_right)
    
%     Width=[3,2];
%     Height=[0.7,0.5];
%     SparCap_Thickness=[0.04,0.01];
%     SparWeb_Thickness=[0.002,0.002];
%     Skin_Thickness=[0.003,0.002];
%     Strg_Area=[1e-4,1e-4];
%     
%     [Area,Iyy, Izz, J]=Beam_proeprties(Width,Height,SparCap_Thickness,SparWeb_Thickness, Skin_Thickness, Strg_Area)
    
    
    