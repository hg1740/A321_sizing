%% SOL 144
function [sigmaf1, sigmaf2]=CalcStress(t1,t2)

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
    E  = 76e9; %[N/m^2], typical YM of aluminium
    nu = 0.333;
    Mat = awi.model.Material;
    Mat.E  = E;
    Mat.Nu = nu;
    Mat.G  = E / (2 * (1 + nu));
    Wingbox_right.Material_eta = [0, 1];
    Wingbox_right.Material     = [Mat, Mat];

    % Define box beam corss section 1
    mybox_right1=awi.model.BoxBeam;
    mybox_right1.BoxType='SymmetricBox';
    mybox_right1.Height=0.78;
    mybox_right1.Width=3;
    
    %% sizing parameters-----------------------------------------------
    mybox_right1.CoverThickness= t2;
    mybox_right1.SparThickness = t1;
    % -----------------------------------------------------------------
    
    getGeometricProps(mybox_right1)
    
    
    % Define box beam corss section 2
    mybox_right2=awi.model.BoxBeam;
    mybox_right2.BoxType='SymmetricBox';
    mybox_right2.Height=0.6;
    mybox_right2.Width=1.8;
    
    %% sizing parameters-----------------------------------------------
    mybox_right2.CoverThickness= t2;
    mybox_right2.SparThickness = t1;
    % -----------------------------------------------------------------
    
    getGeometricProps(mybox_right2)
    
    mybox_right3=awi.model.BoxBeam;
    mybox_right3.BoxType='SymmetricBox';
    mybox_right3.Height=0.18;
    mybox_right3.Width=1;
    
    %% sizing parameters-----------------------------------------------
    mybox_right3.CoverThickness= t2;
    mybox_right3.SparThickness = t1;
    % -----------------------------------------------------------------
    
    getGeometricProps(mybox_right3)
    
    
%     Wingbox_right.BoxBeam = mybox_right;
    Wingbox_right.A   = [mybox_right1.Abb,mybox_right2.Abb,mybox_right3.Abb];
    Wingbox_right.A_eta=[0,0.3,1];
    
    Wingbox_right.I11 = [mybox_right1.Ixx,mybox_right2.Ixx,mybox_right3.Ixx];
    Wingbox_right.I11_eta=[0,0.3,1];
    
    Wingbox_right.I22 = [mybox_right1.Izz, mybox_right2.Izz, mybox_right3.Izz];
    Wingbox_right.I22_eta=[0,0.3,1];
    
    Wingbox_right.J   = [mybox_right1.Jbb,mybox_right2.Jbb,mybox_right3.Jbb];
    Wingbox_right.J_eta=[0,0.3,1];
    
    Wingbox_right.NSM = [mybox_right1.NSM,mybox_right2.NSM,mybox_right3.NSM];
    Wingbox_right.NSM_eta=[0,0.3,1];
    Wingbox_right.NSI = [mybox_right1.NSI,mybox_right2.NSI,mybox_right3.NSI];
    Wingbox_right.NSI_eta=[0,0.3,1];
    

    % Aeropanel definition
    %
    % AeroPanelLength
    %     NumAeroPanel
    Wingbox_right.NumAeroPanel=20;
    % Wingbox_right.AeroPanelLength=0.2;

    % %Control surfaces - flaps
    % myflap_right=awi.model.Flap;
    %
    % myflap_right.Eta=[0.0, 0.3];
    % myflap_right.xLE=[0.9,0.9];
    % myflap_right.xTE=[1,1];
    % myflap_right.Max_def=0.1;
    % myflap_right.Max_rate=0.1;
    % myflap_right.HingeLine='LE';
    % myflap_right.Label='flapR';
    % build(myflap_right)

    %Control surfaces - Aileron
    % myAileron_right=awi.model.Aileron;
    %
    % myAileron_right.Eta=[0.75, 0.92];
    % myAileron_right.xLE=[12/16,12/16];
    % myAileron_right.xTE=[1,1];
    % myAileron_right.Max_def=0.1;
    % myAileron_right.Max_rate=0.1;
    % myAileron_right.HingeLine='LE';
    % myAileron_right.Label='aileronR';
    % myAileron_right.NumAeroPanel=4;
    % build(myAileron_right)
    %
    % Wingbox_right.add([myAileron_right]);
    % Wingbox_right.ModelControlSurf = 1;

    % add point masses
    for i=1:1:5
        handle=strcat('PM_right','i');
        handle=awi.model.PointMass;
        handle.SOffset=0+i*0.2;
        handle.Mass=600;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle.MassGroup='Group1';
        Wingbox_right.add(handle);

    end
    
    
    engine_mass=awi.model.PointMass;
    
    engine_mass.SOffset=0.3;
    engine_mass.Mass=3440;
    Wingbox_right.add(engine_mass);


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
    E  = 76e19; %[N/m^2], typical YM of aluminium
    nu = 0.333;
    Mat1 = awi.model.Material;
    Mat1.E  = E;
    Mat1.Nu = nu;
    Mat1.G  = E / (2 * (1 + nu));

    
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
    Body.NSM = Bodybox.NSM;
    Body.NSI = Bodybox.NSI;


    % add point masses
    for i=1:1:10
        handle=strcat('PM_body','i');
        handle=awi.model.PointMass;
        handle.SOffset=-0.1+i*0.1;
        handle.Mass=2912;
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
    Tailwing_right.Material     = [Mat, Mat];

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
    Tailwing_right.NSM = tailbox_right.NSM;
    Tailwing_right.NSI = tailbox_right.NSI;

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
    Verticalwing.Material     = [Mat, Mat];

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
    Verticalwing.NSM = Verticalbox.NSM;
    Verticalwing.NSI = Verticalbox.NSI;


    %Control surfaces - elevators
    % myrudder=awi.model.ControlSurface;
    % myrudder.Eta=[0.2, 1];
    % myrudder.xLE=[0.6,0.6];
    % myrudder.xTE=[1,1];
    % myrudder.Max_def=0.1;
    % myrudder.Max_rate=0.1;
    % myrudder.HingeLine='LE';
    % myrudder.Label='rudder1';
    % myrudder.FaceColor='r';
    % build(myrudder)
    % Verticalwing.add(myrudder);
    %
    % Verticalwing.ModelControlSurf = 1;

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
    TrimLoadcase.LoadFactor = 1;
    build(TrimLoadcase)

    %% Generate the FEM - SOL 103 Modal analysis

    run_folder = [
        'D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1']; %[-], folder for exporting the NASTRAN model


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
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.xdb');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.h5');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.log');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.f06');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.f04');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.h5.*');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.log.*');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.f06.*');
    delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144*.f04.*');
      
    NastranMethods1.runNastran(trimFile);
    
    
    %% find nodes Y-coordinates on the wing
    
%     Wing_nodes=FEM_half.Children.Children(1, 1).Nodes;
%     
%     Coords=zeros(3,numel(Wing_nodes));
%     
%     for i=1:numel(Wing_nodes)
%         Coords(:,i)=Wing_nodes(i).X;
%     end
%     %  Coords(:,1)=Wing_nodes(1).X
%     Y=Coords(2,1:24);
    
    %% extract results

    F061=NastranMethods1.extractNastranResults('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\A320_half_model_SOL144.f06','ReadF06',true,'ReadHDF5',false);

    Max_momentX=F061.f06data.Bendingmoment.UMPLN2(69);
    Max_shearY=F061.f06data.Bendingmoment.USPLN2(69);
    Root_torque=F061.f06data.Bendingmoment.UTORQUE1(69);
    
    H=mybox_right1.Height;
    W=mybox_right1.Width;
    Ixx=mybox_right1.Ixx;
    
    % web: shear stress
    sigma12=Max_shearY/(2*H*t1) + abs(Root_torque)/(2*H*W*t1);
    sigmaf1=2.5*sqrt(3)*sigma12;
    
    % top: bending stress
    sigma11=0.5*Max_momentX*H/Ixx;
    sigmaf2=2.5*sigma11;
   
   
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



























