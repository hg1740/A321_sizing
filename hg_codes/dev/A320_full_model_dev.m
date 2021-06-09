  
%% sizing parameters thickness1 = spar, thickness2 = skin

%     x=[0.0184549632155560,0.0194332413084522,0.0207559132418991,0.0209074925629807,0.0231661967456849,0.0254828164202534,0.0264270050000000,0.0378563183137105,0.0384112626298514,0.0373093915208931,0.0364906994983588,0.0389743420000000,0.0395598526378275,0.0375818600059361,0.0370256249000000,0.0364906994983588,0.0362261419269957,0.0349193296635013,0.0362261419269957,0.0351743436550000,0.0305996900000000,0.0252980660012066,0.0153689282363496,0.000576499130128597,0.000697563947455602,0.00690405505625000,0.00700777451556970,0.00711050306418304,0.00726742637500000,0.00748474006756110,0.00770855196712667,0.00805255000000000,0.00811426522855439,0.00793621128715938,0.00770855196712667,0.00753940072280141,0.00726742637500000,0.00705895191696771,0.00675497791097389,0.00651130062423758,0.00614088900997626,0.00583384455947745,0.00542053909375000,0.00500179247918198,0.00444730593828125,0.00384085512851562,0.00319858816010073,0.00223369340704222,0.00180616690062519,0.00191511288486283,5.93348577610401e-05,5.93348577610401e-05,6.20049263602869e-05,6.47951480464998e-05,6.82054189963156e-05,7.17951778908585e-05,7.89746956799444e-05,7.84021291362648e-05,7.44820226794515e-05,7.12746628511498e-05,6.52683435371441e-05,6.20049263602869e-05,5.93348577610401e-05,5.39407797827637e-05,4.90370725297852e-05,4.45791568452593e-05,3.85001809118148e-05,3.34929803495562e-05,2.76801490492200e-05,2.17323484270735e-05,1.56247225182875e-05,9.21663720956289e-06,3.02451924944460e-06,2.23953621083449e-06,2.49960268549141e-06];

% x=[0.0201448760556999,0.0203464985649330,0.0210229435034526,0.0222854859088772,0.0239061265538741,0.0265003277217347,0.0296120063126451,0.0394300593581781,0.0395249693722001,0.0394643929529135,0.0390221955517316,0.0400055327454177,0.0401719190684001,0.0401311748430096,0.0400535888056135,0.0390779667807618,0.0392030563485772,0.0389275506463870,0.0385535305378241,0.0330889827800858,0.0295630826409310,0.0279177146013568,0.0231779673976581,0.0179446172737886,0.0116423894419668,0.00727239291671781,0.00738019012612562,0.00750789604354385,0.00765964663576572,0.00785097991358997,0.00806405639546637,0.00834195207736516,0.00849015733642170,0.00829297848323423,0.00806541884830507,0.00784522946941331,0.00759890220763609,0.00732641218925976,0.00703584338851878,0.00672735302799279,0.00639093094110575,0.00600154963916326,0.00558510407596679,0.00514700727696986,0.00482123945396051,0.00438326029543238,0.00381990809439375,0.00305816332295933,0.00227464932524216,0.000893351127432521,6.55062855293365e-05,6.73415238257884e-05,6.95304111960910e-05,7.21599182724826e-05,7.55660221360774e-05,7.97106432366773e-05,8.59133001855090e-05,8.88626004171758e-05,8.35398331898288e-05,7.81703733625587e-05,7.36519159672993e-05,6.90330954154385e-05,6.42077062678873e-05,5.91443879079160e-05,5.40136684723793e-05,4.84956809189784e-05,4.27774738720632e-05,3.69561684117309e-05,3.11605487896153e-05,2.70542660775571e-05,2.23060880812024e-05,1.68343413504158e-05,1.10387742468360e-05,6.31015569155482e-06,9.42950851226249e-07];
%     

% x=[0.0189596521054333,0.0191837427745745,0.0198539169260718,0.0210684229892251,0.0226246281693717,0.0250916230994881,0.0280488941974237,0.0367041175372602,0.0367553345728376,0.0366612768523305,0.0362145073754782,0.0370868191031210,0.0372021678219872,0.0371206354501631,0.0369988383777017,0.0360401301612646,0.0361049634501639,0.0356335403745682,0.0348916331838040,0.0296388702926122,0.0286701277141633,0.0277597217308701,0.0237803010813388,0.0191787360923511,0.0121862598372589,0.00699956238223479,0.00709931837280801,0.00721779831240336,0.00735837047471829,0.00753590622003866,0.00773387497843092,0.00799253834763674,0.00813544045254821,0.00794177830530960,0.00771903722029878,0.00750366043193169,0.00726346792293347,0.00699820361649310,0.00671555090260363,0.00641531700369381,0.00608847909600772,0.00571105500399573,0.00537365397847822,0.00508789516567232,0.00481404400255434,0.00438627803345660,0.00385033301066460,0.00312019520069395,0.00235313284362065,0.000845505438493128,6.07225713152272e-05,6.23599038943347e-05,6.43148014859593e-05,6.66620770111903e-05,6.96900056705333e-05,7.31110690050333e-05,7.77506255979472e-05,8.00010789470476e-05,7.59042231708472e-05,7.15881361765266e-05,6.75093353966714e-05,6.31981619690458e-05,5.87027230664615e-05,5.39922987204662e-05,4.92202772665568e-05,4.40999407554106e-05,3.88105719549454e-05,3.42919819496080e-05,3.05707906441678e-05,2.71569808553286e-05,2.25104864554471e-05,1.72665450978065e-05,1.15467937795578e-05,6.73738284228341e-06,8.31725206912222e-07];

x=[0.0160369400000000,0.0161179750000000,0.0165944330000000,0.0175266290000000,0.0187239040000000,0.0207096620000000,0.0230583810000000,0.0320089960000000,0.0320392420000000,0.0319427880000000,0.0315322530000000,0.0322935890000000,0.0323900670000000,0.0323115470000000,0.0321955620000000,0.0313429260000000,0.0311508860000000,0.0304592680000000,0.0297595010000000,0.0272724890000000,0.0271122090000000,0.0264611570000000,0.0230559960000000,0.0187199650000000,0.0107524740000000,0.00647614100000000,0.00657527000000000,0.00669537400000000,0.00683776700000000,0.00701687400000000,0.00721883300000000,0.00748254600000000,0.00755456300000000,0.00737352600000000,0.00716564300000000,0.00696475900000000,0.00674100500000000,0.00649426000000000,0.00623149000000000,0.00595214700000000,0.00564842500000000,0.00541126400000000,0.00522201100000000,0.00497677000000000,0.00469043300000000,0.00427877200000000,0.00376751800000000,0.00307441800000000,0.00231764800000000,0.000789440000000000,5.19000000000000e-05,5.34000000000000e-05,5.52000000000000e-05,5.74000000000000e-05,6.03000000000000e-05,6.35000000000000e-05,6.78000000000000e-05,6.88000000000000e-05,6.55000000000000e-05,6.18000000000000e-05,5.83000000000000e-05,5.45000000000000e-05,5.06000000000000e-05,4.66000000000000e-05,4.24000000000000e-05,3.80000000000000e-05,3.50000000000000e-05,3.25000000000000e-05,2.94000000000000e-05,2.60000000000000e-05,2.16000000000000e-05,1.67000000000000e-05,1.13000000000000e-05,6.56000000000000e-06,7.23000000000000e-07];
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
      
    tan_LE_sweep=(0.73*Semi_span*tan(LE_sweep*pi/180)-0.412*Root_chord)/(0.73*Semi_span);
    TE_sweep=atan(tan_LE_sweep)*180/pi; % deg
    
    Tip_chord=Root_chord-(Semi_span*tan(LE_sweep*pi/180)-Semi_span*0.73*tan(TE_sweep*pi/180));
    Taper_ratio=Tip_chord/Root_chord;
    
    Mean_cord_coefficient=(2/3)*(1+Taper_ratio+Taper_ratio^2)/(1+Taper_ratio);
    
    %% Mass configurations

    Payload_max=19900; % kg
%     Fuel_fraction=0.723; % percentage of fuel in the tank
    Fuel_fraction=0.56; % percentage of fuel in the tank
    Fuel_capacity=32940; % L

    MTOW=78000; % maximum take off mass
    OWE=42600;  % Operating empty mass
    MWE=44057;  % Manufacture's empty mass
    MZF=66000;  % Maxumum zero fuel mass
    

    Engine_mass=7362/2; % kg
    Pylon=1239/2; % kg
    Horizontal_tail=682; % kg
    Vertical_tail=522; % kg

    
    % estimate the wing weight: structure + secondary masss using
    % statistically 
%     Estimated_wing_weight=0.86*Wing_span*(Total_area*MZF*MTOW)^0.25; % give an initial guess of wing box weight
%     Estimated_wing_mass=Estimated_wing_weight/10; % give an initial guess of wing box weight. For A320, the wing mass = 8097 kg

    % For sized result, estimated wing mass = totol wing mass 2701.3 kg for
    % A 321 model
    
    Wing_total_mass=2510.9; %kg
    
    % The tube shell is assumed to be made by 4mm thick Aluminium skin. 
    
    Fuselage_shell_mass= 2*pi*2*0.004*2800*37.5; %kg
    
    Fuselage_structure_mass=(OWE - 2*Wing_total_mass - Horizontal_tail - Vertical_tail - Engine_mass*2 - Pylon*2 - Fuselage_shell_mass);
    
    % obtain lumped masses: Secondar_mass and fuel mass are the masses for
    % each wingbox. Fuselage_total_mass includs the structure mass and
    % payload
    
    % the sizing is carroied out using max. payload with its corresponding
    % max fuel mass 
    
    [Secondary_mass,Fuel_mass,Fuselage_total_mass]= Lumped_masses(MTOW,Wing_area,LE_sweep,Semi_span,Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);    
    
%     Secondary_mass=1;
%     Fuel_mass=1;
%     Fuselage_total_mass=1;
       
    %% Right root connector

    Connector_right = awi.model.LiftingSurface;
    Connector_right.Name = 'Connector_Right';
    Connector_right.Origin=[15,0,0]; %15
    

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
    Connector_box_right.CoverThickness=0.08;
    Connector_box_right.SparThickness=0.08;
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
    Connector_right.AeroPanelLength=0.5;

    build(Connector_right);
    
%     draw(Connector_right)
    
    
    %% Left root connector
    
    Connector_left = awi.model.LiftingSurface;
    Connector_left.Name = 'Connector_Right';
    Connector_left.Origin=[15,0,0]; %15
    
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
    Connector_box_left.CoverThickness=0.08;
    Connector_box_left.SparThickness=0.08;
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
    Connector_left.AeroPanelLength=0.5;

    build(Connector_left);
    
    draw(Connector_left)
   
    
  %% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[15,2,0]; %15
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';
    
    % Num of element 
    Wingbox_right.NumBeamElem = 23;

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = Semi_span;   %34.1/2;
    Wingbox_right.LESweep     = [LE_sweep, LE_sweep];
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [0, TE_sweep, TE_sweep];
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
    
    
    %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=Wingbox_right.NumBeamElem+2;
    
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
    Num_seg1=ceil(elnum*0.27); % number of elements in the inboard section
    Num_seg2=elnum - Num_seg1; % number of elements in the outboard section
    
    Num_sec1=Num_seg1+1;    % number of sections in the inboard section
    Num_sec2=Num_seg2+1;    % number of sections in the outboard section
    
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
    X_eng=(Wingbox_right.PanelCoords.LE.X(2)+Wingbox_right.PanelCoords.TE.X(2))/2;
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
    flap_R=awi.model.ControlSurface;
    flap_R.Eta=[0, 0.24];
    flap_R.xLE=[0.8,0.8];
    flap_R.xTE=[1,1];
    flap_R.Max_def=0.1;
    flap_R.Max_rate=0.1;
    flap_R.HingeLine='LE';
    flap_R.Label='FlapR';
    flap_R.FaceColor='r';
    
%     flap_R.NumAeroPanel=10;
    flap_R.AeroPanelLength=0.4;
    
    build(flap_R)
    Wingbox_right.add(flap_R);
    
    Wingbox_right.ModelControlSurf = 1;
    
    
    build(Wingbox_right);
    
%     FEM_test=convertToFE(Wingbox_right);
%     export(FEM_test, run_folder);
%     
%     draw(FEM_test)


%% Wingbox 2 - left and control surf.

    Wingbox_left = awi.model.LiftingSurface;
    Wingbox_left.Name = 'A320Wing_left';
    Wingbox_left.Origin=[15,-2,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_left.ActiveSet = 'sSet';

    %Wing dimensions
    Wingbox_left.SpanVector  = 'Y';
    Wingbox_left.Span        = -Semi_span;  
    Wingbox_left.LESweep     = [-LE_sweep, -LE_sweep];
    Wingbox_left.LESweep_eta = [0, 1];
    Wingbox_left.TESweep     = [0, -TE_sweep, -TE_sweep];
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
    
       
    % Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
%     Wingbox_left.NumAeroPanel=20;
    Wingbox_left.AeroPanelLength=0.4;
    
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
    X_eng=(Wingbox_left.PanelCoords.LE.X(2)+Wingbox_left.PanelCoords.TE.X(2))/2;
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
    flap_L=awi.model.ControlSurface;
    flap_L.Eta=[0, 0.24];
    flap_L.xLE=[0.8,0.8];
    flap_L.xTE=[1,1];
    flap_L.Max_def=0.1;
    flap_L.Max_rate=0.1;
    flap_L.HingeLine='LE';
    flap_L.Label='FlapL';
    flap_L.FaceColor='r';
    
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
    Body.Length=38;
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
    Tailwing_right.XOffset = 35;
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
    Verticalwing.XOffset=34;


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
%     Aircraft.RefChord = Aircraft.RefArea/Aircraft.RefSpan; 


    %% Generate the loadcase object
    TrimLoadcase1 = awi.model.LoadCase;
    


    acMass = 94000;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = mach_number*340;
    flap_angle=0;
    
%     FlightPoint=awi.model.FlightPoint;
%     FlightPoint.Mach=0.78;
%     FlightPoint.AcVelocity=FlightPoint.Mach*340;
%     FlightPoint.Altitude = 36000;
%     getFlightPointData(FlightPoint,'ISA');
    
    

    TrimLoadcase1.Name = 'A320_cruise_g';
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
    TrimLoadcase1.LoadFactor = 1;
    
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

    TrimLoadcase.Name = 'A320_cruise_1g';
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
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\full_A320_sizing_test2\range_calc']; %[-], folder for exporting the NASTRAN model

%     run_folder = [
%         'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\NSM_test']; %[-], folder for exporting the NASTRAN model


    % Convert to a finite element model and draw it
    FEM_full = convertToFE(Aircraft);

    % %Export it to a file
    export(FEM_full, run_folder);


  %% NASTRAN method - RUN SOL 144
  
%     % NastranExe: 'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe'
% 
%     NastranMethods1 = awi.methods.Nastran;
%     NastranMethods1.AnalysisModel = FEM_half;
%     MassCases=awi.model.MassCases.empty;
%     ID0=200;
%     % trimdata=NastranMethods1.getTrimData(FEM_half, Aircraft, TrimLoadcase, MassCases, ID0, NastranMethods1.RefNode);
    
    
     %% Run the analysis- SOL 144 static trim analysis 
    
%      NastranExe: 'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe'
     
     NastranMethods1 = awi.methods.Nastran;
     NastranMethods1.AnalysisModel = FEM_full;
     MassCases=awi.model.MassCases.empty;
     ID0=200;    
     
    trimFile1 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase1, MassCases,run_folder,'DatFilename','A320_cruise_g');
    
%     trimFile2 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase2, MassCases,run_folder,'DatFilename','A321_cruise_1g');
    
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\A320_half_model_SOL144*.*')
   
    delete(strcat(run_folder, '\A320_cruise_g*.xdb.*'));
    delete(strcat(run_folder, '\A320_cruise_g*.h5.*'));
    delete(strcat(run_folder, '\A320_cruise_g*.log.*'));
    delete(strcat(run_folder, '\A320_cruise_g*.f06.*'));
    delete(strcat(run_folder, '\A320_cruise_g*.f04.*'));
    delete(strcat(run_folder, '\A320_cruise_g*.xdb'));
    delete(strcat(run_folder, '\A320_cruise_g*.h5'));
    delete(strcat(run_folder, '\A320_cruise_g*.log'));
    delete(strcat(run_folder, '\A320_cruise_g*.f06'));
    delete(strcat(run_folder, '\A320_cruise_g*.f04'));
    delete(strcat(run_folder, '\ajj.op4'));
    delete(strcat(run_folder, '\ffaj.op4'));
    
    delete(strcat(run_folder, '\A320_cruise_1g*.xdb.*'));
    delete(strcat(run_folder, '\A320_cruise_1g*.h5.*'));
    delete(strcat(run_folder, '\A320_cruise_1g*.log.*'));
    delete(strcat(run_folder, '\A320_cruise_1g*.f06.*'));
    delete(strcat(run_folder, '\A320_cruise_1g*.f04.*'));
      
    NastranMethods1.runNastran(trimFile1);
%     NastranMethods1.runNastran(trimFile2);
    
    
 %% lift and drag calculateion   
   
Ajj=mni.result.op4(strcat(run_folder,'/ajj.op4')).read_matrix();
FFaj=mni.result.op4(strcat(run_folder,'/ffaj.op4')).read_matrix();

res_aeroF = mni.result.f06(strcat(run_folder,'/a320_cruise_g.f06')).read_aeroF;
% 
FlightPoint=awi.model.FlightPoint;
FlightPoint.Mach=TrimLoadcase1.Mach;
FlightPoint.AcVelocity=FlightPoint.Mach*340;
FlightPoint.Altitude = TrimLoadcase1.Altitude;
getFlightPointData(FlightPoint,'ISA');

q = FlightPoint.DynPressure;

WJ = Ajj*(FFaj./q);
% surface_normal = model.fwt_normal_vector();
% drag_mag = dot([1 0 0],surface_normal);

idx = 415:906; % wing
idx_c=1:52; %conn
idx_t=105:195; %tail wing 

idx_lw=907:1398;% left wing
idx_lc=53:104;
idx_lt=196:286;


lift_wing=res_aeroF.aeroFz(idx)';
lift_conn=res_aeroF.aeroFz(idx_c)';
lift_tail=res_aeroF.aeroFz(idx_t)';

lift_left_wing=res_aeroF.aeroFz(idx_lw)';
lift_left_conn=res_aeroF.aeroFz(idx_lc)';
lift_left_tail=res_aeroF.aeroFz(idx_lt)';


all_lift1=sum(lift_wing) + sum(lift_conn) + sum(lift_tail);

all_lift2=sum(lift_left_wing) + sum(lift_left_conn) + sum(lift_left_tail);

% CL CD

Drag_wing = sin(WJ(idx)-0.075).*res_aeroF.aeroFz(idx)';
Drag_conn = sin(WJ(idx_c)-0.075).*res_aeroF.aeroFz(idx_c)';
Drag_tail = sin(WJ(idx_t)-0.075).*res_aeroF.aeroFz(idx_t)';

% t=atan(WJ/5);
% Drag_wing = sin(t(idx)).*res_aeroF.aeroFz(idx)';
% Drag_conn = sin(t(idx_c)).*res_aeroF.aeroFz(idx_c)';
% Drag_tail = sin(t(idx_t)).*res_aeroF.aeroFz(idx_t)';

% figure 
% plot(WJ(idx),'r.')


Total_drag=sum(Drag_wing)+sum(Drag_conn)+sum(Drag_tail);

Ref_surf=Wingbox_right.SurfaceArea+Connector_right.SurfaceArea+Tailwing_right.SurfaceArea;

Cd=Total_drag/(q*Ref_surf);

Cl=all_lift1/(q*Ref_surf); 

%% test
% a=inv(Ajj);
% vor=Ajj\WJ;
% 
% alpha=Ajj\(res_aeroF.aeroFz'./q);
% 
% % alpha=atan(WJ/230);
% 
% Drag1 = sin(0.09).*res_aeroF.aeroFz(idx)';
% Drag2 = sin(0.09).*res_aeroF.aeroFz(idx_c)';
% Drag3 = sin(0.09).*res_aeroF.aeroFz(idx_t)';
% 
% Total=sum(Drag1)+sum(Drag2)+sum(Drag3);
% Cdd=Total/(q*Ref_surf);


%     result plotting for SOL 144

    
    Trim1_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A320_cruise_g.f06'),'ReadF06',true,'ReadHDF5',false);
%     Trim2_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_cruise_1g.f06'),'ReadF06',true,'ReadHDF5',false);
%     
%     WingNodes=NastranMethods1.WingNode;
%     
%     index1=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[WingNodes(1:end).GID]);
%     index2=ismember([Trim1_res.f06data.Bendingmoment.LGRID],[WingNodes(end).GID]);
%     
%    
%     % extract forces:
%     
%     % forces from trim1
%     Trim1.M_P1=[Trim1_res.f06data.Bendingmoment.UMPLN1(index1),Trim1_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
%     Trim1.M_P2=[Trim1_res.f06data.Bendingmoment.UMPLN2(index1),Trim1_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
%     Trim1.T=[Trim1_res.f06data.Bendingmoment.UTORQUE1(index1),Trim1_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
% 
%     Trim1.S_P1=[Trim1_res.f06data.Bendingmoment.USPLN1(index1),Trim1_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
%     Trim1.S_P2=[Trim1_res.f06data.Bendingmoment.USPLN2(index1),Trim1_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
%     
%     % forces from trim2
%     Trim2.M_P1=[Trim2_res.f06data.Bendingmoment.UMPLN1(index1),Trim2_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
%     Trim2.M_P2=[Trim2_res.f06data.Bendingmoment.UMPLN2(index1),Trim2_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
%     Trim2.T=[Trim2_res.f06data.Bendingmoment.UTORQUE1(index1),Trim2_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
% 
%     Trim2.S_P1=[Trim2_res.f06data.Bendingmoment.USPLN1(index1),Trim2_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
%     Trim2.S_P2=[Trim2_res.f06data.Bendingmoment.USPLN2(index1),Trim2_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
%        
%     % Maximum forces used for sizing
%     M_P1=max([abs(Trim1.M_P1); abs(Trim2.M_P1)]);
%     M_P2=max([abs(Trim1.M_P2); abs(Trim2.M_P2)]);
%     T=max([abs(Trim1.T); abs(Trim2.T)]);
%     
%     S_P1=max([abs(Trim1.S_P1); abs(Trim2.S_P1)]);
%     S_P2=max([abs(Trim1.S_P2); abs(Trim2.S_P2)]);  
%     
%     Grid_coord = h5read(strcat(run_folder,'\A321_cruise_2p5g.h5'),'/NASTRAN/INPUT/NODE/GRID');
%     Displacement = h5read(strcat(run_folder,'\A321_cruise_2p5g.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');
%     
%     index3=ismember([Grid_coord.ID],[WingNodes(1:end).GID]);
%     
%     Y=Grid_coord.X(2,index3);
%     Displacement_Z=Displacement.Z(index3);
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


    
%     NastranMethods1 = awi.methods.Nastran;
%     NastranMethods1.AnalysisModel = FEM_half;
%     MassCases=awi.model.MassCases.empty;
%     RefGrid=NastranMethods1.RefNode;
% 
%     fid = fopen(strcat(run_folder,'\NastranHeaderFile.dat'),'r');
%     
%     cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true );
%     cac=cac{1};
%     fclose( fid );
%     
%     mid_line0=22;
%     mid_line=length(cac)-1; 
%     fid = fopen( strcat(run_folder,'\full_model_103.dat'), 'w' );
%     
%     % write existing head file
%     for jj = 1 : mid_line0
%         fprintf(fid,'%s\n',cac{jj})       
%     end
%     
% %     fprintf(fid,'SUPORT = %i\r\n',201)
%     fprintf(fid,'SPC = %i\r\n',202)
%     
%     for ii =  mid_line0 : mid_line 
%         fprintf(fid,'%s\n',cac{ii})       
%     end
%     
%     
%     %write spc part
%     line='$.1.....2.......3.......4.......5.......6.......7.......8.......9.......10......\r\n';
%     fprintf(fid,line);
%     ID0=250;
%     dof = 35;
%     spc = 246;
%     SPC_id = ID0;
%     
%     TrimData.RefGrid = RefGrid;
%     TrimData.SPC_id  = SPC_id;
%     TrimData.SPCdof  = spc;
%     TrimData.SUPdof  = dof;
%                         
% %     spc_format='%-8s%-8i%-8i%-8i\r\n';
% %     fprintf(fid,spc_format,'SPC1',1,123456,1005);
%     
%     fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPCADD', TrimData.SPC_id+2, TrimData.SPC_id, TrimData.SPC_id+1);
%     %   - SPC
%     rows=floor((numel(TrimData.RefGrid)-3)/5);
%     remd=rem((numel(TrimData.RefGrid)-3),5);
%     fprintf(fid, '%-8s%-8i%-8i%-8i%-8i%-8i\r\n', 'SPC1', TrimData.SPC_id, TrimData.SPCdof, TrimData.RefGrid(1:3).GID);
%     fprintf(fid, '        %-8i%-8i%-8i%-8i%-8i\r\n', TrimData.RefGrid(4:4+rows*5-1).GID);
%     format_last=strcat(repmat('%-8i',1,remd+1),'\r\n');
%     fprintf(fid, format_last, '', TrimData.RefGrid(4+rows*5:end).GID);
%     
%     %charles added SPC 2
%     % -SPC 2 at COG
%     fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', TrimData.SPC_id+1, 1246, TrimData.RefGrid(34).GID);
%     fprintf(fid, '%-8s%-8i%-8i\r\n', 'SUPORT', TrimData.RefGrid(34).ID, TrimData.SUPdof);
%     
%     % write the end 
%     fprintf(fid,'%s\n',cac{end});


    %% Run the analysis - SOL 145 flutter 
    
%         % NastranExe: 'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe'
% 
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
    
    %% Run the analysis - SOL 146 gust 
        GustLoadcase = awi.model.LoadCase;
    %     GustLoadcase.LoadCaseType = 'Pratt Gust';
        GustLoadcase.Altitude   = 36000;
        GustLoadcase.AcVelocity = 0.78*340;
        GustLoadcase.AcMass = 500;
        GustLoadcase.Mach = 0.78;
        GustLoadcase.GustLength = linspace(18,214,10);%[18:20:214];
%         GustLoadcase.GustLength=20;
        GustLoadcase.GustDirection=-1;
    
        FlightPoint=awi.model.FlightPoint;
        FlightPoint.Mach=0.78;
        FlightPoint.AcVelocity=FlightPoint.Mach*340;
        FlightPoint.Altitude = 36000;
        getFlightPointData(FlightPoint,'ISA');
    
        NastranMethods1 = awi.methods.Nastran;
        NastranMethods1.AnalysisModel = FEM_full;
        MassCases=awi.model.MassCases.empty;
    
        gustfile=NastranMethods1.writeGustFile(Aircraft, GustLoadcase, MassCases, FlightPoint, run_folder,'DatFilename','gust_analysis_new');
    %     gustfile=NastranMethods1.writeGustFile(Aircraft, GustLoadcase, MassCases, FlightPoint, run_folder,'DatFilename','gust_analysis_add_negative');
        NastranMethods1.runNastran(gustfile);
%     
%     
% %     code test 
%     
%     L=GustLoadcase.GustLength; % gust length
% %     L=[200];
%     FlightPoints=FlightPoint;
%     N=numel(L);
%     Ng=64; %number of points
%     bulk_data=zeros(N,Ng*2);
%     Uds=zeros(1,N); % data record for gust amplitude
%     
%     for n=1:N
%         H=L(n)*0.5;
%         Kg=1; % gust alleviation factor
%         
%         
%         tmax=L(n)/FlightPoints.AcVelocity;
%         ts=linspace(0,tmax,Ng);
% %         t0=linspace(tmax,tmax+6,64);
%         
%         
%         %Calculate gust velocity
%         Uref = [ ...
%             0    , 15000, 60000 ; ...
%             17.07, 13.41, 6.36 ];
%         Ug = interp1(Uref(1, :), Uref(2, :),FlightPoints.Altitude);  %[m/s], EAS
%         U_ref = Ug ./ sqrt(FlightPoints.DensityRatio);  %[m/s], TAS
%         
%         Uds(n)=U_ref*Kg*(H/107)^(1/6);
%         
%         gust_1mc=GustLoadcase.GustDirection*Uds(n)*(1-cos(2*pi*FlightPoints.AcVelocity*ts/L(n)))*0.5;
%         
% %         gust_1mc_non=zeros(1,numel(t0));
% %         
% %         gust_1mc_opp=-Uds(n)*(1-cos(2*pi*FlightPoints.AcVelocity*(ts)/L(n)))*0.5;
% %         
% %         t_total=[ts,t0,ts+6];
% %         gust_total=[gust_1mc,gust_1mc_non,gust_1mc_opp];
%         
%         
%         
% %         plot(ts*FlightPoints.AcVelocity,gust_1mc,'b.')
%         plot(ts*FlightPoints.AcVelocity,gust_1mc,'b.')
%         xlabel('Gust length (m)','FontSize',12,'Interpreter','latex')
%         ylabel('Gust velocity (m/s)','FontSize',12,'Interpreter','latex')
%         set(gcf,'color','w')
%         hold on
%    
%         
%         
%     end
% % % %     
% % % % %     Fs=20;
% % % % %     n = 2^nextpow2(length(t_total));
% % % % %     Y = fft(gust_total,n);
% % % % %     f = Fs*(0:(n/2))/n;
% % % % %     P = abs(Y/n).^2;
% % % % %     
% % % % %     figure  
% % % % %     plot(f,P(1:n/2+1))
% % % % %     title('Gaussian Pulse in Frequency Domain')
% % % % %     xlabel('Frequency (f)')
% % % % %     ylabel('|P(f)|^2')
% % % 
% % % 
% % for n=1:N
% %     H=L(n)*0.5;
% %     Kg=1; % gust alleviation factor
% %     
% %     
% %     
% %     tmax=L(n)/FlightPoints.AcVelocity;
% %     delta=8-tmax;
% %     
% %     ts=linspace(0,tmax,Ng); % positive gust 1
% %     t0=linspace(tmax+0.05,tmax+0.05+delta,Ng); % space 1
% %     t1=linspace(tmax+0.1+delta,2*tmax+0.1+delta,Ng);% negative gust 1
% %     t2=linspace(2*tmax+0.15+delta,2*tmax+0.15+2*delta,Ng);% space 2
% %     t3=linspace(2*tmax+0.2+2*delta,3*tmax+0.2+2*delta,Ng);% positive gust 2
% %     t4=linspace(3*tmax+0.25+2*delta,3*tmax+0.25+3*delta,Ng);% space 3
% %     
% %     
% %     %Calculate gust velocity
% %     Uref = [ ...
% %         0    , 15000, 60000 ; ...
% %         17.07, 13.41, 6.36 ];
% %     Ug = interp1(Uref(1, :), Uref(2, :),FlightPoints.Altitude);  %[m/s], EAS
% %     U_ref = Ug ./ sqrt(FlightPoints.DensityRatio);  %[m/s], TAS
% %     
% %     Uds(n)=U_ref*Kg*(H/107)^(1/6);
% %     
% %     gust_1mc=Uds(n)*(1-cos(2*pi*FlightPoints.AcVelocity*ts/L(n)))*0.5;
% %     
% %     gust_1mc_non=zeros(1,numel(t0));
% %     
% %     gust_1mc_opp=-2*Uds(n)*(1-cos(2*pi*FlightPoints.AcVelocity*(ts)/L(n)))*0.5;
% %     
% %     t_total=[ts,t0,t1,t2,t3,t4];
% %     gust_total=[gust_1mc,gust_1mc_non,gust_1mc_opp,gust_1mc_non,gust_1mc,gust_1mc_non];
% %     
% %     plot(t_total,gust_total,'b.')
% %     xlabel('Gust length (m)','FontSize',12,'Interpreter','latex')
% %     ylabel('Gust velocity (m/s)','FontSize',12,'Interpreter','latex')
% %     set(gcf,'color','w')
% %     hold on
% %     
% % %     for i=0:6*numel(ts)-1
% % %         
% % %         index1=1+2*i;
% % %         index2=2+2*i;
% % %         bulk_data(n,index1)=t_total(i+1);
% % %         bulk_data(n,index2)=gust_total(i+1)/max(gust_total);
% %         
% % %     end
% %     
% % end
%     
% 
%     %% Result extraction
%     
%     % Bending moment & torque 
% %     Gust_res=Gust_extract(strcat(run_folder,'\gust_analysis.f06'));
% %     Steps=201;
% 
    Gust_res=Gust_extract(strcat(run_folder,'\gust_analysis_new.f06'));
    Steps=201;
     
    NumGust=numel(GustLoadcase.GustLength);

    t=linspace(0,9,Steps);
    
    M_root=zeros(NumGust,Steps);
    T_root=zeros(NumGust,Steps);

    peaks=zeros(NumGust,2);
    peaks_torque=zeros(NumGust,2);
    
    WingNodes=NastranMethods1.WingNode;
    % find the nodes at the root
    indexR=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[WingNodes(1).GID]);
    % find all nodes on the wing
    indexAll=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[WingNodes(1:end).GID]);


    
    for j=1:NumGust
        
        for i=1:Steps
%             f_name=fieldnames(F);
            
            M_root(j,i)=Gust_res.Bendingmoment(i+(j-1)*Steps).UMPLN2(indexR); % out of plane moment 69 at root
            T_root(j,i)=Gust_res.Bendingmoment(i+(j-1)*Steps).UTORQUE1(indexR); % out of plane moment 69 at root
            
            
        end
        
        peaks(j,:)=[max(M_root(j,:)),min(M_root(j,:))];
        peaks_torque(j,:)=[max(T_root(j,:)),min(T_root(j,:))];
        
    end
    
figure % Time history of the Root bending moment  
plot(t,M_root(1,:),'b.')
hold on 
plot(t,M_root(2,:),'r.')
hold on 
plot(t,M_root(3,:),'k.')
hold on 
plot(t,M_root(4,:),'m.')
hold on 
plot(t,M_root(5,:),'g.')
hold on 
plot(t,M_root(6,:),'b.')
hold on 
plot(t,M_root(7,:),'r.')
hold on 
plot(t,M_root(8,:),'k.')
hold on 
plot(t,M_root(9,:),'m.')
hold on 
plot(t,M_root(10,:),'g.')
% 
% xlabel('Time (s)','FontSize',12,'Interpreter','latex')
% ylabel('Root Bending moment (Nm)','FontSize',12,'Interpreter','latex')
% set(gcf,'color','w')
% legend('L = 20m','L = 40m','L = 60m','L = 80m','L = 100m','L = 120m','L = 140m','L = 160m','L = 180m','L = 200m','Interpreter','latex')
% 
% 
% figure % Time history of the Root torque
% 
% plot(t,T_root(1,:),'b.')
% hold on 
% plot(t,T_root(2,:),'r.')
% hold on 
% plot(t,T_root(3,:),'k.')
% hold on 
% plot(t,T_root(4,:),'m.')
% hold on 
% plot(t,T_root(5,:),'g.')
% hold on 
% plot(t,T_root(6,:),'b.')
% hold on 
% plot(t,T_root(7,:),'r.')
% hold on 
% plot(t,T_root(8,:),'k.')
% hold on 
% plot(t,T_root(9,:),'m.')
% hold on 
% plot(t,T_root(10,:),'g.')
% 
% xlabel('Time (s)','FontSize',12,'Interpreter','latex')
% ylabel('Root torque (Nm)','FontSize',12,'Interpreter','latex')
% set(gcf,'color','w')
% legend('L = 20m','L = 40m','L = 60m','L = 80m','L = 100m','L = 120m','L = 140m','L = 160m','L = 180m','L = 200m','Interpreter','latex')
% % legend('L = 20m','L = 60m','L = 100m','L = 140m','L = 180m','Interpreter','latex')
% 
% 
% 
% %  peak plot 
% figure 
% plot(GustLoadcase.GustLength,peaks(:,1),'b-s')
% hold on 
% plot(GustLoadcase.GustLength,peaks(:,2),'r-s')
% 
% xlabel('Gust length (m)','FontSize',12,'Interpreter','latex')
% ylabel('Root Bending moment (Nm)','FontSize',12,'Interpreter','latex')
% set(gcf,'color','w')
% % 
% % legend('L = 10 m','L = 90 m','L = 190 m')
% 
% %  peak plot 
% figure 
% plot(GustLoadcase.GustLength,peaks_torque(:,1),'b-s')
% hold on 
% plot(GustLoadcase.GustLength,peaks_torque(:,2),'r-s')
% 
% xlabel('Gust length (m)','FontSize',12,'Interpreter','latex')
% ylabel('Root torque (Nm)','FontSize',12,'Interpreter','latex')
% set(gcf,'color','w')
% 
% 
% % potato plot
% M_total=[M_root(1,:)';M_root(2,:)';M_root(3,:)';M_root(4,:)';M_root(5,:)';M_root(6,:)';M_root(7,:)';M_root(8,:)';M_root(9,:)';M_root(10,:)'];
% T_total=[-T_root(1,:)';-T_root(2,:)';-T_root(3,:)';-T_root(4,:)';-T_root(5,:)';-T_root(6,:)';-T_root(7,:)';-T_root(8,:)';-T_root(9,:)';-T_root(10,:)'];
% P=[M_total,T_total];
% 
% [k,av] = convhull(P);
% 
% figure % Time history of the Root bending moment  
% plot(M_root(1,:),-T_root(1,:),'b.')
% hold on
% plot(M_root(2,:),-T_root(2,:),'r.')
% hold on
% plot(M_root(3,:),-T_root(3,:),'m.')
% hold on
% plot(M_root(4,:),-T_root(4,:),'g.')
% hold on
% plot(M_root(5,:),-T_root(5,:),'k.')
% hold on
% plot(M_root(6,:),-T_root(6,:),'b.')
% hold on
% plot(M_root(7,:),-T_root(7,:),'r.')
% hold on
% plot(M_root(8,:),-T_root(8,:),'m.')
% hold on
% plot(M_root(9,:),-T_root(9,:),'g.')
% hold on
% plot(M_root(10,:),-T_root(10,:),'k.')
% hold on
% plot(P(k,1),P(k,2),'k-','LineWidth',1)
% 
% xlabel('Root bending moment (Nm)','FontSize',12,'Interpreter','latex')
% ylabel('Root torque (Nm)','FontSize',12,'Interpreter','latex')
% set(gcf,'color','w')
% 
% 
% %% code test - bending moment distribution 
% 
%     NumGust=numel(GustLoadcase.GustLength);
% 
%     t=linspace(0,5,Steps);
%     
%     M_wing=zeros(NumGust,Steps);
%     T_wing=zeros(NumGust,Steps);
% 
%     peaks_wing=zeros(NumGust,2);
%     peaks_torque_wing=zeros(NumGust,2);
%     
%     wing_pt=75:1:100;
%     
%     Moment_total=[];
%     Torque_total=[];
%     Shear_total=[];
%     
%     for j=1:NumGust
%         
%         for k=1:24
%             
%             for i=1:Steps
%                 
%                 % find the index of the nodes 
%                 temp_index=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[WingNodes(k).GID]);
%                                
%                 M_wing(k,i)=Gust_res.Bendingmoment(i+(j-1)*Steps).UMPLN2(temp_index); 
%                 T_wing(k,i)=Gust_res.Bendingmoment(i+(j-1)*Steps).UTORQUE1(temp_index);
%                 V_wing(k,i)=Gust_res.Bendingmoment(i+(j-1)*Steps).USPLN2(temp_index);
%                 
%             end
%             
%             
%                         
%         end
%         
%         M_matrix(j).moment=M_wing;
%         T_matrix(j).torque=T_wing;
%         V_matrix(j).shear=V_wing;
%         
%         Moment_total=horzcat(Moment_total,M_matrix(j).moment);
%         Torque_total=horzcat(Torque_total,T_matrix(j).torque);
%         Shear_total=horzcat(Shear_total,V_matrix(j).shear);
%                
%     end
%     
%     
%     max_M= max(Moment_total, [], 2);
%     min_M= min(Moment_total, [], 2);
%     
%     max_T= max(Torque_total, [], 2);
%     min_T= min(Torque_total, [], 2);
%     
%     max_V= max(Shear_total, [], 2);
%     min_V= min(Shear_total, [], 2);
%     
%     xnum=1:1:24;
%     
%     figure
%     
%     plot(xnum, max_M,'b.')
%     hold on
%     plot(xnum, min_M,'b.')
%     
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Bending moment (Nm)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
%     
%     figure
%     
%     plot(xnum, max_T,'b.')
%     hold on
%     plot(xxnum, min_T,'b.')
%     
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Torque (Nm)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
%     
%     figure
%     
%     plot(xnum, max_V,'b.')
%     hold on
%     plot(xnum, min_V,'b.')
%     
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Torque (Nm)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
%     
% % legend('L = 20m','L = 40m','L = 60m','L = 80m','L = 100m','L = 120m','L = 140m','L = 160m','L = 180m','L = 200m','Interpreter','latex')
%     
% 
% 
% 
% %%  z-motion of cog and wing tip
% 
% % Gust_data = h5read(strcat(run_folder,'\gust_analysis.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');
% % Gust_data = h5read(strcat(run_folder,'\gust_analysis_long.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');
% Gust_data = h5read(strcat(run_folder,'\gust_analysis_new.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');
% Z_disp_result=ones(Steps*10,numel(WingNodes));
% 
% 
% %group 1: wing Z displacement ------------------------------------------------------------------
% 
% WingNodes=NastranMethods1.WingNode;
% % ID=WingNodes(1:25).GID;
% 
% % Z_disp_wing=ones
% for i=1:numel(WingNodes)
%     
%     [index,~]=find(Gust_data.ID==WingNodes(i).GID); 
%     
%     Z_disp_result(:,i)=Gust_data.Z([index]);
%     
% end
% 
% gust1= Z_disp_result(1:Steps,:);
% gust2= Z_disp_result(Steps+1:Steps*2,:);
% gust3= Z_disp_result(Steps*2+1:Steps*3,:);
% gust4= Z_disp_result(Steps*3+1:Steps*4,:);
% gust5= Z_disp_result(Steps*4+1:Steps*5,:);
% gust6= Z_disp_result(Steps*5+1:Steps*6,:);
% gust7= Z_disp_result(Steps*6+1:Steps*7,:);
% gust8= Z_disp_result(Steps*7+1:Steps*8,:);
% gust9= Z_disp_result(Steps*8+1:Steps*9,:);
% gust10 = Z_disp_result(Steps*9+1:Steps*10,:);
% 
% % curve = animatedline;
% 
% % 
% % for i=1:201
% %     clearpoints(curve)
% %     x=1:1:25;
% %     addpoints(curve,x,gust1(i,:))
% %     
% %    drawnow
% %     
% %     
% % end
% close all
% x=[1:1:25]*16/25;
% 
% y=gust1;
% 
% figure 
% ph=plot(x,y(100,:),'b-s');
% set(gca,'Xlim',[0 25],'Ylim',[-0.5 0.5])
% set(gcf,'Color','w')
% xlabel('Wing span','Interpreter','latex')
% ylabel('Vertical displacement','Interpreter','latex')
% 
% for i=1:Steps
%     ph.XData=x;
%     ph.YData=y(i,:);
%     
%     drawnow
%     pause(0.05)
%     
% end
%     
% %group 1: fuselage Z displacement ------------------------------------------------------------------
% 
% FuselageNodes=NastranMethods1.RefNode;
% num=numel(FuselageNodes);
% Z_disp_fuselage=ones(Steps*10,numel(FuselageNodes));
% 
% 
% % Z_disp_wing=ones
% for i=1:numel(FuselageNodes)
%     
%     [index,~]=find(Gust_data.ID==FuselageNodes(i).GID); 
%     
%     Z_disp_fuselage(:,i)=Gust_data.Z([index]);
%     
% end
% 
% gust1_Z_fuselage = Z_disp_fuselage(1:Steps,:);
% gust2_Z_fuselage= Z_disp_fuselage(Steps+1:Steps*2,:);
% gust3_Z_fuselage= Z_disp_fuselage(Steps*2+1:Steps*3,:);
% gust4_Z_fuselage= Z_disp_fuselage(Steps*3+1:Steps*4,:);
% gust5_Z_fuselage= Z_disp_fuselage(Steps*4+1:Steps*5,:);
% gust6_Z_fuselage= Z_disp_fuselage(Steps*5+1:Steps*6,:);
% gust7_Z_fuselage= Z_disp_fuselage(Steps*6+1:Steps*7,:);
% gust8_Z_fuselage= Z_disp_fuselage(Steps*7+1:Steps*8,:);
% gust9_Z_fuselage= Z_disp_fuselage(Steps*8+1:Steps*9,:);
% gust10_Z_fuselage = Z_disp_fuselage(Steps*9+1:Steps*10,:);
% 
% x=1:1:num;
% y=gust6_Z_fuselage;
% 
% figure 
% ph=plot(x,y(100,:),'bs-');
% set(gca,'Xlim',[0 75],'Ylim',[-2 2])
% set(gcf,'Color','w')
% xlabel('Fuselage','Interpreter','latex')
% ylabel('Vertical displacement','Interpreter','latex')
% 
% for i=1:Steps
%     ph.XData=x;
%     ph.YData=y(i,:);
%     
%     drawnow
%     pause(0.05)
%     
% end
% 
% 
% 
% % group 3: wing tip twist  ------------------------------------------------------------------
% 
% WingNodes=NastranMethods1.WingNode;
% % ID=WingNodes(1:25).GID;
% Wing_twist=ones(Steps*10,numel(WingNodes));
% 
% % Z_disp_wing=ones
% for i=1:numel(WingNodes)
%     
%     [index,~]=find(Gust_data.ID==WingNodes(i).GID); 
%     
%     Wing_twist(:,i)=Gust_data.RX([index]);
%     
% end
% 
% gust1t= Wing_twist(1:Steps,:);
% gust2t= Wing_twist(Steps+1:Steps*2,:);
% gust3t= Wing_twist(Steps*2+1:Steps*3,:);
% gust4t= Wing_twist(Steps*3+1:Steps*4,:);
% gust5t= Wing_twist(Steps*4+1:Steps*5,:);
% gust6t= Wing_twist(Steps*5+1:Steps*6,:);
% gust7t= Wing_twist(Steps*6+1:Steps*7,:);
% gust8t= Wing_twist(Steps*7+1:Steps*8,:);
% gust9t= Wing_twist(Steps*8+1:Steps*9,:);
% gust10t = Wing_twist(Steps*9+1:Steps*10,:);
% 
% % curve = animatedline;
% 
% % 
% % for i=1:201
% %     clearpoints(curve)
% %     x=1:1:25;
% %     addpoints(curve,x,gust1(i,:))
% %     
% %    drawnow
% %     
% %     
% % end
% 
% x=1:1:25;
% y=gust6t;
% 
% figure 
% ph=plot(x,y(100,:),'b-s');
% set(gca,'Xlim',[0 25],'Ylim',[-0.01 0.01])
% set(gcf,'Color','w')
% xlabel('Wing span','Interpreter','latex')
% ylabel('Twist','Interpreter','latex')
% 
% for i=1:Steps
%     ph.XData=x;
%     ph.YData=y(i,:);
%     
%     drawnow
%     pause(0.05)
%     
% end    
%-------------------------------------------------------------------------

% [ind1,~]=find(Gust_data.ID==1396); % cog
% [ind21,~]=find(Gust_data.ID==1404); % wing monitor 1
% [ind22,~]=find(Gust_data.ID==1409); % wing monitor 2
% [ind23,~]=find(Gust_data.ID==1414); % wing monitor 3
% [ind2,~]=find(Gust_data.ID==1420); % wing tip 
% 
% temp1=Gust_data.Z([ind1]);
% temp21=Gust_data.Z([ind21]);
% temp22=Gust_data.Z([ind22]);
% temp23=Gust_data.Z([ind23]);
% temp2=Gust_data.Z([ind2]);
% 
% Z_disp_cog=reshape(temp1,[201,NumGust]);
% Z_disp_tip=reshape(temp2,[201,NumGust]);
% Z_disp_w1=reshape(temp21,[201,NumGust]); 
% Z_disp_w2=reshape(temp22,[201,NumGust]);
% Z_disp_w3=reshape(temp23,[201,NumGust]);
% 
% % group 2: fuselage
% % ---------------------------------------------------------
% [ind3,~]=find(Gust_data.ID==1021); % head of fuselage
% [ind4,~]=find(Gust_data.ID==1095); % tail of the fuselage
% 
% temp3=Gust_data.Z([ind3]);
% temp4=Gust_data.Z([ind4]);
% 
% Z_disp_head=reshape(temp3,[201,NumGust]);
% Z_disp_tail=reshape(temp4,[201,NumGust]);
% 
% % plot ----------------------------------------------------------------
% 
% figure  % cog displacement 
% plot(t,Z_disp_cog(:,1))
% hold on 
% plot(t,Z_disp_cog(:,2))
% hold on 
% plot(t,Z_disp_cog(:,3))
% hold on 
% plot(t,Z_disp_cog(:,4))
% hold on 
% plot(t,Z_disp_cog(:,5))
% hold on 
% plot(t,Z_disp_cog(:,6))
% hold on 
% plot(t,Z_disp_cog(:,7))
% hold on 
% plot(t,Z_disp_cog(:,8))
% 
% 
% figure  % Tip displacment 
% plot(t,Z_disp_tip(:,1))
% hold on 
% plot(t,Z_disp_tip(:,2))
% hold on 
% plot(t,Z_disp_tip(:,3))
% hold on 
% plot(t,Z_disp_tip(:,4))
% hold on 
% plot(t,Z_disp_tip(:,5))
% hold on 
% plot(t,Z_disp_tip(:,6))
% hold on 
% plot(t,Z_disp_tip(:,7))
% hold on 
% plot(t,Z_disp_tip(:,8))
% hold on 
% plot(t,Z_disp_tip(:,9))
% hold on 
% plot(t,Z_disp_tip(:,10))
% 
% figure % wing shape
% plot(t,Z_disp_cog(:,1))
% hold on 
% plot(t,Z_disp_w1(:,1))
% hold on 
% plot(t,Z_disp_w2(:,1))
% hold on 
% plot(t,Z_disp_w3(:,1))
% hold on 
% plot(t,Z_disp_tip(:,1))
% 
% 
% 
% 
% 
% figure 
% plot(t,Z_disp_cog(:,2))
% 
% hold on 
% plot(t,Z_disp_tip(:,2))
% % Wing tip twist
% 
% 
% % Root torque
% 
% 
% % potato plot


%% old code
    

%     GustLoadcase = awi.model.LoadCase;
%     GustLoadcase.LoadCaseType = 'Pratt Gust';
%     GustLoadcase.Altitude   = 3000;
%     GustLoadcase.AcVelocity = 0.70*340;
%     GustLoadcase.AcMass = 70000;
%     GustLoadcase.Mach = 0.7;
%     GustLoadcase.GustLength = [50,200];
% %     obj=GustLoadcase;
% 
% [a,b]=calculateGustLoadFactor(GustLoadcase,Aircraft);
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
%     NastranMethods1.writeGustFile(Aircraft, GustLoadcase, MassCases, FlightPoint, run_folder)
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
%         M_root10(i)=F06_corrected.f06data.Bendingmoment(i+201).UMPLN2(69); % out of plane moment 69 at root
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
% plot(t,M_root10,'r.')
% hold on 
% plot(t,M_root11,'b.')
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
 
 
 
%  % Result extract
%  F061=NastranMethods1.extractNastranResults('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\gust_test\guest_analysis.f06','ReadF06',true,'ReadHDF5',false);
%  
%  % extract forces
%  M_P1=[F061.f06data.Bendingmoment.UMPLN1(69:89),F061.f06data.Bendingmoment.LMPLN1(89)]; % in-plane moment
%  M_P2=[F061.f06data.Bendingmoment.UMPLN2(69:89),F061.f06data.Bendingmoment.LMPLN2(89)]; % out of plane moment
%  T=[F061.f06data.Bendingmoment.UTORQUE1(69:89),F061.f06data.Bendingmoment.LTORQUE1(89)];% torque
%  
%  S_P1=[F061.f06data.Bendingmoment.USPLN1(69:89),F061.f06data.Bendingmoment.LSPLN1(89)]; % in plane shear
%  S_P2=[F061.f06data.Bendingmoment.USPLN2(69:89),F061.f06data.Bendingmoment.LSPLN2(89)]; % out of plane shear
%  
%  
%  Grid_coord = h5read('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\NSM_test\A320_half_model_SOL144.h5','/NASTRAN/INPUT/NODE/GRID');
%  Displacement = h5read('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\NSM_test\A320_half_model_SOL144.h5','/NASTRAN/RESULT/NODAL/DISPLACEMENT');
%  
%  Y=Grid_coord.X(2,346:367);
%  Displacement_Z=Displacement.Z(346:367);
%  
%  figure % bending moment
%  plot(Y,M_P2,'b-s','MarkerFaceColor','b')
%  
%  xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%  ylabel('Bending moment (Nm)','FontSize',12,'Interpreter','latex')
%  set(gcf,'color','w')
%  
%  figure %shear force
%  plot(Y, S_P2,'b-s','MarkerFaceColor','b')
%  xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%  ylabel('Vertical shear force (N)','FontSize',12,'Interpreter','latex')
%  set(gcf,'color','w')
%  
%  figure % torque
%  plot(Y, T,'b-s','MarkerFaceColor','b')
%  xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%  ylabel('Torque (Nm)','FontSize',12,'Interpreter','latex')
%  set(gcf,'color','w')
 
 
 
 
 
     
     
     
     
    