

% initial guess

% x=[0.0184549632155560,0.0194332413084522,0.0207559132418991,0.0209074925629807,0.0231661967456849,0.0254828164202534,0.0264270050000000,0.0378563183137105,0.0384112626298514,0.0373093915208931,0.0364906994983588,0.0389743420000000,0.0395598526378275,0.0375818600059361,0.0370256249000000,0.0364906994983588,0.0362261419269957,0.0349193296635013,0.0362261419269957,0.0351743436550000,0.0305996900000000,0.0252980660012066,0.0153689282363496,0.000576499130128597,0.000697563947455602,0.00690405505625000,0.00700777451556970,0.00711050306418304,0.00726742637500000,0.00748474006756110,0.00770855196712667,0.00805255000000000,0.00811426522855439,0.00793621128715938,0.00770855196712667,0.00753940072280141,0.00726742637500000,0.00705895191696771,0.00675497791097389,0.00651130062423758,0.00614088900997626,0.00583384455947745,0.00542053909375000,0.00500179247918198,0.00444730593828125,0.00384085512851562,0.00319858816010073,0.00223369340704222,0.00180616690062519,0.00191511288486283,5.93348577610401e-05,5.93348577610401e-05,6.20049263602869e-05,6.47951480464998e-05,6.82054189963156e-05,7.17951778908585e-05,7.89746956799444e-05,7.84021291362648e-05,7.44820226794515e-05,7.12746628511498e-05,6.52683435371441e-05,6.20049263602869e-05,5.93348577610401e-05,5.39407797827637e-05,4.90370725297852e-05,4.45791568452593e-05,3.85001809118148e-05,3.34929803495562e-05,2.76801490492200e-05,2.17323484270735e-05,1.56247225182875e-05,9.21663720956289e-06,3.02451924944460e-06,2.23953621083449e-06,2.49960268549141e-06];
% x=[0.023166197	0.023514222	0.023858922	0.025298066	0.027626121	0.03061066	0.033173363	0.043200348	0.043833632	0.044492143	0.043833632	0.044476199	0.044817067	0.045144364	0.044476199	0.043833632	0.042887146	0.043200348	0.042887146	0.037856318	0.032932856	0.030168414	0.024217354	0.018488067	0.010905653	0.007594461	0.007708552	0.007821553	0.007994169	0.008173523	0.008417931	0.008729832	0.00886098	0.008666541	0.008417931	0.008233214	0.007936211	0.007652665	0.007376605	0.007058952	0.006706004	0.006278664	0.005833845	0.005344141	0.004856569	0.004415063	0.003842232	0.003061947	0.002251619	0.000934287	7.18E-05	7.18E-05	7.50E-05	7.84E-05	8.19E-05	9.08E-05	9.49E-05	9.91E-05	9.42E-05	8.56E-05	7.84E-05	7.45E-05	7.08E-05	6.48E-05	5.85E-05	5.35E-05	4.76E-05	3.97E-05	3.28E-05	2.73E-05	2.25E-05	1.69E-05	1.14E-05	6.07E-06	1.02E-06
% ];
% x=[0.023166197	0.023514222	0.023858922	0.025298066	0.027626121	0.03061066	0.033173363	0.043200348	0.043833632	0.044492143	0.043833632	0.044476199	0.044817067	0.045144364	0.044476199	0.043833632	0.042887146	0.043200348	0.042887146	0.037856318	0.032932856	0.030168414	0.024217354	0.018488067	0.010905653	0.007594461	0.007708552	0.007821553	0.007994169	0.008173523	0.008417931	0.008729832	0.00886098	0.008666541	0.008417931	0.008233214	0.007936211	0.007652665	0.007376605	0.007058952	0.006706004	0.006278664	0.005833845	0.005344141	0.004856569	0.004415063	0.003842232	0.003061947	0.002251619	0.000934287	7.18E-05	7.18E-05	7.50E-05	7.84E-05	8.19E-05	9.08E-05	9.49E-05	9.91E-05	9.42E-05	8.56E-05	7.84E-05	7.45E-05	7.08E-05	6.48E-05	5.85E-05	5.35E-05	4.76E-05	3.97E-05	3.28E-05	2.73E-05	2.25E-05	1.69E-05	1.14E-05	6.07E-06	1.02E-06
% ];

% x=[0.021828831	0.022111285	0.022887508	0.0242794	0.02605858	0.0288485	0.032216892	0.041899822	0.0419988	0.041941755	0.041493825	0.04249848	0.042657263	0.042605098	0.042517685	0.041506143	0.041608794	0.041304963	0.041016216	0.035575335	0.030881665	0.028483845	0.023391906	0.01752022	0.010532432	0.007541687	0.007649915	0.007777495	0.007929487	0.008121036	0.008333552	0.008611974	0.008776701	0.008573593	0.008339977	0.008113313	0.007859807	0.007579299	0.007280153	0.006963213	0.006616882	0.006217132	0.005789641	0.005285283	0.004872181	0.004420193	0.003847862	0.003067052	0.002247143	0.000925715	7.04E-05	7.22E-05	7.45E-05	7.74E-05	8.13E-05	8.63E-05	9.37E-05	9.77E-05	9.11E-05	8.44E-05	7.87E-05	7.38E-05	6.86E-05	6.33E-05	5.78E-05	5.19E-05	4.59E-05	3.97E-05	3.28E-05	2.75E-05	2.26E-05	1.70E-05	1.11E-05	6.15E-06	1.02E-06
% ];


%A320
% x=[0.0160369401546491,0.0161179748570384,0.0165944334426013,0.0175266294784001,0.0187239043216224,0.0207096617121334,0.0230583809562627,0.0320089957111805,0.0320392423264178,0.0319427876288856,0.0315322529714459,0.0322935890059832,0.0323900672122562,0.0323115468409948,0.0321955616764684,0.0313429255556067,0.0311508861821847,0.0304592682549718,0.0297595013116213,0.0272724891942158,0.0271122094889916,0.0264611574714236,0.0230559955692698,0.0187199645209936,0.0107524736198303,0.00647614132210356,0.00657526982363507,0.00669537365768939,0.00683776660015478,0.00701687360523511,0.00721883252923472,0.00748254587435518,0.00755456273799370,0.00737352569279649,0.00716564274351893,0.00696475875843593,0.00674100541947114,0.00649425980510992,0.00623148973761537,0.00595214725710100,0.00564842457162947,0.00541126353104570,0.00522201059650361,0.00497676988997949,0.00469043302797361,0.00427877225656423,0.00376751794116612,0.00307441811486333,0.00231764849271284,0.000789440488034985,5.19302425632857e-05,5.34268357423746e-05,5.52318694509354e-05,5.74180048238046e-05,6.02580769479868e-05,6.35124694224661e-05,6.77949943208388e-05,6.87713796834609e-05,6.55154629087215e-05,6.17954684314261e-05,5.82619158760843e-05,5.45303477290052e-05,5.06453591835854e-05,4.65755305830927e-05,4.24490132198509e-05,3.80255976131977e-05,3.49566602319797e-05,3.25217162959330e-05,2.94290956549030e-05,2.60051023480533e-05,2.16034838172783e-05,1.67146290536804e-05,1.13000939281217e-05,6.55690899275516e-06,7.23256842198466e-07];


%A321
x=[0.0189596521054333,0.0191837427745745,0.0198539169260718,0.0210684229892251,0.0226246281693717,0.0250916230994881,0.0280488941974237,0.0367041175372602,0.0367553345728376,0.0366612768523305,0.0362145073754782,0.0370868191031210,0.0372021678219872,0.0371206354501631,0.0369988383777017,0.0360401301612646,0.0361049634501639,0.0356335403745682,0.0348916331838040,0.0296388702926122,0.0286701277141633,0.0277597217308701,0.0237803010813388,0.0191787360923511,0.0121862598372589,0.00699956238223479,0.00709931837280801,0.00721779831240336,0.00735837047471829,0.00753590622003866,0.00773387497843092,0.00799253834763674,0.00813544045254821,0.00794177830530960,0.00771903722029878,0.00750366043193169,0.00726346792293347,0.00699820361649310,0.00671555090260363,0.00641531700369381,0.00608847909600772,0.00571105500399573,0.00537365397847822,0.00508789516567232,0.00481404400255434,0.00438627803345660,0.00385033301066460,0.00312019520069395,0.00235313284362065,0.000845505438493128,6.07225713152272e-05,6.23599038943347e-05,6.43148014859593e-05,6.66620770111903e-05,6.96900056705333e-05,7.31110690050333e-05,7.77506255979472e-05,8.00010789470476e-05,7.59042231708472e-05,7.15881361765266e-05,6.75093353966714e-05,6.31981619690458e-05,5.87027230664615e-05,5.39922987204662e-05,4.92202772665568e-05,4.40999407554106e-05,3.88105719549454e-05,3.42919819496080e-05,3.05707906441678e-05,2.71569808553286e-05,2.25104864554471e-05,1.72665450978065e-05,1.15467937795578e-05,6.73738284228341e-06,8.31725206912222e-07];

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
    
    
    %% obtain wingbox geometric properties 

    Wingbox = awi.model.LiftingSurface;
    Wingbox.Origin=[20,2,0];

    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox.ActiveSet = 'sSet';

    % Num of element
    Wingbox.NumBeamElem = 23;

    %Wing dimensions
    Wingbox.SpanVector  = 'Y';
    Wingbox.Span        = Semi_span;   %34.1/2;
    Wingbox.LESweep     = [LE_sweep, LE_sweep];
    Wingbox.LESweep_eta = [0, 1];
    Wingbox.TESweep     = [0, TE_sweep, TE_sweep];
    Wingbox.TESweep_eta = [0, 0.27, 1];
    Wingbox.RootChord   = Root_chord;
    
    build(Wingbox)

    NumSec=Wingbox.NumBeamElem+2;

    YData=Wingbox.YData;
    SparWidth=Wingbox.Chord*0.5;

    RootH=Wingbox.Chord(1)*0.15; % root thickness/chord = 0.15
    MidH=Wingbox.Chord(2)*0.12;  % middle thickness/chord = 0.12
    TipH=Wingbox.Chord(end)*0.11;% tip thickness/chord = 0.11


    % set up eta values
    elnum=Wingbox.NumBeamElem + 1; % total number of beam elements along the wing
    Num_seg1=ceil(elnum*0.27); % number of elements in the inboard section
    Num_seg2=elnum - Num_seg1; % number of elements in the outboard section

    Num_sec1=Num_seg1+1;    % number of sections in the inboard section
    Num_sec2=Num_seg2+1;    % number of sections in the outboard section

    eta1_=linspace(0,0.27, Num_sec1);
    eta2_=linspace(0.27,1,Num_sec2);
    etaS=[eta1_(1:end-1),eta2_(1:end)];

    RData=Wingbox.RData;
    eta_R=RData/RData(end);
    eta_Y=YData/YData(end);
    etaRS=interp1(eta_Y,eta_R,etaS);
    
    etaRL=etaRS*RData(end);

    Width_var=interp1(RData/RData(end),SparWidth,etaRS);
    Height_var=interp1(RData/RData(end),0.79*[RootH,MidH,TipH],etaRS);
  
    
    %% A321 Fuselage length and wing positions 
    Fuselage_length=45;
    Wing_position=20;
    Horizontal_tail_position=42;
    Vertical_tail_position=41;
    
    %% A320 Fuselage length and wing positions
%     Fuselage_length=38;
%     Wing_position=15;
%     Horizontal_tail_position=35;
%     Vertical_tail_position=34;

  %% A321 mass configurations 

    Payload_max=25000; % kg
    Fuel_fraction=0.723; % percentage of fuel in the tank
    
    Fuel_capacity=32940; % L

    MTOW=93500; % maximum take off mass
    OWE=48500;  % Operating empty mass
    MWE=44057;  % Manufacture's empty mass
    MZF=73000;  % Maxumum zero fuel mass
    Fuselage_shell_mass= 2*pi*2*0.004*2800*44.5;

  
    Engine_mass=7362/2; % kg
    Pylon=1239/2; % kg
    Horizontal_tail=682; % kg
    Vertical_tail=522; % kg
    
    Guess_wing_weight=0.86*Wing_span*(Total_area*MZF*MTOW)^0.25; % give an initial guess of wing box weight
    Guess_wing_mass=Guess_wing_weight/10; % give an initial guess of wing box weight. For A320, the wing mass = 8097 kg
    
    Fuselage_structure_mass=OWE - 2*Guess_wing_mass - Horizontal_tail - Vertical_tail -Engine_mass*2 - Pylon*2-Fuselage_shell_mass;
    
%% A320 mass configurations

%     Payload_max=19900; % kg
%     Fuel_fraction=0.56; % percentage of fuel in the tank
%     
%     Fuel_capacity=32940; % L
% 
%     MTOW=78000; % maximum take off mass
%     OWE=42600;  % Operating empty mass
%     MWE=44057;  % Manufacture's empty mass
%     MZF=66000;  % Maxumum zero fuel mass
% 
%   
%     Engine_mass=7362/2; % kg
%     Pylon=1239/2; % kg
%     Horizontal_tail=682; % kg
%     Vertical_tail=522; % kg
%     Fuselage_shell_mass= 2*pi*2*0.004*2800*37.5;
%     
%     Guess_wing_weight=0.86*Wing_span*(Total_area*MZF*MTOW)^0.25; % give an initial guess of wing box weight
%     Guess_wing_mass=Guess_wing_weight/10; % give an initial guess of wing box weight. For A320, the wing mass = 8097 kg
%     
%     Fuselage_structure_mass=OWE - 2*Guess_wing_mass - Horizontal_tail - Vertical_tail -Engine_mass*2 - Pylon*2-Fuselage_shell_mass;
    
    %% Material properties : Al7075
    
    Yield_strength = 5.2e8;
      
    %% start the loop
 
    Wing_total_mass=100; % Triger the sizing loop
       
    counter=1;
    
    % mass loop
    while abs(Wing_total_mass/Guess_wing_mass - 1)>=0.05
        
        Wingbox_mass=Wing_masses_v1(x,Height_var,Width_var,etaRL);
        
        [Secondary_mass,~,~]= Lumped_masses(MTOW,Wing_area,LE_sweep,Semi_span,Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
        
        Guess_wing_mass=Wingbox_mass+Secondary_mass;
        
        Fuselage_structure_mass=OWE - 2*Guess_wing_mass - Horizontal_tail - Vertical_tail -Engine_mass*2 - Pylon*2 - Fuselage_shell_mass;
        
        [Secondary_mass,Fuel_mass,Fuselage_total_mass]= Lumped_masses(MTOW,Wing_area,LE_sweep,Semi_span,Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
        
      
        % size loop
        %initial condition
        cond_set1=[1.1,1.1,1.1,1.1];
        cond_set2=[1.1,1.1,1.1,1.1];
        
        Numsec=25;
        
        record=zeros(100,75);
        
        
        while max(cond_set1)>1.005 || min(cond_set2)<0.995
            
            record(counter,:)=x(1,:);
                
            [RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=Stress_Calc_LoadCases_v2(x,Semi_span,Root_chord,LE_sweep,TE_sweep,BeamLoc,Fuselage_length, Wing_position, Horizontal_tail_position, Vertical_tail_position, Secondary_mass,Fuel_mass,Fuselage_total_mass);

%             run_folder = [
%              'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\full_A321_sizing_test4']; %[-], folder for exporting the NASTRAN model
%          
%          Engine_position=4.2;
%          
% [Connector_right, Connector_left, Wingbox_right, Wingbox_left, Tailwing_right, Tailwing_left, Verticalwing, Aircraft, FEM_full]=Aircraft_Model(run_folder,x,Semi_span,Root_chord,LE_sweep,TE_sweep,BeamLoc,Fuselage_length, Wing_position, Engine_position, Horizontal_tail_position, Vertical_tail_position, Secondary_mass,Fuel_mass,Fuselage_total_mass)

            % skin check
            Check_von_skin1=RVon_skn/Yield_strength;
            Check_von_skin2=Rsigma_pp./Rsigmab_skn;
            
            V_skin_change=step_size(Check_von_skin1);
            B_skin_change=step_size(Check_von_skin2);
            
            skin_coefficient=max([V_skin_change;B_skin_change]);
                
            x(1,26:50)=x(1,26:50).*skin_coefficient;

            % spar check
            Check_von_spar=RVon_spr/Yield_strength;
            
            Spar_coefficient=step_size(Check_von_spar);
           
            x(1,1:25)=x(1,1:25).*Spar_coefficient;
            
            % stringers check
            Check_strg=Rsigma_strg./Rsigma_col;
            Stringer_coefficient=step_size(Check_strg);
            
            x(1,51:75)=x(1,51:75).*Stringer_coefficient;

            % end condition check
            cond1=max(Check_von_skin1);
            cond2=max(Check_von_skin2);
            cond3=max(Check_von_spar);
            cond4=max(Check_strg);
            
            cond5=min(Check_von_skin1);
            cond6=min(Check_von_skin2);
            cond7=min(Check_von_spar);
            cond8=min(Check_strg);
            
            cond_set1=[cond1,cond2,cond3,cond4];
            cond_set2=[cond6,cond7,cond8];
            
            disp(cond_set1);
            disp(cond_set2);
            
            counter=counter+1;
            
            
        end
        
        
        Wingbox_mass=Wing_masses_v1(x,Height_var,Width_var,etaRL);
        
        [Secondary_mass,Fuel_mass,Fuselage_total_mass]= Lumped_masses(MTOW,Wing_area,LE_sweep,Semi_span,Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
        
        Wing_total_mass=Wingbox_mass+Secondary_mass;
        
        progress_indicator=abs(Wing_total_mass/Guess_wing_mass - 1);
        
        disp(progress_indicator);
        

        
    end
    
    % obtain lumped masses: Secondar_mass and fuel mass are the masses for
    % each wingbox. Fuselage_total_mass includs the structure mass and
    % payload
    
    % the sizing is carroied out using max. payload with its corresponding
    % max fuel mass 
    
    %% figure 
    
%     NastranMethods1 = awi.methods.Nastran;
%     NastranMethods1.AnalysisModel = FEM_full;
%     
%     run_folder = [
%         'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\full_A320_sizing_test2']; %[-], folder for exporting the NASTRAN model
%     
%     WingNodes=NastranMethods1.WingNode;
%     
%     Y_data = h5read(strcat(run_folder,'\gust_analysis.h5'),'/NASTRAN/INPUT/NODE/GRID');
%     
%     Yindex=ismember([Y_data.ID],[WingNodes(1:end).GID]);
%     
%     Y=Y_data.X(2,Yindex)-2;
%     
%     figure % thickness1: spar
%     plot(Y, x(1:25),'bs','MarkerFaceColor','b')
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Spar thickness (m)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
%     
%     figure % thickness1: skin
%     plot(Y, x(26:50),'bs','MarkerFaceColor','b')
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Skin thickness (m)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
%     
%     figure % area: stringers
%     plot(Y, x(51:75),'bs','MarkerFaceColor','b')
%     xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
%     ylabel('Stringers Area (m$^2$)','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
%     
%     figure
%     rt=linspace(0.1,10,100);
%     delta=step_size(rt);
%     
%     loglog(rt,delta,'b-')
%     xlabel('$\frac{\sigma}{\sigma_c}$','FontSize',18,'Interpreter','latex')
%     ylabel('Thickness correction facor','FontSize',12,'Interpreter','latex')
%     set(gcf,'color','w')
%     axis([0.1 10 0.2 2])
    
    
    
    


