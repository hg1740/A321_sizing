

% initial guess

% x=[0.0184549632155560,0.0194332413084522,0.0207559132418991,0.0209074925629807,0.0231661967456849,0.0254828164202534,0.0264270050000000,0.0378563183137105,0.0384112626298514,0.0373093915208931,0.0364906994983588,0.0389743420000000,0.0395598526378275,0.0375818600059361,0.0370256249000000,0.0364906994983588,0.0362261419269957,0.0349193296635013,0.0362261419269957,0.0351743436550000,0.0305996900000000,0.0252980660012066,0.0153689282363496,0.000576499130128597,0.000697563947455602,0.00690405505625000,0.00700777451556970,0.00711050306418304,0.00726742637500000,0.00748474006756110,0.00770855196712667,0.00805255000000000,0.00811426522855439,0.00793621128715938,0.00770855196712667,0.00753940072280141,0.00726742637500000,0.00705895191696771,0.00675497791097389,0.00651130062423758,0.00614088900997626,0.00583384455947745,0.00542053909375000,0.00500179247918198,0.00444730593828125,0.00384085512851562,0.00319858816010073,0.00223369340704222,0.00180616690062519,0.00191511288486283,5.93348577610401e-05,5.93348577610401e-05,6.20049263602869e-05,6.47951480464998e-05,6.82054189963156e-05,7.17951778908585e-05,7.89746956799444e-05,7.84021291362648e-05,7.44820226794515e-05,7.12746628511498e-05,6.52683435371441e-05,6.20049263602869e-05,5.93348577610401e-05,5.39407797827637e-05,4.90370725297852e-05,4.45791568452593e-05,3.85001809118148e-05,3.34929803495562e-05,2.76801490492200e-05,2.17323484270735e-05,1.56247225182875e-05,9.21663720956289e-06,3.02451924944460e-06,2.23953621083449e-06,2.49960268549141e-06];
% x=[0.023166197	0.023514222	0.023858922	0.025298066	0.027626121	0.03061066	0.033173363	0.043200348	0.043833632	0.044492143	0.043833632	0.044476199	0.044817067	0.045144364	0.044476199	0.043833632	0.042887146	0.043200348	0.042887146	0.037856318	0.032932856	0.030168414	0.024217354	0.018488067	0.010905653	0.007594461	0.007708552	0.007821553	0.007994169	0.008173523	0.008417931	0.008729832	0.00886098	0.008666541	0.008417931	0.008233214	0.007936211	0.007652665	0.007376605	0.007058952	0.006706004	0.006278664	0.005833845	0.005344141	0.004856569	0.004415063	0.003842232	0.003061947	0.002251619	0.000934287	7.18E-05	7.18E-05	7.50E-05	7.84E-05	8.19E-05	9.08E-05	9.49E-05	9.91E-05	9.42E-05	8.56E-05	7.84E-05	7.45E-05	7.08E-05	6.48E-05	5.85E-05	5.35E-05	4.76E-05	3.97E-05	3.28E-05	2.73E-05	2.25E-05	1.69E-05	1.14E-05	6.07E-06	1.02E-06
% ];
x=[0.023166197	0.023514222	0.023858922	0.025298066	0.027626121	0.03061066	0.033173363	0.043200348	0.043833632	0.044492143	0.043833632	0.044476199	0.044817067	0.045144364	0.044476199	0.043833632	0.042887146	0.043200348	0.042887146	0.037856318	0.032932856	0.030168414	0.024217354	0.018488067	0.010905653	0.007594461	0.007708552	0.007821553	0.007994169	0.008173523	0.008417931	0.008729832	0.00886098	0.008666541	0.008417931	0.008233214	0.007936211	0.007652665	0.007376605	0.007058952	0.006706004	0.006278664	0.005833845	0.005344141	0.004856569	0.004415063	0.003842232	0.003061947	0.002251619	0.000934287	7.18E-05	7.18E-05	7.50E-05	7.84E-05	8.19E-05	9.08E-05	9.49E-05	9.91E-05	9.42E-05	8.56E-05	7.84E-05	7.45E-05	7.08E-05	6.48E-05	5.85E-05	5.35E-05	4.76E-05	3.97E-05	3.28E-05	2.73E-05	2.25E-05	1.69E-05	1.14E-05	6.07E-06	1.02E-06
];

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
    
    
  %% Mass configurations

    Payload_max=25000; % kg
    Fuel_fraction=0.9; % percentage of fuel in the tank
    
    Fuel_capacity=25000; % L

    MTOW=93500; % maximum take off mass
    OWE=48500;  % Operating empty mass
    MWE=44057;  % Manufacture's empty mass
    MZF=73000;  % Maxumum zero fuel mass
  
    Engine_mass=7362/2; % kg
    Pylon=1239/2; % kg
    Horizontal_tail=682; % kg
    Vertical_tail=522; % kg
    
    Guess_wing_weight=0.86*Wing_span*(Total_area*MZF*MTOW)^0.25; % give an initial guess of wing box weight
    Guess_wing_mass=1.1*Guess_wing_weight/10; % give an initial guess of wing box weight. For A320, the wing mass = 8097 kg
    
    %% material property : Al 7075
    
    Yield_Strength=5.2e8;
      
    %% start the loop
    
    Wingbox_mass=Wing_masses_v1(x,Height_var,Width_var,etaRL);
    
    [Secondary_mass,Fuel_mass,Fuselage_total_mass]= Lumped_masses(MTOW,Wing_area,LE_sweep,Semi_span,Fuel_fraction,Payload_max,Fuselage_structure_mass);
    
    Wing_total_mass=Wingbox_mass+Secondary_mass;
       
    counter=1;
    
    % mass loop
    while abs(Wing_total_mass/Guess_wing_mass - 1)>=0.05
        
        Wingbox_mass=Wing_masses_v1(x,Height_var,Width_var,etaRL);
        
        [Secondary_mass,Fuel_mass,Fuselage_total_mass]= Lumped_masses(MTOW,Wing_area,LE_sweep,Semi_span,Fuel_fraction,Payload_max,Fuselage_structure_mass);
        
        Guess_wing_mass=Wingbox_mass+Secondary_mass;
        
        Fuselage_structure_mass=OWE - Guess_wing_mass - Horizontal_tail - Vertical_tail;
       
        % size loop
        %initial condition
        cond_set1=[1.1,1.1,1.1,1.1];
        cond_set2=[1.1,1.1,1.1,1.1];
        
        Numsec=25;
        
        record=zeros(100,75);
        
        
        while max(cond_set1)>1 || min(cond_set2)<0.95
            
            increase_co=1.1;
            decrease_co=0.95;
            record(counter,:)=x(1,:);
                
            [RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=Stress_Calc_LoadCases_v1(x,Semi_span,Root_chord,LE_sweep,TE_sweep,BeamLoc,Secondary_mass,Fuel_mass,Fuselage_total_mass);
            
            % skin check
            Check_von_skin1=RVon_skn/Yield_Strength;
            Check_von_skin2=Rsigma_pp./Rsigmab_skn;
            
            [~,index1]=find(Check_von_skin1>1);
            [~,index2]=find(Check_von_skin1<0.95);
            [~,index3]=find(Check_von_skin2>1);
            [~,index4]=find(Check_von_skin2<0.95);
            
            skin_inc=unique([index1,index3]); % one of the condition suggest increase then need to increase t
            
            skin_dec=intersect(index2,index4); % only both condition suggest to decrease then decrease
            
            x(1,Numsec+skin_inc)=x(1,Numsec+skin_inc)*increase_co;
            x(1,Numsec+skin_dec)=x(1,Numsec+skin_dec)*decrease_co;
            
            % spar check
            Check_von_spar=RVon_spr/Yield_Strength;
            
            [~,sp_index1]=find(Check_von_spar>1);
            [~,sp_index2]=find(Check_von_spar<0.95);
            
            spar_inc=sp_index1;
            spar_dec=sp_index2;
            
            
            x(1,spar_inc)=x(1,spar_inc)*increase_co;
            x(1,spar_dec)=x(1,spar_dec)*decrease_co;
            
            % stringers check
            Check_strg=Rsigma_strg./Rsigma_col;
            
            [~,sg_index1]=find(Check_strg>1);
            [~,sg_index2]=find(Check_strg<0.95);
            
            strg_inc=sg_index1;
            strg_dec=sg_index2;
            
            
            x(1,Numsec*2+strg_inc)=x(1,Numsec*2+strg_inc)*increase_co;
            x(1,Numsec*2+strg_dec)=x(1,Numsec*2+strg_dec)*decrease_co;
            
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
        
        [Secondary_mass,Fuel_mass,Fuselage_total_mass]= Lumped_masses(MTOW,Wing_area,LE_sweep,Semi_span,Fuel_fraction,Payload_max,Fuselage_structure_mass);
        
        Wing_total_mass=Wingbox_mass+Secondary_mass;
        
        progress_indicator=abs(Wing_total_mass/Guess_wing_mass - 1);
        
        disp(progress_indicator);
        

        
    end
    
    % obtain lumped masses: Secondar_mass and fuel mass are the masses for
    % each wingbox. Fuselage_total_mass includs the structure mass and
    % payload
    
    % the sizing is carroied out using max. payload with its corresponding
    % max fuel mass 
    


