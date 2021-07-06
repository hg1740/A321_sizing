

%% initial guesses

%A320
x=1.2*[0.0160369401546491,0.0161179748570384,0.0165944334426013,0.0175266294784001,0.0187239043216224,0.0207096617121334,0.0230583809562627,0.0320089957111805,0.0320392423264178,0.0319427876288856,0.0315322529714459,0.0322935890059832,0.0323900672122562,0.0323115468409948,0.0321955616764684,0.0313429255556067,0.0311508861821847,0.0304592682549718,0.0297595013116213,0.0272724891942158,0.0271122094889916,0.0264611574714236,0.0230559955692698,0.0187199645209936,0.0107524736198303,0.00647614132210356,0.00657526982363507,0.00669537365768939,0.00683776660015478,0.00701687360523511,0.00721883252923472,0.00748254587435518,0.00755456273799370,0.00737352569279649,0.00716564274351893,0.00696475875843593,0.00674100541947114,0.00649425980510992,0.00623148973761537,0.00595214725710100,0.00564842457162947,0.00541126353104570,0.00522201059650361,0.00497676988997949,0.00469043302797361,0.00427877225656423,0.00376751794116612,0.00307441811486333,0.00231764849271284,0.000789440488034985,5.19302425632857e-05,5.34268357423746e-05,5.52318694509354e-05,5.74180048238046e-05,6.02580769479868e-05,6.35124694224661e-05,6.77949943208388e-05,6.87713796834609e-05,6.55154629087215e-05,6.17954684314261e-05,5.82619158760843e-05,5.45303477290052e-05,5.06453591835854e-05,4.65755305830927e-05,4.24490132198509e-05,3.80255976131977e-05,3.49566602319797e-05,3.25217162959330e-05,2.94290956549030e-05,2.60051023480533e-05,2.16034838172783e-05,1.67146290536804e-05,1.13000939281217e-05,6.55690899275516e-06,7.23256842198466e-07];

%A321
x=[0.0189596521054333,0.0191837427745745,0.0198539169260718,0.0210684229892251,0.0226246281693717,0.0250916230994881,0.0280488941974237,0.0367041175372602,0.0367553345728376,0.0366612768523305,0.0362145073754782,0.0370868191031210,0.0372021678219872,0.0371206354501631,0.0369988383777017,0.0360401301612646,0.0361049634501639,0.0356335403745682,0.0348916331838040,0.0296388702926122,0.0286701277141633,0.0277597217308701,0.0237803010813388,0.0191787360923511,0.0121862598372589,0.00699956238223479,0.00709931837280801,0.00721779831240336,0.00735837047471829,0.00753590622003866,0.00773387497843092,0.00799253834763674,0.00813544045254821,0.00794177830530960,0.00771903722029878,0.00750366043193169,0.00726346792293347,0.00699820361649310,0.00671555090260363,0.00641531700369381,0.00608847909600772,0.00571105500399573,0.00537365397847822,0.00508789516567232,0.00481404400255434,0.00438627803345660,0.00385033301066460,0.00312019520069395,0.00235313284362065,0.000845505438493128,6.07225713152272e-05,6.23599038943347e-05,6.43148014859593e-05,6.66620770111903e-05,6.96900056705333e-05,7.31110690050333e-05,7.77506255979472e-05,8.00010789470476e-05,7.59042231708472e-05,7.15881361765266e-05,6.75093353966714e-05,6.31981619690458e-05,5.87027230664615e-05,5.39922987204662e-05,4.92202772665568e-05,4.40999407554106e-05,3.88105719549454e-05,3.42919819496080e-05,3.05707906441678e-05,2.71569808553286e-05,2.25104864554471e-05,1.72665450978065e-05,1.15467937795578e-05,6.73738284228341e-06,8.31725206912222e-07];
% x=1.5*[0.0191135905034748,0.0193953734978214,0.0201137604804691,0.0213823646563527,0.0229999860459657,0.0255251238127354,0.0285624381189545,0.0372113474624173,0.0373207109980618,0.0372799899057456,0.0368848456350822,0.0378039760454905,0.0379687298578435,0.0379303106594230,0.0378464169968381,0.0368959067148779,0.0365928771944430,0.0358740997358953,0.0351157988905984,0.0325204105754142,0.0321366355707771,0.0317223020101411,0.0280889429523174,0.0226448795839794,0.0118821391075936,0.00707914701666205,0.00718361444866246,0.00730636454511555,0.00745189886461910,0.00763438438081004,0.00783787156936675,0.00810313558972640,0.00823100835918603,0.00803824902557083,0.00781619282888881,0.00760074969073965,0.00736034765742680,0.00709415722939268,0.00680991625431645,0.00650766890180177,0.00618441788926231,0.00597655499191214,0.00574729769090399,0.00545697336175644,0.00516124669194006,0.00473319707846889,0.00418339089852517,0.00343190389276599,0.00256919715920123,0.000896999342517673,6.20362573091953e-05,6.37511021940462e-05,6.57911120554895e-05,6.82358984870429e-05,7.13734187868434e-05,7.49159107475837e-05,7.98695370113729e-05,8.20549947219400e-05,7.77731126376305e-05,7.33381910989789e-05,6.92092572133319e-05,6.48380207016889e-05,6.02644951660675e-05,5.54647083143445e-05,5.05953480091047e-05,4.54572672185587e-05,4.24462089148306e-05,3.91789917569685e-05,3.51963688050882e-05,3.13298871582112e-05,2.63163732081641e-05,2.05279963925629e-05,1.40184642076611e-05,8.01533883220779e-06,9.15451870145967e-07];
%% run_folder

    run_folder = [
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\Full_Ref_AC_Sizing_SF']; %[-], folder for exporting the NASTRAN model
    
%% Safty Factor

    SF=1.5;
    
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
    
    Mid_chord=0.63685*Root_chord;
    Tip_chord=0.2248*Root_chord;
    
    X0=Root_chord; 
    X1=0.27*Semi_span*tan(27*pi/180) + Mid_chord;
    X2=Semi_span*tan(27*pi/180) + Tip_chord;
    
    tan_TE_sweep1=(X1-X0)/(0.27*Semi_span);
    tan_TE_sweep2=(X2-X1)/(0.73*Semi_span);
    
    TE_sweep1=atan(tan_TE_sweep1)*180/pi; % deg
    TE_sweep2=atan(tan_TE_sweep2)*180/pi; % deg
      
 
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
    Wingbox.TESweep     = [TE_sweep1, TE_sweep2, TE_sweep2];
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
    
    draw(Wingbox)
  
    
    %% A321 Fuselage length, wing positions, engine position
    Fuselage_length=45;
    Wing_position=20;
    Horizontal_tail_position=42;
    Vertical_tail_position=41;
    Engine_position=4.29;
    
    %% A320 Fuselage length and wing positions
%     Fuselage_length=38;
%     Wing_position=15;
%     Horizontal_tail_position=35;
%     Vertical_tail_position=34;
%     Engine_position=4.29;       

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
        
        [Wingbox_mass, Secondary_mass, Total_mass]=Wing_masses_v1(x,Height_var,Width_var,etaRL, MTOW, Wing_area, LE_sweep, Semi_span);
        
        Guess_wing_mass=Total_mass;
        
        Fuselage_structure_mass=OWE - 2*Guess_wing_mass - Horizontal_tail - Vertical_tail -Engine_mass*2 - Pylon*2 - Fuselage_shell_mass;
        
        [Fuel_mass,Fuselage_total_mass]= Lumped_masses(Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
        
      
        % size loop
        %initial condition
        cond_set1=[1.1,1.1,1.1,1.1];
        cond_set2=[1.1,1.1,1.1,1.1];
        
        Numsec=25;
        
        record=zeros(100,75);
        
        
        while max(cond_set1)>1.01 || min(cond_set2)<0.992
            
            record(counter,:)=x(1,:);
                
            [Aircraft,FEM_full,Y_Data, Load_distribution, Delta, RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=Stress_Calc_LoadCases_v3(run_folder,x,Semi_span,Root_chord,LE_sweep,TE_sweep1,TE_sweep2, BeamLoc,Fuselage_length, Wing_position, Engine_position, Horizontal_tail_position, Vertical_tail_position, Secondary_mass,Fuel_mass,Fuselage_total_mass);
          
            % skin check
            Check_von_skin1=SF*RVon_skn/Yield_strength;
            Check_von_skin2=SF*Rsigma_pp./Rsigmab_skn;
            
            V_skin_change=step_size(Check_von_skin1);
            B_skin_change=step_size(Check_von_skin2);
            
            skin_coefficient=max([V_skin_change;B_skin_change]);
                
            x(1,26:50)=x(1,26:50).*skin_coefficient;

            % spar check
            Check_von_spar=SF*RVon_spr/Yield_strength;
            
            Spar_coefficient=step_size(Check_von_spar);
           
            x(1,1:25)=x(1,1:25).*Spar_coefficient;
            
            % stringers check
            Check_strg=SF*Rsigma_strg./Rsigma_col;
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
        
        
        [Wingbox_mass, Secondary_mass, Total_mass]=Wing_masses_v1(x,Height_var,Width_var,etaRL,MTOW,Wing_area,LE_sweep,Semi_span);
        
        [Fuel_mass,Fuselage_total_mass]= Lumped_masses(Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
        
        Wing_total_mass=Wingbox_mass+Secondary_mass;
        
        Aircraft_total_mass=Fuselage_total_mass+Fuselage_shell_mass+Wing_total_mass*2+Horizontal_tail+ Vertical_tail + Engine_mass*2 + Pylon*2 + Fuel_mass*2;
        
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
    
    
    
    


