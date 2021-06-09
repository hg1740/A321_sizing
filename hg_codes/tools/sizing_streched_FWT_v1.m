

%% initial guesses

%A320
% x=[0.0160369401546491,0.0161179748570384,0.0165944334426013,0.0175266294784001,0.0187239043216224,0.0207096617121334,0.0230583809562627,0.0320089957111805,0.0320392423264178,0.0319427876288856,0.0315322529714459,0.0322935890059832,0.0323900672122562,0.0323115468409948,0.0321955616764684,0.0313429255556067,0.0311508861821847,0.0304592682549718,0.0297595013116213,0.0272724891942158,0.0271122094889916,0.0264611574714236,0.0230559955692698,0.0187199645209936,0.0107524736198303,0.00647614132210356,0.00657526982363507,0.00669537365768939,0.00683776660015478,0.00701687360523511,0.00721883252923472,0.00748254587435518,0.00755456273799370,0.00737352569279649,0.00716564274351893,0.00696475875843593,0.00674100541947114,0.00649425980510992,0.00623148973761537,0.00595214725710100,0.00564842457162947,0.00541126353104570,0.00522201059650361,0.00497676988997949,0.00469043302797361,0.00427877225656423,0.00376751794116612,0.00307441811486333,0.00231764849271284,0.000789440488034985,5.19302425632857e-05,5.34268357423746e-05,5.52318694509354e-05,5.74180048238046e-05,6.02580769479868e-05,6.35124694224661e-05,6.77949943208388e-05,6.87713796834609e-05,6.55154629087215e-05,6.17954684314261e-05,5.82619158760843e-05,5.45303477290052e-05,5.06453591835854e-05,4.65755305830927e-05,4.24490132198509e-05,3.80255976131977e-05,3.49566602319797e-05,3.25217162959330e-05,2.94290956549030e-05,2.60051023480533e-05,2.16034838172783e-05,1.67146290536804e-05,1.13000939281217e-05,6.55690899275516e-06,7.23256842198466e-07];

%A321
% x=[0.0191135905034748,0.0193953734978214,0.0201137604804691,0.0213823646563527,0.0229999860459657,0.0255251238127354,0.0285624381189545,0.0372113474624173,0.0373207109980618,0.0372799899057456,0.0368848456350822,0.0378039760454905,0.0379687298578435,0.0379303106594230,0.0378464169968381,0.0368959067148779,0.0365928771944430,0.0358740997358953,0.0351157988905984,0.0325204105754142,0.0321366355707771,0.0317223020101411,0.0280889429523174,0.0226448795839794,0.0118821391075936,0.00707914701666205,0.00718361444866246,0.00730636454511555,0.00745189886461910,0.00763438438081004,0.00783787156936675,0.00810313558972640,0.00823100835918603,0.00803824902557083,0.00781619282888881,0.00760074969073965,0.00736034765742680,0.00709415722939268,0.00680991625431645,0.00650766890180177,0.00618441788926231,0.00597655499191214,0.00574729769090399,0.00545697336175644,0.00516124669194006,0.00473319707846889,0.00418339089852517,0.00343190389276599,0.00256919715920123,0.000896999342517673,6.20362573091953e-05,6.37511021940462e-05,6.57911120554895e-05,6.82358984870429e-05,7.13734187868434e-05,7.49159107475837e-05,7.98695370113729e-05,8.20549947219400e-05,7.77731126376305e-05,7.33381910989789e-05,6.92092572133319e-05,6.48380207016889e-05,6.02644951660675e-05,5.54647083143445e-05,5.05953480091047e-05,4.54572672185587e-05,4.24462089148306e-05,3.91789917569685e-05,3.51963688050882e-05,3.13298871582112e-05,2.63163732081641e-05,2.05279963925629e-05,1.40184642076611e-05,8.01533883220779e-06,9.15451870145967e-07];

% x=[0.0263919314434434,0.0270287037846719,0.0283640513655332,0.0290476918785440,0.0305000319425819,0.0311454414799324,0.0396804988555224,0.0425444017731640,0.0441525651871870,0.0482605196982675,0.0464068437203061,0.0457224486348197,0.0456279338158660,0.0448628529788929,0.0435569212840316,0.0417731169704868,0.0396815956770828,0.0380480440213013,0.0358316633566355,0.0339687923155953,0.0337268684564342,0.0323902190802610,0.0291451850232813,0.0234232299918240,0.0191524074792969,0.00893809271271453,0.00898021206205229,0.00902935314915270,0.00908668154345703,0.00917673153841885,0.00928838436391813,0.00921342865699722,0.00927695325733132,0.00935527519615850,0.00946264343682268,0.00916930010549819,0.00885042818227504,0.00850732408723034,0.00814899922819781,0.00777157655127306,0.00739408745863920,0.00696700999391705,0.00650587525589280,0.00613500379948089,0.00572350270401927,0.00518759260675683,0.00454772718759639,0.00375606295531985,0.00275656826201626,0.00207551938332536,0.000107769528835642,0.000109277971748372,0.000111098275638144,0.000113312890132058,0.000117207101164729,0.000122294284870810,0.000117327222040188,0.000119760317903371,0.000122453750364758,0.000125623995549039,0.000110452705915810,9.76732081323682e-05,8.71238682953567e-05,7.85664528001014e-05,7.12931319423352e-05,6.44508042026831e-05,5.72739316098663e-05,4.99574540056908e-05,4.44387491003240e-05,3.86595497293785e-05,3.17529525926490e-05,2.43758225041782e-05,1.65445216648620e-05,8.66912571133985e-06,4.60089946878572e-06];
x=[0.027271727	0.028170607	0.029030873	0.028870509	0.030299687	0.031496228	0.037214804	0.0402045	0.042006046	0.045425073	0.045066543	0.044433855	0.044061764	0.043657371	0.042976798	0.04090401	0.037839473	0.036353178	0.034373859	0.034904096	0.035869532	0.035481756	0.032934135	0.026823544	0.01575832	0.009067587	0.009076604	0.009115128	0.00916492	0.009227016	0.009302398	0.009159568	0.009178833	0.009206502	0.009262611	0.008962744	0.008651363	0.008316654	0.00795962	0.007575908	0.007172649	0.006811135	0.006602081	0.006327993	0.005931433	0.00542962	0.004795852	0.004005556	0.002994765	0.001739658	0.00011314	0.000113382	0.000114783	0.000116651	0.000119396	0.000123002	0.00011471	0.000114644	0.000114304	0.000115537	0.000102167	9.15E-05	8.25E-05	7.49E-05	6.78E-05	6.07E-05	5.47E-05	5.13E-05	4.71E-05	4.14E-05	3.46E-05	2.69E-05	1.85E-05	9.75E-06	3.55E-06
];
%% run_folder

    run_folder = [
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\AR19_FWT_eta85']; %[-], folder for exporting the NASTRAN model

     
%% Safty Factor

%     SF=1.5;

%% Tip configuration 

    fold_angle  = -10;   %[deg],
    flare_angle = 25;   %[deg],
    fold_eta=0.85;
    hinge_stiffness=1e-4;
    
%% Wing configurations for starboard wing
  
    Aspect_ratio=19; % Aspect ratio = 10.172 for A321 model
    
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
    
    
%     fold_angle  = -10;   %[deg],
%     flare_angle = 25;   %[deg],
%     fold_eta=0.75;
    
    
    FWT = insertWingFold(Wingbox, 'FlareAngle', flare_angle, 'FoldAngle', fold_angle,'EtaFold',fold_eta);
    FWT.HingeStiffness = [1e14 1e14 1e14 1e14 hinge_stiffness 1e14];
    

    NumSec=Wingbox.NumBeamElem+2;
    
    kink_eta=0.27*(1/fold_eta);

    YData=Wingbox.YData;
    SparWidth=Wingbox.Chord*0.5;
    
    End_H=0.12-(1-0.27-(1-fold_eta))*0.01/(1-0.27);

    RootH=Wingbox.Chord(1)*0.15; % root thickness/chord = 0.15
    MidH=Wingbox.Chord(2)*0.12;  % middle thickness/chord = 0.12
    TipH=Wingbox.Chord(end)*End_H;% tip thickness/chord = 0.11


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
    
    Wing_mass0=1912.1;
    
    Secondary_mass0=835.2;
    
    Wing_total_mass0=Wing_mass0+Secondary_mass0;
    
    Fuselage_structure_mass=OWE - 2*Wing_total_mass0 - Horizontal_tail - Vertical_tail - Engine_mass*2 - Pylon*2 - Fuselage_shell_mass;
    
    [Fuel_mass,Fuselage_total_mass]= Lumped_masses(Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
    
    [~, Secondary_mass, ~]=Wing_masses_v1(x,Height_var,Width_var,etaRL,MTOW,Wing_area,LE_sweep,Semi_span);
    
    OWE_no_wing=OWE- 2*Wing_total_mass0-Engine_mass*2 - Pylon*2;
    
    total_attachment_mass=Pylon*2 + Engine_mass*2;

    
    %% Material properties : Al7075
    
    Yield_strength = 5.2e8;
    
    %% Folding wing tip 
    
    FWT_root_width=0.5*Wingbox.Chord(end);
    FWT_root_height=TipH;
    
    FWT_tip_width=0.5*Tip_chord;
    FWT_tip_height=Tip_chord*0.11;
    
    FWT_thickness = 0.002;
    
    %% start the loop
 
    % data record
%     result_counter=1;
%     X_result=zeros(10,75);
%     Mass_result=zeros(10,3);
       
     
    % start size loop
    %initial condition
    cond_set1=[1.1,1.1,1.1,1.1];
    cond_set2=[1.1,1.1,1.1,1.1];
    
    Numsec=25;
    counter=1;
    
    record=zeros(100,75); 
    
    while max(cond_set1)>1.01 || min(cond_set2)<0.99
        
        record(counter,:)=x(1,:);
        
        [Aircraft, FEM_full, Y_all, Load_distribution, Delta, RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=Stress_Calc_LoadCases_v3_FWT(run_folder,x,Semi_span,Root_chord,LE_sweep,TE_sweep1,TE_sweep2,BeamLoc,Fuselage_length, Wing_position, Engine_position, Horizontal_tail_position, Vertical_tail_position, Secondary_mass,Fuel_mass,Fuselage_total_mass,fold_angle,flare_angle,fold_eta,hinge_stiffness,FWT_root_width,FWT_root_height,FWT_tip_width,FWT_tip_height, FWT_thickness);

        
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
        
   
    
    [Wingbox_mass, Secondary_mass, Total_mass]=Wing_masses_v1(x,Height_var,Width_var,etaRL,MTOW,Wing_area,LE_sweep,Semi_span);
    
    [Fuel_mass,Fuselage_total_mass]= Lumped_masses(Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
    
    Wing_total_mass=Wingbox_mass+Secondary_mass;
    
    % calculate the mass of FWT

    FWT_RootArea= (FWT_root_width + FWT_root_height)*2*FWT_thickness;
    FWT_TipArea= (FWT_tip_width + FWT_tip_height)*2*FWT_thickness;
    
    FWT_mass = (FWT_RootArea + FWT_TipArea)*FWT.Span*2800/2;
    
    Aircraft_total_mass=Fuselage_total_mass+Fuselage_shell_mass+Wing_total_mass*2+Horizontal_tail+ Vertical_tail + Engine_mass*2 + Pylon*2 + Fuel_mass*2;
    
    
    %% vertical wing displacement 
    
     NastranMethods1 = awi.methods.Nastran;
     NastranMethods1.AnalysisModel=FEM_full;
     
    
    WingNodes=NastranMethods1.WingNode;
    FWTNodes=NastranMethods1.FWTNode;
    WingNodes_All=[WingNodes(1:end-1),FWTNodes];
    
    
    Grids = h5read(strcat(run_folder,'\a321_36000ft_1g.h5'),'/NASTRAN/INPUT/NODE/GRID');
    
    Index=ismember([Grids.ID],[WingNodes_All(1:end).GID]);
    
    Displacement=h5read(strcat(run_folder,'\a321_36000ft_1g.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');
    
    Wing_displacement=Displacement.Z(Index);
    

    figure 
    plot(Y_all,Wing_displacement,'bs-')

    Wing_displacement=Wing_displacement';


% figure 
% 
% for i= 24:34
%     plot(Y_all(1:25),record(i,1:25))
%     hold on 
%     set(gcf,'Color','w')
%     xlabel('Span Distance (m)', 'Interpreter','latex','FontSize',12)
%     ylabel('Spar thickness (mm)', 'Interpreter','latex','FontSize',12)
%     
% end
% 
% 
% figure 
% 
% for i= 24:34
%     plot(Y_all(1:25),record(i,26:50))
%     hold on
%     set(gcf,'Color','w')
%     xlabel('Span Distance (m)', 'Interpreter','latex','FontSize',12)
%     ylabel('Skin thickness (mm)', 'Interpreter','latex','FontSize',12)
%     
%     
% end
    
%     %% data logger
%     X_result(result_counter,:)=x;
%     
%     Mass_result(result_counter,:)=[Wingbox_mass,Secondary_mass,Aircraft_total_mass];
%     
%     result_counter=result_counter+1;
%     
%     
% 
% end    

        

    
    % obtain lumped masses: Secondar_mass and fuel mass are the masses for
    % each wingbox. Fuselage_total_mass includs the structure mass and
    % payload
    
    % the sizing is carroied out using max. payload with its corresponding
    % max fuel mass 
    
%     figure 
%     
%     plot(Y_all(1:end-1),Delta.Moment_P2_Max.LC1,'b.')
    
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
    
    
    
    


