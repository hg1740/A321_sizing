

%% initial guesses

%A320
% x=[0.0160369401546491,0.0161179748570384,0.0165944334426013,0.0175266294784001,0.0187239043216224,0.0207096617121334,0.0230583809562627,0.0320089957111805,0.0320392423264178,0.0319427876288856,0.0315322529714459,0.0322935890059832,0.0323900672122562,0.0323115468409948,0.0321955616764684,0.0313429255556067,0.0311508861821847,0.0304592682549718,0.0297595013116213,0.0272724891942158,0.0271122094889916,0.0264611574714236,0.0230559955692698,0.0187199645209936,0.0107524736198303,0.00647614132210356,0.00657526982363507,0.00669537365768939,0.00683776660015478,0.00701687360523511,0.00721883252923472,0.00748254587435518,0.00755456273799370,0.00737352569279649,0.00716564274351893,0.00696475875843593,0.00674100541947114,0.00649425980510992,0.00623148973761537,0.00595214725710100,0.00564842457162947,0.00541126353104570,0.00522201059650361,0.00497676988997949,0.00469043302797361,0.00427877225656423,0.00376751794116612,0.00307441811486333,0.00231764849271284,0.000789440488034985,5.19302425632857e-05,5.34268357423746e-05,5.52318694509354e-05,5.74180048238046e-05,6.02580769479868e-05,6.35124694224661e-05,6.77949943208388e-05,6.87713796834609e-05,6.55154629087215e-05,6.17954684314261e-05,5.82619158760843e-05,5.45303477290052e-05,5.06453591835854e-05,4.65755305830927e-05,4.24490132198509e-05,3.80255976131977e-05,3.49566602319797e-05,3.25217162959330e-05,2.94290956549030e-05,2.60051023480533e-05,2.16034838172783e-05,1.67146290536804e-05,1.13000939281217e-05,6.55690899275516e-06,7.23256842198466e-07];

%A321
% x=[0.0191135905034748,0.0193953734978214,0.0201137604804691,0.0213823646563527,0.0229999860459657,0.0255251238127354,0.0285624381189545,0.0372113474624173,0.0373207109980618,0.0372799899057456,0.0368848456350822,0.0378039760454905,0.0379687298578435,0.0379303106594230,0.0378464169968381,0.0368959067148779,0.0365928771944430,0.0358740997358953,0.0351157988905984,0.0325204105754142,0.0321366355707771,0.0317223020101411,0.0280889429523174,0.0226448795839794,0.0118821391075936,0.00707914701666205,0.00718361444866246,0.00730636454511555,0.00745189886461910,0.00763438438081004,0.00783787156936675,0.00810313558972640,0.00823100835918603,0.00803824902557083,0.00781619282888881,0.00760074969073965,0.00736034765742680,0.00709415722939268,0.00680991625431645,0.00650766890180177,0.00618441788926231,0.00597655499191214,0.00574729769090399,0.00545697336175644,0.00516124669194006,0.00473319707846889,0.00418339089852517,0.00343190389276599,0.00256919715920123,0.000896999342517673,6.20362573091953e-05,6.37511021940462e-05,6.57911120554895e-05,6.82358984870429e-05,7.13734187868434e-05,7.49159107475837e-05,7.98695370113729e-05,8.20549947219400e-05,7.77731126376305e-05,7.33381910989789e-05,6.92092572133319e-05,6.48380207016889e-05,6.02644951660675e-05,5.54647083143445e-05,5.05953480091047e-05,4.54572672185587e-05,4.24462089148306e-05,3.91789917569685e-05,3.51963688050882e-05,3.13298871582112e-05,2.63163732081641e-05,2.05279963925629e-05,1.40184642076611e-05,8.01533883220779e-06,9.15451870145967e-07];

% x=[0.0263919314434434,0.0270287037846719,0.0283640513655332,0.0290476918785440,0.0305000319425819,0.0311454414799324,0.0396804988555224,0.0425444017731640,0.0441525651871870,0.0482605196982675,0.0464068437203061,0.0457224486348197,0.0456279338158660,0.0448628529788929,0.0435569212840316,0.0417731169704868,0.0396815956770828,0.0380480440213013,0.0358316633566355,0.0339687923155953,0.0337268684564342,0.0323902190802610,0.0291451850232813,0.0234232299918240,0.0191524074792969,0.00893809271271453,0.00898021206205229,0.00902935314915270,0.00908668154345703,0.00917673153841885,0.00928838436391813,0.00921342865699722,0.00927695325733132,0.00935527519615850,0.00946264343682268,0.00916930010549819,0.00885042818227504,0.00850732408723034,0.00814899922819781,0.00777157655127306,0.00739408745863920,0.00696700999391705,0.00650587525589280,0.00613500379948089,0.00572350270401927,0.00518759260675683,0.00454772718759639,0.00375606295531985,0.00275656826201626,0.00207551938332536,0.000107769528835642,0.000109277971748372,0.000111098275638144,0.000113312890132058,0.000117207101164729,0.000122294284870810,0.000117327222040188,0.000119760317903371,0.000122453750364758,0.000125623995549039,0.000110452705915810,9.76732081323682e-05,8.71238682953567e-05,7.85664528001014e-05,7.12931319423352e-05,6.44508042026831e-05,5.72739316098663e-05,4.99574540056908e-05,4.44387491003240e-05,3.86595497293785e-05,3.17529525926490e-05,2.43758225041782e-05,1.65445216648620e-05,8.66912571133985e-06,4.60089946878572e-06];
x=[0.020125557	0.019401577	0.021958091	0.022539262	0.025767917	0.025652551	0.032682526	0.035610407	0.037135288	0.036060835	0.036340768	0.036871667	0.035677056	0.035171832	0.03289602	0.032989368	0.03199809	0.029865572	0.029124862	0.027113727	0.025639198	0.025349052	0.022339044	0.017219893	0.015604291	0.007430958	0.007570034	0.007675861	0.007863728	0.008039501	0.008326083	0.008388571	0.008592074	0.008327159	0.008068052	0.007797718	0.007519421	0.007227152	0.006914659	0.006600855	0.006252395	0.005858025	0.005527251	0.005183017	0.004799195	0.004359102	0.003862057	0.003188274	0.002330943	0.001702933	6.82E-05	7.06E-05	7.25E-05	7.57E-05	7.90E-05	8.50E-05	8.67E-05	9.15E-05	8.42E-05	7.83E-05	7.29E-05	6.79E-05	6.27E-05	5.72E-05	5.21E-05	4.69E-05	4.13E-05	3.65E-05	3.19E-05	2.77E-05	2.26E-05	1.70E-05	1.16E-05	6.12E-06	3.21E-06
];
%% run_folder

    run_folder = [
        'C:\Git\A321_sizing\hg_codes\results\Noextension1']; %[-], folder for exporting the NASTRAN model

     
%% Safty Factor

%     SF=1.5;


    
%% Wing configurations for starboard wing
  
    Wing_area=126;
    
    LE_sweep=27;            % deg
    
    TE_sweep1=0; % deg
    TE_sweep2=21; % deg
    
    Wing_span=35.8;
    
    Semi_span=(Wing_span-4)/2;
    Root_chord=6;
    BeamLoc = 0.4; 
    
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
    
    while max(cond_set1)>1.04 || min(cond_set2)<0.96
        
        record(counter,:)=x(1,:);
        
        [Aircraft, FEM_full, Y_Nodes, Load_distribution, Delta, RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=Stress_Calc_LoadCases_v3_Noextension(run_folder,x,Semi_span,Root_chord,LE_sweep,TE_sweep1,TE_sweep2,BeamLoc,Fuselage_length, Wing_position, Engine_position, Horizontal_tail_position, Vertical_tail_position, Secondary_mass,Fuel_mass,Fuselage_total_mass);

        
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
    
    
    Aircraft_total_mass=Fuselage_total_mass+Fuselage_shell_mass+Wing_total_mass*2+Horizontal_tail+ Vertical_tail + Engine_mass*2 + Pylon*2 + Fuel_mass*2;
    
    
    %% Drag 
    
    % create an object array 
    All_FEM=flatlist(FEM_full);
    
    Wing_conn=All_FEM(ismember({All_FEM.Name}, 'Connector_Right'));
    Wing_inboard = All_FEM(ismember({All_FEM.Name}, 'A320Wing_right'));
    
    Wing_conn_ID=Wing_conn.AeroSets.E;
    Wing_inboard_ID=Wing_inboard.AeroSets(1).E;
    Wing_outboard_ID=Wing_inboard.AeroSets(2).E;
    
%     total_ID=[Wing_conn_ID,Wing_inboard_ID,FWT_ID];
    
    % aeroforces at each segment
    res_aeroF = mni.result.f06(strcat('C:\Git\A321_sizing\hg_codes\results\Noextension1','/A321_36000ft_1g.f06')).read_aeroF;
    
    Index_conn=ismember([res_aeroF.PanelID(1,:)],Wing_conn_ID);
    Index_wing_inboard=ismember([res_aeroF.PanelID(1,:)],Wing_inboard_ID);
    Index_wing_outboard=ismember([res_aeroF.PanelID(1,:)],Wing_outboard_ID);
    
    Lift_conn=res_aeroF.aeroFz(Index_conn);
    Lift_wing_inboard=res_aeroF.aeroFz(Index_wing_inboard);
    Lift_wing_outboard=res_aeroF.aeroFz(Index_wing_outboard);
    
    
    % calculate panel width for each segment: conn + wing_bef_kink +
    % wing_aft_kink 
    
    num_panel_conn=numel(Wing_conn_ID)/10;
    num_panel_inboard=numel(Wing_inboard_ID)/10;
    num_panel_outboard=numel(Wing_outboard_ID)/10;
    
    R_chord=Wingbox.Chord(1);
    semi_span=Wingbox.Span;
    conn=2;
    kink_pos=Wingbox.YData(2);
    
    width_conn=2/num_panel_conn;
    
    width_wing1=kink_pos/num_panel_inboard;
    width_wing2=(Wingbox.YData(3)-Wingbox.YData(2))/num_panel_outboard;

    panel_width =[width_conn*ones(1,num_panel_conn),width_wing1*ones(1,num_panel_inboard),width_wing2*ones(1,num_panel_outboard)];

    
    % Y-Data
    % find corresponding Y positions
    Y_conn=0.5*width_conn:width_conn*0.999:2;
    
    Y_wing1=2+0.5*width_wing1:width_wing1*0.999:2+Wingbox.YData(2);
    
    Y_wing2=2+Wingbox.YData(2)+ 0.5*width_wing2:width_wing2*0.999:2+Wingbox.YData(3);
    
    Y=[Y_conn,Y_wing1,Y_wing2];
    
    % normalise lift by panel width to obtain lift per unit span
    lift_wing=[Lift_conn';Lift_wing_inboard';Lift_wing_outboard'];
    
    lift_wing_matrix=reshape(lift_wing, 10, numel(lift_wing)/10);
    
    % abs. value of the lift
    lift_wing_abs=sum(lift_wing_matrix);
    
    Cum_num=cumsum([num_panel_conn,num_panel_inboard,num_panel_outboard]);
    
    lift_wing_matrix(:,1:Cum_num(1))=lift_wing_matrix(:,1:Cum_num(1))/width_conn;
    lift_wing_matrix(:,Cum_num(1)+1:Cum_num(2))=lift_wing_matrix(:,Cum_num(1)+1:Cum_num(2))/width_wing1;
    lift_wing_matrix(:,Cum_num(2)+1:Cum_num(3))=lift_wing_matrix(:,Cum_num(2)+1:Cum_num(3))/width_wing2;
    
    
    % lift per unit span
    lift_wing_var=sum(lift_wing_matrix);
    
    % whole wing span lift distribution
    Y_left=sort(-Y);
    lift_wing_var_left=flip(lift_wing_var);
    lift_wing_abs_left=flip(lift_wing_abs);
    panel_width_left=flip(panel_width);
    
    Y_all=[Y_left,Y]';
    Lift_all=[lift_wing_var_left,lift_wing_var];
    Lift_all_abs=[lift_wing_abs_left,lift_wing_abs];
    
    panel_width_all=[panel_width_left,panel_width];
    
    
    % calculate vorticity
    % Air conditions
    FlightPoint=awi.model.FlightPoint;
    FlightPoint.Mach=0.78;
    FlightPoint.AcVelocity=0.78*340;
    FlightPoint.Altitude = 36000;
    getFlightPointData(FlightPoint,'ISA');
    DynPressure = FlightPoint.DynPressure;
    
    Gamma_all=Lift_all/(FlightPoint.AirDensity*FlightPoint.AcVelocity);
    
    % calculate the derivative
    dGdy= gradient(Gamma_all(:)) ./ gradient(Y_all(:));
    
    % downwash
    wj=zeros(1,length(Y_all));
    alphai=zeros(1,length(Y_all));
    
    
    for i=1:length(Y_all)
        
        w=-(1/(4*pi))*dGdy* panel_width_all(i)./(Y_all(i)-Y_all);
        
        w=w( ~any( isnan( w ) | isinf( w ), 2 ),: );
        
        wj(i)=sum(w);
        
        alphai(i)=wj(i)/FlightPoint.AcVelocity;
        
    end
    
    Dragi_var=Lift_all_abs.*sin(alphai);
    Drag_force=sum([Dragi_var]);
    Total_area=(Wingbox.SurfaceArea)*2+Wingbox.Chord(1)*4;
    Cdi=Drag_force/(DynPressure*Total_area);
    
    
    figure
    Y=Y';
    lift_wing_var=lift_wing_var';
    plot(Y,lift_wing_var,'bs-','LineWidth',1)
    
    hold on 
    
     plot(Y1,lift_wing_var1,'rs-','LineWidth',1)
    
    set(gcf,'Color','w')
    
    xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
    
    ylabel('Lift per unit span (N/m)','Interpreter','latex','FontSize',12)
    legend('No extension','FWT extension','Interpreter','latex','FontSize',12)
    
    figure
    Y=Y';
    lift_wing_var=lift_wing_var';
    plot(Y_all,alphai,'rs-','LineWidth',1)
    set(gcf,'Color','w')
    
    xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
    
    ylabel('$\alpha_i$','Interpreter','latex','FontSize',12)
    
    
    % spar thickness
    
    figure
    
    plot(Y2(1:25),x2(1:25),'bs','MarkerFaceColor','b')
    hold on
    plot(Y1(1:25),x1(1:25),'rs','MarkerFaceColor','r')
    
    set(gcf,'Color','w')
    
    xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
    
    ylabel('Spar thickness (m)','Interpreter','latex','FontSize',12)
    
    legend('No extension','FWT extension','Interpreter','latex','FontSize',12)
    
    % skin thickness
    figure
    
    plot(Y2(1:25),x2(26:50),'bs','MarkerFaceColor','b')
    hold on
    plot(Y1(1:25),x1(26:50),'rs','MarkerFaceColor','r')
    
    set(gcf,'Color','w')
    
    xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
    
    ylabel('Skin thickness (m)','Interpreter','latex','FontSize',12)
    
    legend('No extension','FWT extension','Interpreter','latex','FontSize',12)
    
    % stringer
    figure
    
    plot(Y2(1:25),x2(51:75),'bs','MarkerFaceColor','b')
    hold on
    plot(Y1(1:25),x1(51:75),'rs','MarkerFaceColor','r')
    
    set(gcf,'Color','w')
    
    xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
    
    ylabel('Stringer area (m$^{2}$)','Interpreter','latex','FontSize',12)
    
    legend('No extension','FWT extension','Interpreter','latex','FontSize',12)
    
    
    
    
    
    
    

    
    
    
    
    %% vertical wing displacement 
    
%      NastranMethods1 = awi.methods.Nastran;
%      NastranMethods1.AnalysisModel=FEM_full;
%      
%     
%     WingNodes=NastranMethods1.WingNode;
%     FWTNodes=NastranMethods1.FWTNode;
%     WingNodes_All=[WingNodes(1:end-1),FWTNodes];
%     
%     
%     Grids = h5read(strcat(run_folder,'\a321_36000ft_1g.h5'),'/NASTRAN/INPUT/NODE/GRID');
%     
%     Index=ismember([Grids.ID],[WingNodes_All(1:end).GID]);
%     
%     Displacement=h5read(strcat(run_folder,'\a321_36000ft_1g.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');
%     
%     Wing_displacement=Displacement.Z(Index);
%     
% 
%     figure 
%     plot(Y_all,Wing_displacement,'bs-')
% 
%     Wing_displacement=Wing_displacement';


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
    
    
    
    


