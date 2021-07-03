
%% initial guesses

%A320
% x=[0.0160369401546491,0.0161179748570384,0.0165944334426013,0.0175266294784001,0.0187239043216224,0.0207096617121334,0.0230583809562627,0.0320089957111805,0.0320392423264178,0.0319427876288856,0.0315322529714459,0.0322935890059832,0.0323900672122562,0.0323115468409948,0.0321955616764684,0.0313429255556067,0.0311508861821847,0.0304592682549718,0.0297595013116213,0.0272724891942158,0.0271122094889916,0.0264611574714236,0.0230559955692698,0.0187199645209936,0.0107524736198303,0.00647614132210356,0.00657526982363507,0.00669537365768939,0.00683776660015478,0.00701687360523511,0.00721883252923472,0.00748254587435518,0.00755456273799370,0.00737352569279649,0.00716564274351893,0.00696475875843593,0.00674100541947114,0.00649425980510992,0.00623148973761537,0.00595214725710100,0.00564842457162947,0.00541126353104570,0.00522201059650361,0.00497676988997949,0.00469043302797361,0.00427877225656423,0.00376751794116612,0.00307441811486333,0.00231764849271284,0.000789440488034985,5.19302425632857e-05,5.34268357423746e-05,5.52318694509354e-05,5.74180048238046e-05,6.02580769479868e-05,6.35124694224661e-05,6.77949943208388e-05,6.87713796834609e-05,6.55154629087215e-05,6.17954684314261e-05,5.82619158760843e-05,5.45303477290052e-05,5.06453591835854e-05,4.65755305830927e-05,4.24490132198509e-05,3.80255976131977e-05,3.49566602319797e-05,3.25217162959330e-05,2.94290956549030e-05,2.60051023480533e-05,2.16034838172783e-05,1.67146290536804e-05,1.13000939281217e-05,6.55690899275516e-06,7.23256842198466e-07];

%A321
% x=[0.0191135905034748,0.0193953734978214,0.0201137604804691,0.0213823646563527,0.0229999860459657,0.0255251238127354,0.0285624381189545,0.0372113474624173,0.0373207109980618,0.0372799899057456,0.0368848456350822,0.0378039760454905,0.0379687298578435,0.0379303106594230,0.0378464169968381,0.0368959067148779,0.0365928771944430,0.0358740997358953,0.0351157988905984,0.0325204105754142,0.0321366355707771,0.0317223020101411,0.0280889429523174,0.0226448795839794,0.0118821391075936,0.00707914701666205,0.00718361444866246,0.00730636454511555,0.00745189886461910,0.00763438438081004,0.00783787156936675,0.00810313558972640,0.00823100835918603,0.00803824902557083,0.00781619282888881,0.00760074969073965,0.00736034765742680,0.00709415722939268,0.00680991625431645,0.00650766890180177,0.00618441788926231,0.00597655499191214,0.00574729769090399,0.00545697336175644,0.00516124669194006,0.00473319707846889,0.00418339089852517,0.00343190389276599,0.00256919715920123,0.000896999342517673,6.20362573091953e-05,6.37511021940462e-05,6.57911120554895e-05,6.82358984870429e-05,7.13734187868434e-05,7.49159107475837e-05,7.98695370113729e-05,8.20549947219400e-05,7.77731126376305e-05,7.33381910989789e-05,6.92092572133319e-05,6.48380207016889e-05,6.02644951660675e-05,5.54647083143445e-05,5.05953480091047e-05,4.54572672185587e-05,4.24462089148306e-05,3.91789917569685e-05,3.51963688050882e-05,3.13298871582112e-05,2.63163732081641e-05,2.05279963925629e-05,1.40184642076611e-05,8.01533883220779e-06,9.15451870145967e-07];

x=[0.0272031060000000,0.0284343580000000,0.0295254410000000,0.0303307180000000,0.0313258800000000,0.0398912860000000,0.0445268050000000,0.0466773110000000,0.0458190000000000,0.0455495620000000,0.0436103050000000,0.0430748590000000,0.0411122050000000,0.0409634860000000,0.0403045640000000,0.0374371890000000,0.0353579660000000,0.0314095960000000,0.0331796210000000,0.0379988910000000,0.0410566640000000,0.0452771000000000,0.0457390880000000,0.0398047130000000,0.0106772500000000,0.00927907500000000,0.00936425700000000,0.00946086900000000,0.00961691300000000,0.00979752600000000,0.00982986500000000,0.00997293200000000,0.0102379670000000,0.00993122700000000,0.00960289300000000,0.00926602000000000,0.00890772000000000,0.00855306000000000,0.00815705500000000,0.00773865200000000,0.00772948200000000,0.00776624800000000,0.00774718100000000,0.00759222600000000,0.00729932500000000,0.00688257200000000,0.00630456700000000,0.00538514500000000,0.00410838700000000,0.00168449400000000,0.000125803000000000,0.000130510000000000,0.000136150000000000,0.000146897000000000,0.000161351000000000,0.000161664000000000,0.000172723000000000,0.000198329000000000,0.000159929000000000,0.000132782000000000,0.000113342000000000,9.84000000000000e-05,8.74000000000000e-05,7.79000000000000e-05,7.00000000000000e-05,6.96000000000000e-05,6.99000000000000e-05,6.93000000000000e-05,6.63000000000000e-05,6.11000000000000e-05,5.42000000000000e-05,4.55000000000000e-05,3.34000000000000e-05,1.97000000000000e-05,2.69000000000000e-06];
%% run_folder

    run_folder = [
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\AR19_ref_gust']; %[-], folder for exporting the NASTRAN model

     
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
    
    % reference ac A321 wing box weight, subtracted from OWE tp generate a
    % constant mass value for the sizing of new configurations.
    Wingbox_mass0=1912.1;
    
    Secondary_mass0=835.2;
    
    Wing_total_mass0=Wingbox_mass0+Secondary_mass0;
    
    Fuselage_structure_mass=OWE - 2*Wing_total_mass0 - Horizontal_tail - Vertical_tail -Engine_mass*2 - Pylon*2-Fuselage_shell_mass;
    
    [Fuel_mass,Fuselage_total_mass]= Lumped_masses(Fuel_fraction,Fuel_capacity,Payload_max,Fuselage_structure_mass);
    
    [~, Secondary_mass, ~]=Wing_masses_v1(x,Height_var,Width_var,etaRL,MTOW,Wing_area,LE_sweep,Semi_span);


    %% Run analysis
    
[Aircraft,FEM_full,Y_Data, Load_distribution, Delta, RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=Stress_Calc_LoadCases_v3(run_folder,x,Semi_span,Root_chord,LE_sweep,TE_sweep1,TE_sweep2, BeamLoc,Fuselage_length, Wing_position, Engine_position, Horizontal_tail_position, Vertical_tail_position, Secondary_mass,Fuel_mass,Fuselage_total_mass);
        

 %%   Reasult plot
 
    Y_Data=Y_all;
    
    figure 
    plot(Y_Data-2,Load_distribution.Shear_P2.LC1,'rs','MarkerFaceColor','r')
    hold on
    plot(Y_Data-2,Load_distribution.Shear_P2.LC2,'bs','MarkerFaceColor','b')
    hold on
    plot(Y_Data-2,Load_distribution.Shear_P2.LC3,'ks','MarkerFaceColor','k')
    hold on
    plot(Y_Data-2,Load_distribution.Shear_P2.LC4,'r-','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Shear_P2.LC5,'b--','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Shear_P2.LC6,'k-.','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Shear_P2.LC7,'ro','MarkerFaceColor','r')
    hold on
    plot(Y_Data-2,Load_distribution.Shear_P2.LC8,'bo','MarkerFaceColor','b')
    hold on
    plot(Y_Data-2,Load_distribution.Shear_P2.LC9,'ko','MarkerFaceColor','k')
%     hold on
%     plot(Y-2,S_P2,'r--','LineWidth',1.2)
    
%     legend('2.5 g','2.5 g','-g ','g + 1MC','g + 1MC', 'g + 1MC','g - 1MC','g - 1MC', 'g - 1MC','Interpreter','latex','FontSize',10)
    legend('2.5 g (36 kft)','2.5 g (3 kft)','-g (3 kft)','g + 1MC (36 kft)','g + 1MC (20 kft)', 'g + 1MC (3 kft)','g - 1MC (36 kft)','g - 1MC (20 kft)', 'g - 1MC (3 kft)','Interpreter','latex','FontSize',10)
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('Vertical shear force (N)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')
    
    
    
    figure 
    plot(Y_Data-2,Load_distribution.Moment_P2.LC1,'rs','MarkerFaceColor','r')
    hold on
    plot(Y_Data-2,Load_distribution.Moment_P2.LC2,'bs','MarkerFaceColor','b')
    hold on
    plot(Y_Data-2,Load_distribution.Moment_P2.LC3,'ks','MarkerFaceColor','k')
    hold on
    plot(Y_Data-2,Load_distribution.Moment_P2.LC4,'r-','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Moment_P2.LC5,'b--','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Moment_P2.LC6,'k-.','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Moment_P2.LC7,'ro','MarkerFaceColor','r')
    hold on
    plot(Y_Data-2,Load_distribution.Moment_P2.LC8,'bo','MarkerFaceColor','b')
    hold on
    plot(Y_Data-2,Load_distribution.Moment_P2.LC9,'ko','MarkerFaceColor','k')
%     hold on
%     plot(Y-2,S_P2,'r--','LineWidth',1.2)
    
%     legend('2.5 g','2.5 g','-g ','g + 1MC','g + 1MC', 'g + 1MC','g - 1MC','g - 1MC', 'g - 1MC','Interpreter','latex','FontSize',10)
    legend('2.5 g (36 kft)','2.5 g (3 kft)','-g (3 kft)','g + 1MC (36 kft)','g + 1MC (20 kft)', 'g + 1MC (3 kft)','g - 1MC (36 kft)','g - 1MC (20 kft)', 'g - 1MC (3 kft)','Interpreter','latex','FontSize',10)
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('Bending Moment (Nm)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')
    

    figure 
    plot(Y_Data-2,Load_distribution.Torque.LC1,'r-s','MarkerFaceColor','r')
    hold on
    plot(Y_Data-2,Load_distribution.Torque.LC2,'b-s','MarkerFaceColor','b')
    hold on
    plot(Y_Data-2,Load_distribution.Torque.LC3,'k-s','MarkerFaceColor','k')
    hold on
    plot(Y_Data-2,Load_distribution.Torque.LC4,'r-','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Torque.LC5,'b--','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Torque.LC6,'k-.','LineWidth',0.8)
    hold on
    plot(Y_Data-2,Load_distribution.Torque.LC7,'--ro','MarkerFaceColor','r')
    hold on
    plot(Y_Data-2,Load_distribution.Torque.LC8,'--bo','MarkerFaceColor','b')
    hold on
    plot(Y_Data-2,Load_distribution.Torque.LC9,'--ko','MarkerFaceColor','k')
%     hold on
%     plot(Y-2,S_P2,'r--','LineWidth',1.2)
    
%     legend('2.5 g','2.5 g','-g ','g + 1MC','g + 1MC', 'g + 1MC','g - 1MC','g - 1MC', 'g - 1MC','Interpreter','latex','FontSize',10)
    legend('2.5 g (36 kft)','2.5 g (3 kft)','-g (3 kft)','g + 1MC (36 kft)','g + 1MC (20 kft)', 'g + 1MC (3 kft)','g - 1MC (36 kft)','g - 1MC (20 kft)', 'g - 1MC (3 kft)','Interpreter','latex','FontSize',10)
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('Torque (Nm)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')

    
    figure 
    plot(Y_Data(1:end-1)-2,Delta.Moment_P2_Max.LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y_Data(1:end-1)-2,Delta.Moment_P2_Max.LC2,'b.','MarkerSize',10)
    hold on
    plot(Y_Data(1:end-1)-2,Delta.Moment_P2_Max.LC3,'r.','MarkerSize',10)
    
    hold on 
    
    plot(Y_Data(1:end-1)-2,Delta.Moment_P2_Min.LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y_Data(1:end-1)-2,Delta.Moment_P2_Min.LC2,'b.','MarkerSize',10)
    hold on
    plot(Y_Data(1:end-1)-2,Delta.Moment_P2_Min.LC3,'r.','MarkerSize',10)
 
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('$\Delta$ Bending Moment (Nm)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')
    legend('36 kft +1MC','20 kft +1MC','3 ~ kft +1MC','Interpreter','latex','FontSize',12)
    

    figure 
    plot(Y_Data(1:end-1)-2,Delta.Torque_Max.LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y_Data(1:end-1)-2,Delta.Torque_Max.LC2,'b.','MarkerSize',10)
    hold on
    plot(Y_Data(1:end-1)-2,Delta.Torque_Max.LC3,'r.','MarkerSize',10)
 
    hold on
    
    plot(Y_Data(1:end-1)-2,Delta.Torque_Min.LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y_Data(1:end-1)-2,Delta.Torque_Min.LC2,'b.','MarkerSize',10)
    hold on
    plot(Y_Data(1:end-1)-2,Delta.Torque_Min.LC3,'r.','MarkerSize',10)
    
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('$\Delta$ Torque (Nm)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')
    legend('36 kft +1MC','20 kft +1MC','3 ~ kft +1MC','Interpreter','latex','FontSize',12)


    figure 
    plot(Y_Data(1:end-1)-2,Delta.Shear_Max.LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y_Data(1:end-1)-2,Delta.Shear_Max.LC2,'b.','MarkerSize',10)
    hold on
    plot(Y_Data(1:end-1)-2,Delta.Shear_Max.LC3,'r.','MarkerSize',10)
 
    hold on

    plot(Y_Data(1:end-1)-2,Delta.Shear_Min.LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y_Data(1:end-1)-2,Delta.Shear_Min.LC2,'b.','MarkerSize',10)
    hold on
    plot(Y_Data(1:end-1)-2,Delta.Shear_Min.LC3,'r.','MarkerSize',10)
    
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('$\Delta$ Vertical shear force (Nm)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')
    legend('36 kft +1MC','20 kft +1MC','3 ~ kft +1MC','Interpreter','latex','FontSize',12)
    
    
    

Gust_data = h5read(strcat(run_folder,'\gust_analysis_g_36000ft_pos_1MC.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');

NastranMethods1 = awi.methods.Nastran;
NastranMethods1.AnalysisModel = FEM_full;

WingNodes=NastranMethods1.WingNode;

FWTNodes=NastranMethods1.FWTNode;

WingNodes=[WingNodes(1:end-1),FWTNodes];

Steps=201;
Z_disp_result=ones(Steps*7,numel(WingNodes));


%group 1: wing Z displacement ------------------------------------------------------------------

% ID=WingNodes(1:25).GID;

% Z_disp_wing=ones
for i=1:numel(WingNodes)
    
    [index,~]=find(Gust_data.ID==WingNodes(i).GID); 
    
    Z_disp_result(:,i)=Gust_data.Z([index]);
    
end

gust1= Z_disp_result(1:Steps,:);
gust2= Z_disp_result(Steps+1:Steps*2,:);
gust3= Z_disp_result(Steps*2+1:Steps*3,:);
gust4= Z_disp_result(Steps*3+1:Steps*4,:);
gust5= Z_disp_result(Steps*4+1:Steps*5,:);
gust6= Z_disp_result(Steps*5+1:Steps*6,:);
gust7= Z_disp_result(Steps*6+1:Steps*7,:);


% curve = animatedline;

% 
% for i=1:201
%     clearpoints(curve)
%     x=1:1:25;
%     addpoints(curve,x,gust1(i,:))
%     
%    drawnow
%     
%     
% end
% close all
x=[1:1:35]*22/35;

y=gust6;

figure 
ph=plot(x,y(100,:),'b-s');
set(gca,'Xlim',[0 25],'Ylim',[-2 2])
set(gcf,'Color','w')
xlabel('Wing span','Interpreter','latex')
ylabel('Vertical displacement','Interpreter','latex')

for i=1:Steps
    ph.XData=x;
    ph.YData=y(i,:);
    
    drawnow
    pause(0.05)
    
end
    
%group 1: fuselage Z displacement ------------------------------------------------------------------

FuselageNodes=NastranMethods1.RefNode;
num=numel(FuselageNodes);
Z_disp_fuselage=ones(Steps*10,numel(FuselageNodes));


% Z_disp_wing=ones
for i=1:numel(FuselageNodes)
    
    [index,~]=find(Gust_data.ID==FuselageNodes(i).GID); 
    
    Z_disp_fuselage(:,i)=Gust_data.Z([index]);
    
end

gust1_Z_fuselage = Z_disp_fuselage(1:Steps,:);
gust2_Z_fuselage= Z_disp_fuselage(Steps+1:Steps*2,:);
gust3_Z_fuselage= Z_disp_fuselage(Steps*2+1:Steps*3,:);
gust4_Z_fuselage= Z_disp_fuselage(Steps*3+1:Steps*4,:);
gust5_Z_fuselage= Z_disp_fuselage(Steps*4+1:Steps*5,:);
gust6_Z_fuselage= Z_disp_fuselage(Steps*5+1:Steps*6,:);
gust7_Z_fuselage= Z_disp_fuselage(Steps*6+1:Steps*7,:);
gust8_Z_fuselage= Z_disp_fuselage(Steps*7+1:Steps*8,:);
gust9_Z_fuselage= Z_disp_fuselage(Steps*8+1:Steps*9,:);
gust10_Z_fuselage = Z_disp_fuselage(Steps*9+1:Steps*10,:);

x=1:1:num;
y=gust6_Z_fuselage;

figure 
ph=plot(x,y(100,:),'bs-');
set(gca,'Xlim',[0 75],'Ylim',[-2 2])
set(gcf,'Color','w')
xlabel('Fuselage','Interpreter','latex')
ylabel('Vertical displacement','Interpreter','latex')

for i=1:Steps
    ph.XData=x;
    ph.YData=y(i,:);
    
    drawnow
    pause(0.05)
    
end



% group 3: wing tip twist  ------------------------------------------------------------------

WingNodes=NastranMethods1.WingNode;
% ID=WingNodes(1:25).GID;
Wing_twist=ones(Steps*10,numel(WingNodes));

% Z_disp_wing=ones
for i=1:numel(WingNodes)
    
    [index,~]=find(Gust_data.ID==WingNodes(i).GID); 
    
    Wing_twist(:,i)=Gust_data.RX([index]);
    
end

gust1t= Wing_twist(1:Steps,:);
gust2t= Wing_twist(Steps+1:Steps*2,:);
gust3t= Wing_twist(Steps*2+1:Steps*3,:);
gust4t= Wing_twist(Steps*3+1:Steps*4,:);
gust5t= Wing_twist(Steps*4+1:Steps*5,:);
gust6t= Wing_twist(Steps*5+1:Steps*6,:);
gust7t= Wing_twist(Steps*6+1:Steps*7,:);
gust8t= Wing_twist(Steps*7+1:Steps*8,:);
gust9t= Wing_twist(Steps*8+1:Steps*9,:);
gust10t = Wing_twist(Steps*9+1:Steps*10,:);



x=1:1:25;
y=gust6t;

figure 
ph=plot(x,y(100,:),'b-s');
set(gca,'Xlim',[0 25],'Ylim',[-0.01 0.01])
set(gcf,'Color','w')
xlabel('Wing span','Interpreter','latex')
ylabel('Twist','Interpreter','latex')

for i=1:Steps
    ph.XData=x;
    ph.YData=y(i,:);
    
    drawnow
    pause(0.05)
    
end 



%% stress plot


    figure 
    plot(Y(1:24)-2,Max_Moment_LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y(1:24)-2,Max_Moment_LC2,'b.','MarkerSize',10)
    hold on
    plot(Y(1:24)-2,Max_Moment_LC3,'r.','MarkerSize',10)
 
%     hold on
%     
%     plot(Y(1:24)-2,Max_Moment_LC4,'k-.','LineWidth',0.8)
%     hold on
%     plot(Y(1:24)-2,Max_Moment_LC5,'b-.','LineWidth',0.8)
%     hold on
%     plot(Y(1:24)-2,Max_Moment_LC6,'r-.','LineWidth',0.8)
    
    hold on 
    
    plot(Y(1:24)-2,Min_Moment_LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y(1:24)-2,Min_Moment_LC2,'b.','MarkerSize',10)
    hold on
    plot(Y(1:24)-2,Min_Moment_LC3,'r.','MarkerSize',10)
 
%     hold on
%     
%     plot(Y(1:24)-2,Min_Moment_LC4,'k-.','LineWidth',0.8)
%     hold on
%     plot(Y(1:24)-2,Min_Moment_LC5,'b-.','LineWidth',0.8)
%     hold on
%     plot(Y(1:24)-2,Min_Moment_LC6,'r-.','LineWidth',0.8)
    
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('$\Delta$ Bending Moment (Nm)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')
    legend('36 kft +1MC','20 kft +1MC','3 ~ kft +1MC','36 kft -1MC','20 kft -1MC','3 ~ kft -1MC','Interpreter','latex','FontSize',12)

    figure 
    plot(Y(1:24)-2,Max_Torque_LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y(1:24)-2,Max_Torque_LC2,'b.','MarkerSize',10)
    hold on
    plot(Y(1:24)-2,Max_Torque_LC3,'r.','MarkerSize',10)
 
    hold on
    
%     plot(Y(1:24)-2,Max_Torque_LC4,'k-.','LineWidth',0.8)
%     hold on
%     plot(Y(1:24)-2,Max_Torque_LC5,'b-.','LineWidth',0.8)
%     hold on
%     plot(Y(1:24)-2,Max_Torque_LC6,'r-.','LineWidth',0.8)
%     
%     hold on 
    
    plot(Y(1:24)-2,Min_Torque_LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y(1:24)-2,Min_Torque_LC2,'b.','MarkerSize',10)
    hold on
    plot(Y(1:24)-2,Min_Torque_LC3,'r.','MarkerSize',10)
 
%     hold on
%     
%     plot(Y(1:24)-2,Min_Torque_LC4,'k-.','LineWidth',0.8)
%     hold on
%     plot(Y(1:24)-2,Min_Torque_LC5,'b-.','LineWidth',0.8)
%     hold on
%     plot(Y(1:24)-2,Min_Torque_LC6,'r-.','LineWidth',0.8)
    
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('$\Delta$ Torque (Nm)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')
    legend('36 kft +1MC','20 kft +1MC','3 ~ kft +1MC','36 kft -1MC','20 kft -1MC','3 ~ kft -1MC','Interpreter','latex','FontSize',12)


    figure 
    plot(Y(1:24)-2,Max_Shear_LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y(1:24)-2,Max_Shear_LC2,'b.','MarkerSize',10)
    hold on
    plot(Y(1:24)-2,Max_Shear_LC3,'r.','MarkerSize',10)
 
    hold on
    
    plot(Y(1:24)-2,Min_Shear_LC1,'k.','MarkerSize',10)
    hold on 
    plot(Y(1:24)-2,Min_Shear_LC2,'b.','MarkerSize',10)
    hold on
    plot(Y(1:24)-2,Min_Shear_LC3,'r.','MarkerSize',10)
 
    
    xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
    ylabel('$\Delta$ Vertical shear force (Nm)','FontSize',12,'Interpreter','latex')
    set(gcf,'color','w')
    legend('36 kft +1MC','20 kft +1MC','3 ~ kft +1MC','36 kft -1MC','20 kft -1MC','3 ~ kft -1MC','Interpreter','latex','FontSize',12)