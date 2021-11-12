%% Wing-box mass vs FWT eta

% AR 10

AR_10_Eta60_Structure=load('D22_results/Res_AR10_Eta_60_Structure');
AR_10_Eta70_Structure=load('D22_results/Res_AR10_Eta_70_Structure');
AR_10_Eta80_Structure=load('D22_results/Res_AR10_Eta_80_Structure');
AR_10_Eta90_Structure=load('D22_results/Res_AR10_Eta_90_Structure');
AR_10_Eta100_Structure=load('D22_results/Res_AR10_Eta_100_Structure');


% AR 13

AR_13_Eta60_Structure=load('D22_results/Res_AR13_Eta_60_Structure');
AR_13_Eta70_Structure=load('D22_results/Res_AR13_Eta_70_Structure');
AR_13_Eta80_Structure=load('D22_results/Res_AR13_Eta_80_Structure');
AR_13_Eta90_Structure=load('D22_results/Res_AR13_Eta_90_Structure');
AR_13_Eta100_Structure=load('D22_results/Res_AR13_Eta_100_Structure');

% AR 16

AR_16_Eta60_Structure=load('D22_results/Res_AR16_Eta_60_Structure');
AR_16_Eta70_Structure=load('D22_results/Res_AR16_Eta_70_Structure');
AR_16_Eta80_Structure=load('D22_results/Res_AR16_Eta_80_Structure');
AR_16_Eta90_Structure=load('D22_results/Res_AR16_Eta_90_Structure');
AR_16_Eta100_Structure=load('D22_results/Res_AR16_Eta_100_Structure');


% AR 19

AR_19_Eta60_Structure=load('D22_results/Res_AR19_Eta_60_Structure');
AR_19_Eta70_Structure=load('D22_results/Res_AR19_Eta_70_Structure');
AR_19_Eta80_Structure=load('D22_results/Res_AR19_Eta_80_Structure');
AR_19_Eta90_Structure=load('D22_results/Res_AR19_Eta_90_Structure');
AR_19_Eta100_Structure=load('D22_results/Res_AR19_Eta_100_Structure');


% AR 22
AR_22_Eta60_Structure=load('D22_results/Res_AR22_Eta_60_Structure');
AR_22_Eta70_Structure=load('D22_results/Res_AR22_Eta_70_Structure');
AR_22_Eta80_Structure=load('D22_results/Res_AR22_Eta_80_Structure');
AR_22_Eta90_Structure=load('D22_results/Res_AR22_Eta_90_Structure');
AR_22_Eta100_Structure=load('D22_results/Res_AR22_Eta_100_Structure');



Eta=[40,30,20,10,0];

Wingbox_masses_AR10=[AR_10_Eta60_Structure.Wingbox_property.Mass,AR_10_Eta70_Structure.Wingbox_property.Mass,AR_10_Eta80_Structure.Wingbox_property.Mass,AR_10_Eta90_Structure.Wingbox_property.Mass,AR_10_Eta100_Structure.Wingbox_property.Mass];

Wingbox_masses_AR13=[AR_13_Eta60_Structure.Wingbox_property.Mass,AR_13_Eta70_Structure.Wingbox_property.Mass,AR_13_Eta80_Structure.Wingbox_property.Mass,AR_13_Eta90_Structure.Wingbox_property.Mass,AR_13_Eta100_Structure.Wingbox_property.Mass];

Wingbox_masses_AR16=[AR_16_Eta60_Structure.Wingbox_property.Mass,AR_16_Eta70_Structure.Wingbox_property.Mass,AR_16_Eta80_Structure.Wingbox_property.Mass,AR_16_Eta90_Structure.Wingbox_property.Mass,AR_16_Eta100_Structure.Wingbox_property.Mass];

Wingbox_masses_AR19=[AR_19_Eta60_Structure.Wingbox_property.Mass,AR_19_Eta70_Structure.Wingbox_property.Mass,AR_19_Eta80_Structure.Wingbox_property.Mass,AR_19_Eta90_Structure.Wingbox_property.Mass,AR_19_Eta100_Structure.Wingbox_property.Mass];

Wingbox_masses_AR22=[AR_22_Eta60_Structure.Wingbox_property.Mass,AR_22_Eta70_Structure.Wingbox_property.Mass,AR_22_Eta80_Structure.Wingbox_property.Mass,AR_22_Eta90_Structure.Wingbox_property.Mass,AR_22_Eta100_Structure.Wingbox_property.Mass];

% wing box mass 1

figure 

plot(Eta,Wingbox_masses_AR10,'bs-')
hold on 
plot(Eta,Wingbox_masses_AR13,'ks-')
hold on 
plot(Eta,Wingbox_masses_AR16,'ms-')
hold on 
plot(Eta,Wingbox_masses_AR19,'rs-')
hold on 
plot(Eta,Wingbox_masses_AR22,'gs-')

set(gcf,'Color','w')
xlabel('Fold percentage  ~$\%$','Interpreter','latex','FontSize',12)
ylabel('Wingbox mass (kg)','Interpreter','latex','FontSize',12)

legend('AR=10','AR=13','AR=16','AR=19','Interpreter','latex','FontSize',12)



% wing box mass 2 **

AR_mass=[10,13,16,19,22];

Wingbox_masses_Eta100=[AR_10_Eta100_Structure.Wingbox_property.Mass, AR_13_Eta100_Structure.Wingbox_property.Mass, AR_16_Eta100_Structure.Wingbox_property.Mass, AR_19_Eta100_Structure.Wingbox_property.Mass,AR_22_Eta100_Structure.Wingbox_property.Mass];
Wingbox_masses_Eta90=[AR_10_Eta90_Structure.Wingbox_property.Mass, AR_13_Eta90_Structure.Wingbox_property.Mass, AR_16_Eta90_Structure.Wingbox_property.Mass, AR_19_Eta90_Structure.Wingbox_property.Mass, AR_22_Eta90_Structure.Wingbox_property.Mass];
Wingbox_masses_Eta80=[AR_10_Eta80_Structure.Wingbox_property.Mass, AR_13_Eta80_Structure.Wingbox_property.Mass, AR_16_Eta80_Structure.Wingbox_property.Mass, AR_19_Eta80_Structure.Wingbox_property.Mass, AR_22_Eta80_Structure.Wingbox_property.Mass];
Wingbox_masses_Eta70=[AR_10_Eta70_Structure.Wingbox_property.Mass, AR_13_Eta70_Structure.Wingbox_property.Mass, AR_16_Eta70_Structure.Wingbox_property.Mass, AR_19_Eta70_Structure.Wingbox_property.Mass, AR_22_Eta70_Structure.Wingbox_property.Mass];
Wingbox_masses_Eta60=[AR_10_Eta60_Structure.Wingbox_property.Mass, AR_13_Eta60_Structure.Wingbox_property.Mass, AR_16_Eta60_Structure.Wingbox_property.Mass, AR_19_Eta60_Structure.Wingbox_property.Mass, AR_22_Eta60_Structure.Wingbox_property.Mass];

figure 

plot(AR_mass,Wingbox_masses_Eta100,'rs','MarkerFaceColor','r')
hold on 
plot(AR_mass,Wingbox_masses_Eta90,'ks','MarkerFaceColor','b')
hold on 
plot(AR_mass,Wingbox_masses_Eta80,'ks','MarkerFaceColor','g')
hold on 
plot(AR_mass,Wingbox_masses_Eta70,'bs','MarkerFaceColor','k')
hold on 
plot(AR_mass,Wingbox_masses_Eta60,'ks','MarkerFaceColor','y')
% hold on 
% plot(AR_mass,Wingbox_masses_Eta60,'ks','MarkerFaceColor','k')

set(gcf,'Color','w')
xlabel('Wing aspect ratio','Interpreter','latex','FontSize',12)
ylabel('Wingbox mass (kg)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)




% stiffness distributions Ixx **

Y_nodes_eta100=AR_19_Eta100_Structure.Wingbox_property.NodesY;

Y_nodes_eta90=AR_19_Eta90_Structure.Wingbox_property.NodesY;

Y_nodes_eta80=AR_19_Eta80_Structure.Wingbox_property.NodesY;

Y_nodes_eta70=AR_19_Eta70_Structure.Wingbox_property.NodesY;

Y_nodes_eta60=AR_19_Eta60_Structure.Wingbox_property.NodesY;


Wingbox_Ixx_AR19_Eta100=AR_19_Eta100_Structure.Wingbox_property.Ixx;

Wingbox_Ixx_AR19_Eta90=AR_19_Eta90_Structure.Wingbox_property.Ixx;

Wingbox_Ixx_AR19_Eta80=AR_19_Eta80_Structure.Wingbox_property.Ixx;

Wingbox_Ixx_AR19_Eta70=AR_19_Eta70_Structure.Wingbox_property.Ixx;

Wingbox_Ixx_AR19_Eta60=AR_19_Eta60_Structure.Wingbox_property.Ixx;

E=70e9;

figure 

plot(Y_nodes_eta100-2,Wingbox_Ixx_AR19_Eta100*E,'ks-','MarkerFaceColor','r')
hold on 
plot(Y_nodes_eta90-2,Wingbox_Ixx_AR19_Eta90*E,'ks-','MarkerFaceColor','b')
hold on 
plot(Y_nodes_eta80-2,Wingbox_Ixx_AR19_Eta80*E,'ks-','MarkerFaceColor','g')
hold on 
plot(Y_nodes_eta70-2,Wingbox_Ixx_AR19_Eta70*E,'ks-','MarkerFaceColor','k')
hold on 
plot(Y_nodes_eta60-2,Wingbox_Ixx_AR19_Eta60*E,'ks-','MarkerFaceColor','y')

set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Out-of-plane bending stiffness (Nm$^2$)','Interpreter','latex','FontSize',12)
axis([0 22.5 0 4e8])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)


% stiffness distributions Izz **

Wingbox_Izz_AR19_Eta100=AR_19_Eta100_Structure.Wingbox_property.Izz;

Wingbox_Izz_AR19_Eta90=AR_19_Eta90_Structure.Wingbox_property.Izz;

Wingbox_Izz_AR19_Eta80=AR_19_Eta80_Structure.Wingbox_property.Izz;

Wingbox_Izz_AR19_Eta70=AR_19_Eta70_Structure.Wingbox_property.Izz;

Wingbox_Izz_AR19_Eta60=AR_19_Eta60_Structure.Wingbox_property.Izz;

E=70e9;

figure 

plot(Y_nodes_eta100-2,Wingbox_Izz_AR19_Eta100*E,'ks-','MarkerFaceColor','r')
hold on 
plot(Y_nodes_eta90-2,Wingbox_Izz_AR19_Eta90*E,'ks-','MarkerFaceColor','b')
hold on 
plot(Y_nodes_eta80-2,Wingbox_Izz_AR19_Eta80*E,'ks-','MarkerFaceColor','g')
hold on 
plot(Y_nodes_eta70-2,Wingbox_Izz_AR19_Eta70*E,'ks-','MarkerFaceColor','k')
hold on 
plot(Y_nodes_eta60-2,Wingbox_Izz_AR19_Eta60*E,'ks-','MarkerFaceColor','y')
axis([0 22.5 0 5e9])
set(gcf,'Color','w')
set(gca,'TickLabelInterpreter','latex')

xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('In-plane bending stiffness (Nm$^2$)','Interpreter','latex','FontSize',12)
set(gca,'FontSize',14)

legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)


%% Polar plot 


% AR 10

AR_10_Eta60_Aero=load('D22_results/Res_AR10_Eta_60_Aero');
AR_10_Eta70_Aero=load('D22_results/Res_AR10_Eta_70_Aero');
AR_10_Eta80_Aero=load('D22_results/Res_AR10_Eta_80_Aero');
AR_10_Eta90_Aero=load('D22_results/Res_AR10_Eta_90_Aero');
AR_10_Eta100_Aero=load('D22_results/Res_AR10_Eta_100_Aero');


% AR 13

AR_13_Eta60_Aero=load('D22_results/Res_AR13_Eta_60_Aero');
AR_13_Eta70_Aero=load('D22_results/Res_AR13_Eta_70_Aero');
AR_13_Eta80_Aero=load('D22_results/Res_AR13_Eta_80_Aero');
AR_13_Eta90_Aero=load('D22_results/Res_AR13_Eta_90_Aero');
AR_13_Eta100_Aero=load('D22_results/Res_AR13_Eta_100_Aero');

% AR 16

AR_16_Eta60_Aero=load('D22_results/Res_AR16_Eta_60_Aero');
AR_16_Eta70_Aero=load('D22_results/Res_AR16_Eta_70_Aero');
AR_16_Eta80_Aero=load('D22_results/Res_AR16_Eta_80_Aero');
AR_16_Eta90_Aero=load('D22_results/Res_AR16_Eta_90_Aero');
AR_16_Eta100_Aero=load('D22_results/Res_AR16_Eta_100_Aero');


% AR 19

AR_19_Eta60_Aero=load('D22_results/Res_AR19_Eta_60_Aero');
AR_19_Eta70_Aero=load('D22_results/Res_AR19_Eta_70_Aero');
AR_19_Eta80_Aero=load('D22_results/Res_AR19_Eta_80_Aero');
AR_19_Eta90_Aero=load('D22_results/Res_AR19_Eta_90_Aero');
AR_19_Eta100_Aero=load('D22_results/Res_AR19_Eta_100_Aero');


% AR 22

AR_22_Eta60_Aero=load('D22_results/Res_AR22_Eta_60_Aero');
AR_22_Eta70_Aero=load('D22_results/Res_AR22_Eta_70_Aero');
AR_22_Eta80_Aero=load('D22_results/Res_AR22_Eta_80_Aero');
AR_22_Eta90_Aero=load('D22_results/Res_AR22_Eta_90_Aero');
AR_22_Eta100_Aero=load('D22_results/Res_AR22_Eta_100_Aero');


% Elliptical 







K=-[AR_10_Eta100_Aero.Aerodynamics.k,AR_13_Eta100_Aero.Aerodynamics.k,AR_16_Eta100_Aero.Aerodynamics.k,AR_19_Eta100_Aero.Aerodynamics.k, AR_22_Eta100_Aero.Aerodynamics.k];

Cl_plot=0:0.04:1.5;
zero_lift_drag=0.02;

Cd_plot_AR10=K(1).* Cl_plot.^2 + zero_lift_drag;
Cd_plot_AR13=K(2).* Cl_plot.^2 + zero_lift_drag;
Cd_plot_AR16=K(3).* Cl_plot.^2 + zero_lift_drag;
Cd_plot_AR19=K(4).* Cl_plot.^2 + zero_lift_drag;
Cd_plot_AR22=K(5).* Cl_plot.^2 + zero_lift_drag;


figure % polar **

plot(Cd_plot_AR10,Cl_plot,'ks','MarkerFaceColor','b')
hold on 
plot(Cd_plot_AR13,Cl_plot,'ks','MarkerFaceColor','y')
hold on 
plot(Cd_plot_AR16,Cl_plot,'ks','MarkerFaceColor','r')
hold on 
plot(Cd_plot_AR19,Cl_plot,'ks','MarkerFaceColor','g')
hold on 
plot(Cd_plot_AR22,Cl_plot,'ks','MarkerFaceColor','k')


set(gcf,'Color','w')
xlabel('Drag coefficient  ~$C_D$','Interpreter','latex','FontSize',12)
ylabel('Lift coefficient ~$C_L$','Interpreter','latex','FontSize',12)

axis([0 0.08 0 1.5])

legend('AR=10','AR=13','AR=16','AR=19','AR=22','Interpreter','latex','FontSize',12)

figure % L/D **

LD_AR10=Cl_plot./Cd_plot_AR10;
LD_AR13=Cl_plot./Cd_plot_AR13;
LD_AR16=Cl_plot./Cd_plot_AR16;
LD_AR19=Cl_plot./Cd_plot_AR19;
LD_AR22=Cl_plot./Cd_plot_AR22;

plot(Cl_plot,LD_AR10,'ks','MarkerFaceColor','b')
hold on
plot(Cl_plot,LD_AR13,'ks','MarkerFaceColor','y')
hold on
plot(Cl_plot,LD_AR16,'ks','MarkerFaceColor','r')
hold on
plot(Cl_plot,LD_AR19,'ks','MarkerFaceColor','g')
hold on
plot(Cl_plot,LD_AR22,'ks','MarkerFaceColor','k')

set(gcf,'Color','w')
xlabel('Lift coefficient  ~$C_L$','Interpreter','latex','FontSize',12)
ylabel('L/D','Interpreter','latex','FontSize',12)

legend('AR=10','AR=13','AR=16','AR=19','AR=22','Interpreter','latex','FontSize',12)




%% Range 

% Mission: full fuel tank with half payload max

% Payload + Fuel 
Mass_payload=25000;
Fuel_Mass=12000;

% OWE AR 10
OEW_AR10_Eta100=40000 + 1200*2 + 2*AR_10_Eta100_Structure.Wingbox_property.Mass;
OEW_AR10_Eta90=40000 + 1200*2 + 2*AR_10_Eta90_Structure.Wingbox_property.Mass;
OEW_AR10_Eta80=40000 + 1200*2 + 2*AR_10_Eta80_Structure.Wingbox_property.Mass;
OEW_AR10_Eta70=40000 + 1200*2 + 2*AR_10_Eta70_Structure.Wingbox_property.Mass;
OEW_AR10_Eta60=40000 + 1200*2 + 2*AR_10_Eta60_Structure.Wingbox_property.Mass;


% OWE AR 13
OEW_AR13_Eta100=40000 + 1200*2 + 2*AR_13_Eta100_Structure.Wingbox_property.Mass;
OEW_AR13_Eta90=40000 + 1200*2 + 2*AR_13_Eta90_Structure.Wingbox_property.Mass;
OEW_AR13_Eta80=40000 + 1200*2 + 2*AR_13_Eta80_Structure.Wingbox_property.Mass;
OEW_AR13_Eta70=40000 + 1200*2 + 2*AR_13_Eta70_Structure.Wingbox_property.Mass;
OEW_AR13_Eta60=40000 + 1200*2 + 2*AR_13_Eta60_Structure.Wingbox_property.Mass;


% OWE AR 16
OEW_AR16_Eta100=40000 + 1200*2 + 2*AR_16_Eta100_Structure.Wingbox_property.Mass;
OEW_AR16_Eta90=40000 + 1200*2 + 2*AR_16_Eta90_Structure.Wingbox_property.Mass;
OEW_AR16_Eta80=40000 + 1200*2 + 2*AR_16_Eta80_Structure.Wingbox_property.Mass;
OEW_AR16_Eta70=40000 + 1200*2 + 2*AR_16_Eta70_Structure.Wingbox_property.Mass;
OEW_AR16_Eta60=40000 + 1200*2 + 2*AR_16_Eta60_Structure.Wingbox_property.Mass;


% OWE AR 19
OEW_AR19_Eta100=40000 + 1200*2 + 2*AR_19_Eta100_Structure.Wingbox_property.Mass;
OEW_AR19_Eta90=40000 + 1200*2 + 2*AR_19_Eta90_Structure.Wingbox_property.Mass;
OEW_AR19_Eta80=40000 + 1200*2 + 2*AR_19_Eta80_Structure.Wingbox_property.Mass;
OEW_AR19_Eta70=40000 + 1200*2 + 2*AR_19_Eta70_Structure.Wingbox_property.Mass;
OEW_AR19_Eta60=40000 + 1200*2 + 2*AR_19_Eta60_Structure.Wingbox_property.Mass;

% OWE AR 22
OEW_AR22_Eta100=40000 + 1200*2 + 2*AR_22_Eta100_Structure.Wingbox_property.Mass;
OEW_AR22_Eta90=40000 + 1200*2 + 2*AR_22_Eta90_Structure.Wingbox_property.Mass;
OEW_AR22_Eta80=40000 + 1200*2 + 2*AR_22_Eta80_Structure.Wingbox_property.Mass;
OEW_AR22_Eta70=40000 + 1200*2 + 2*AR_22_Eta70_Structure.Wingbox_property.Mass;
OEW_AR22_Eta60=40000 + 1200*2 + 2*AR_22_Eta60_Structure.Wingbox_property.Mass;


% TOW AR 10 
TOW_AR10_Eta100=OEW_AR10_Eta100 + Mass_payload + Fuel_Mass;
TOW_AR10_Eta90=OEW_AR10_Eta90 + Mass_payload + Fuel_Mass;
TOW_AR10_Eta80=OEW_AR10_Eta80 + Mass_payload + Fuel_Mass;
TOW_AR10_Eta70=OEW_AR10_Eta70 + Mass_payload + Fuel_Mass;
TOW_AR10_Eta60=OEW_AR10_Eta60 + Mass_payload + Fuel_Mass;

TOW_AR10=[TOW_AR10_Eta100, TOW_AR10_Eta90, TOW_AR10_Eta80, TOW_AR10_Eta70, TOW_AR10_Eta60];

% TOW AR 13 
TOW_AR13_Eta100=OEW_AR13_Eta100 + Mass_payload + Fuel_Mass;
TOW_AR13_Eta90=OEW_AR13_Eta90 + Mass_payload + Fuel_Mass;
TOW_AR13_Eta80=OEW_AR13_Eta80 + Mass_payload + Fuel_Mass;
TOW_AR13_Eta70=OEW_AR13_Eta70 + Mass_payload + Fuel_Mass;
TOW_AR13_Eta60=OEW_AR13_Eta60 + Mass_payload + Fuel_Mass;

TOW_AR13=[TOW_AR13_Eta100, TOW_AR13_Eta90, TOW_AR13_Eta80, TOW_AR13_Eta70, TOW_AR13_Eta60];

% TOW AR 16 
TOW_AR16_Eta100=OEW_AR16_Eta100 + Mass_payload + Fuel_Mass;
TOW_AR16_Eta90=OEW_AR16_Eta90 + Mass_payload + Fuel_Mass;
TOW_AR16_Eta80=OEW_AR16_Eta80 + Mass_payload + Fuel_Mass;
TOW_AR16_Eta70=OEW_AR16_Eta70 + Mass_payload + Fuel_Mass;
TOW_AR16_Eta60=OEW_AR16_Eta60 + Mass_payload + Fuel_Mass;

TOW_AR16=[TOW_AR16_Eta100, TOW_AR16_Eta90, TOW_AR16_Eta80, TOW_AR16_Eta70, TOW_AR16_Eta60];


% TOW AR 19 
TOW_AR19_Eta100=OEW_AR19_Eta100 + Mass_payload + Fuel_Mass;
TOW_AR19_Eta90=OEW_AR19_Eta90 + Mass_payload + Fuel_Mass;
TOW_AR19_Eta80=OEW_AR19_Eta80 + Mass_payload + Fuel_Mass;
TOW_AR19_Eta70=OEW_AR19_Eta70 + Mass_payload + Fuel_Mass;
TOW_AR19_Eta60=OEW_AR19_Eta60 + Mass_payload + Fuel_Mass;

TOW_AR19=[TOW_AR19_Eta100, TOW_AR19_Eta90, TOW_AR19_Eta80, TOW_AR19_Eta70, TOW_AR19_Eta60];

% TOW AR 22 
TOW_AR22_Eta100=OEW_AR22_Eta100 + Mass_payload + Fuel_Mass;
TOW_AR22_Eta90=OEW_AR22_Eta90 + Mass_payload + Fuel_Mass;
TOW_AR22_Eta80=OEW_AR22_Eta80 + Mass_payload + Fuel_Mass;
TOW_AR22_Eta70=OEW_AR22_Eta70 + Mass_payload + Fuel_Mass;
TOW_AR22_Eta60=OEW_AR22_Eta60 + Mass_payload + Fuel_Mass;

TOW_AR22=[TOW_AR22_Eta100, TOW_AR22_Eta90, TOW_AR22_Eta80, TOW_AR22_Eta70, TOW_AR22_Eta60];

fwt_eta=[0,10,20,30,40];

% Take off weight vs folding eta
figure 

plot(fwt_eta,TOW_AR10,'bs-','MarkerFaceColor','b')
hold on 
plot(fwt_eta,TOW_AR13,'ks-','MarkerFaceColor','y')
hold on 
plot(fwt_eta,TOW_AR16,'rs-','MarkerFaceColor','r')
hold on 
plot(fwt_eta,TOW_AR19,'ks-','MarkerFaceColor','g')
hold on 
plot(fwt_eta,TOW_AR22,'ks-','MarkerFaceColor','k')

set(gcf,'Color','w')
xlabel('Fold percentage  ~$\%$','Interpreter','latex','FontSize',12)
ylabel('Take off mass (kg)','Interpreter','latex','FontSize',12)

legend('AR=10','AR=13','AR=16','AR=19','AR=22','Interpreter','latex','FontSize',12)



% Take off weight vs AR ** 
figure 

AR_TOW=[10,13,16,19,22];

plot(AR_TOW,[TOW_AR10_Eta100,TOW_AR13_Eta100,TOW_AR16_Eta100,TOW_AR19_Eta100,TOW_AR22_Eta100],'ks','MarkerFaceColor','r')
hold on 
plot(AR_TOW,[TOW_AR10_Eta90,TOW_AR13_Eta90,TOW_AR16_Eta90,TOW_AR19_Eta90, TOW_AR22_Eta90],'ks','MarkerFaceColor','b')
hold on 
plot(AR_TOW,[TOW_AR10_Eta80,TOW_AR13_Eta80,TOW_AR16_Eta80,TOW_AR19_Eta80, TOW_AR22_Eta80],'ks','MarkerFaceColor','g')
hold on 
plot(AR_TOW,[TOW_AR10_Eta70,TOW_AR13_Eta70,TOW_AR16_Eta70,TOW_AR19_Eta70,TOW_AR22_Eta70],'ks','MarkerFaceColor','k')
hold on 
plot(AR_TOW,[TOW_AR10_Eta60,TOW_AR13_Eta60,TOW_AR16_Eta60,TOW_AR19_Eta60,TOW_AR22_Eta60],'ks','MarkerFaceColor','y')

set(gcf,'Color','w')
xlabel('Wing aspect ratio','Interpreter','latex','FontSize',12)
ylabel('Take off mass (kg)','Interpreter','latex','FontSize',12)

legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)




c_factor=0.95; % cuz of the tail not included

Cl_mission_AR10=c_factor*9.8*TOW_AR10/(10000*126);

Cl_mission_AR13=c_factor*9.8*TOW_AR13/(10000*126);

Cl_mission_AR16=c_factor*9.8*TOW_AR16/(10000*126);

Cl_mission_AR19=c_factor*9.8*TOW_AR19/(10000*126);

Cl_mission_AR22=c_factor*9.8*TOW_AR22/(10000*126);


Cd_mission_AR10=K(1)*Cl_mission_AR10.^2 + 0.02;

Cd_mission_AR13=K(2)*Cl_mission_AR13.^2 + 0.02;

Cd_mission_AR16=K(3)*Cl_mission_AR16.^2 + 0.02;

Cd_mission_AR19=K(4)*Cl_mission_AR19.^2 + 0.02;

Cd_mission_AR22=K(5)*Cl_mission_AR22.^2 + 0.02;


LD_AR10=Cl_mission_AR10./Cd_mission_AR10;

LD_AR13=Cl_mission_AR13./Cd_mission_AR13;

LD_AR16=Cl_mission_AR16./Cd_mission_AR16;

LD_AR19=Cl_mission_AR19./Cd_mission_AR19;

LD_AR22=Cl_mission_AR22./Cd_mission_AR22;


V_cruise=230*3.6;
Isp=1/0.6;

Wi_AR10=TOW_AR10;
Wf_AR10=TOW_AR10-Fuel_Mass;

Wi_AR13=TOW_AR13;
Wf_AR13=TOW_AR13-Fuel_Mass;

Wi_AR16=TOW_AR16;
Wf_AR16=TOW_AR16-Fuel_Mass;

Wi_AR19=TOW_AR19;
Wf_AR19=TOW_AR19-Fuel_Mass;

Wi_AR22=TOW_AR22;
Wf_AR22=TOW_AR22-Fuel_Mass;


Range_AR10=V_cruise.*LD_AR10.*Isp.*log(Wi_AR10./Wf_AR10);
Range_AR13=V_cruise.*LD_AR13.*Isp.*log(Wi_AR13./Wf_AR13);
Range_AR16=V_cruise.*LD_AR16.*Isp.*log(Wi_AR16./Wf_AR16);
Range_AR19=V_cruise.*LD_AR19.*Isp.*log(Wi_AR19./Wf_AR19);
Range_AR22=V_cruise.*LD_AR22.*Isp.*log(Wi_AR22./Wf_AR22);

AR=[10,13,16,19,22];

Range_eta100=[Range_AR10(1),Range_AR13(1),Range_AR16(1),Range_AR19(1),Range_AR22(1)];
Range_eta90=[Range_AR10(2),Range_AR13(2),Range_AR16(2),Range_AR19(2),Range_AR22(2)];
Range_eta80=[Range_AR10(3),Range_AR13(3),Range_AR16(3),Range_AR19(3),Range_AR22(3)];
Range_eta70=[Range_AR10(4),Range_AR13(4),Range_AR16(4),Range_AR19(4),Range_AR22(4)];
Range_eta60=[Range_AR10(5),Range_AR13(5),Range_AR16(5),Range_AR19(5),Range_AR22(5)];

figure % Range **

plot(AR,Range_eta100,'k-s','MarkerFaceColor','r')
hold on 
plot(AR,Range_eta90,'k-s','MarkerFaceColor','b')
hold on 
plot(AR,Range_eta80,'k-s','MarkerFaceColor','g')
hold on 
plot(AR,Range_eta70,'k-s','MarkerFaceColor','k')
hold on 
plot(AR,Range_eta60,'k-s','MarkerFaceColor','y')

set(gcf,'Color','w')
xlabel('Wing aspect ratio','Interpreter','latex','FontSize',12)
ylabel('Breguet range (km)','Interpreter','latex','FontSize',12)

legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)


figure % L_D vs AR

AR_LD=[10,13,16,19,22];

LD1=[LD_AR10(1),LD_AR13(1),LD_AR16(1),LD_AR19(1),LD_AR22(1)];
LD2=[LD_AR10(2),LD_AR13(2),LD_AR16(2),LD_AR19(2),LD_AR22(2)];
LD3=[LD_AR10(3),LD_AR13(3),LD_AR16(3),LD_AR19(3),LD_AR22(3)];
LD4=[LD_AR10(4),LD_AR13(4),LD_AR16(4),LD_AR19(4),LD_AR22(4)];
LD5=[LD_AR10(5),LD_AR13(5),LD_AR16(5),LD_AR19(5),LD_AR22(5)];

plot(AR_LD,LD1,'k-s','MarkerFaceColor','r')
hold on 
plot(AR_LD,LD2,'k-s','MarkerFaceColor','b')
hold on 
plot(AR_LD,LD3,'k-s','MarkerFaceColor','g')
hold on 
plot(AR_LD,LD4,'k-s','MarkerFaceColor','k')
hold on 
plot(AR_LD,LD5,'k-s','MarkerFaceColor','y')

set(gcf,'Color','w')
xlabel('Wing aspect ratio','Interpreter','latex','FontSize',12)
ylabel('Lift to drag ratio (L/D)','Interpreter','latex','FontSize',12)
legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)


%% peak range prediction 




%% Gust responses

% AR 10

AR_10_Eta60_Loads=load('D22_results/Res_AR10_Eta_60_Loads');
AR_10_Eta70_Loads=load('D22_results/Res_AR10_Eta_70_Loads');
AR_10_Eta80_Loads=load('D22_results/Res_AR10_Eta_80_Loads');
AR_10_Eta90_Loads=load('D22_results/Res_AR10_Eta_90_Loads');
AR_10_Eta100_Loads=load('D22_results/Res_AR10_Eta_100_Loads');


% AR 13

AR_13_Eta60_Loads=load('D22_results/Res_AR13_Eta_60_Loads');
AR_13_Eta70_Loads=load('D22_results/Res_AR13_Eta_70_Loads');
AR_13_Eta80_Loads=load('D22_results/Res_AR13_Eta_80_Loads');
AR_13_Eta90_Loads=load('D22_results/Res_AR13_Eta_90_Loads');
AR_13_Eta100_Loads=load('D22_results/Res_AR13_Eta_100_Loads');

% AR 16

AR_16_Eta60_Loads=load('D22_results/Res_AR16_Eta_60_Loads');
AR_16_Eta70_Loads=load('D22_results/Res_AR16_Eta_70_Loads');
AR_16_Eta80_Loads=load('D22_results/Res_AR16_Eta_80_Loads');
AR_16_Eta90_Loads=load('D22_results/Res_AR16_Eta_90_Loads');
AR_16_Eta100_Loads=load('D22_results/Res_AR16_Eta_100_Loads');


% % AR 19
% 
% AR_19_Eta60_Loads=load('D22_results/Res_AR19_Eta_60_Loads');
% AR_19_Eta70_Loads=load('D22_results/Res_AR19_Eta_70_Loads');
% AR_19_Eta80_Loads=load('D22_results/Res_AR19_Eta_80_Loads');
% AR_19_Eta90_Loads=load('D22_results/Res_AR19_Eta_90_Loads');
% AR_19_Eta100_Loads=load('D22_results/Res_AR19_Eta_100_Loads');





% gust wing root 

Gust_length=linspace(18,214,7);


AR13_Eta100_Root_max=AR_13_Eta100_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR13_Eta100_Root_min=AR_13_Eta100_Loads.Loads.Root_Delta.Moment_Min.LC3;

AR13_Eta90_Root_max=AR_13_Eta90_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR13_Eta90_Root_min=AR_13_Eta90_Loads.Loads.Root_Delta.Moment_Min.LC3;

AR13_Eta80_Root_max=AR_13_Eta80_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR13_Eta80_Root_min=AR_13_Eta80_Loads.Loads.Root_Delta.Moment_Min.LC3;

AR13_Eta70_Root_max=AR_13_Eta70_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR13_Eta70_Root_min=AR_13_Eta70_Loads.Loads.Root_Delta.Moment_Min.LC3;

AR13_Eta60_Root_max=AR_13_Eta60_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR13_Eta60_Root_min=AR_13_Eta60_Loads.Loads.Root_Delta.Moment_Min.LC3;




AR16_Eta100_Root_max=AR_16_Eta100_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR16_Eta100_Root_min=AR_16_Eta100_Loads.Loads.Root_Delta.Moment_Min.LC3;

AR16_Eta90_Root_max=AR_16_Eta90_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR16_Eta90_Root_min=AR_16_Eta90_Loads.Loads.Root_Delta.Moment_Min.LC3;

AR16_Eta80_Root_max=AR_16_Eta80_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR16_Eta80_Root_min=AR_16_Eta80_Loads.Loads.Root_Delta.Moment_Min.LC3;

AR16_Eta70_Root_max=AR_16_Eta70_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR16_Eta70_Root_min=AR_16_Eta70_Loads.Loads.Root_Delta.Moment_Min.LC3;

AR16_Eta60_Root_max=AR_16_Eta60_Loads.Loads.Root_Delta.Moment_Max.LC3;
AR16_Eta60_Root_min=AR_16_Eta60_Loads.Loads.Root_Delta.Moment_Min.LC3;



figure 
plot(Gust_length,AR16_Eta100_Root_max,'r-s')
hold on 
plot(Gust_length,AR16_Eta100_Root_min,'b-s')

hold on

plot(Gust_length,AR16_Eta80_Root_max,'r-o')
hold on 
plot(Gust_length,AR16_Eta80_Root_min,'b-o')

hold on

plot(Gust_length,AR16_Eta60_Root_max,'r-v')
hold on 
plot(Gust_length,AR16_Eta60_Root_min,'b-v')


figure 
plot(Gust_length,AR13_Eta100_Root_max,'r-s')
hold on 
plot(Gust_length,AR13_Eta100_Root_min,'b-s')

hold on

plot(Gust_length,AR13_Eta80_Root_max,'r-o')
hold on 
plot(Gust_length,AR13_Eta80_Root_min,'b-o')

hold on

plot(Gust_length,AR13_Eta60_Root_max,'r-v')
hold on 
plot(Gust_length,AR13_Eta60_Root_min,'b-v')


% gust wing delta

AR16_Eta100_Y_load=AR_16_Eta100_Loads.Loads.Y;
AR16_Eta90_Y_load=AR_16_Eta90_Loads.Loads.Y;
AR16_Eta80_Y_load=AR_16_Eta80_Loads.Loads.Y;
AR16_Eta70_Y_load=AR_16_Eta70_Loads.Loads.Y;
AR16_Eta60_Y_load=AR_16_Eta60_Loads.Loads.Y;

% bending moment
AR16_Eta100_Wing_max=AR_16_Eta100_Loads.Loads.Wing_Delta.Moment_Max.LC3;
AR16_Eta100_Wing_min=AR_16_Eta100_Loads.Loads.Wing_Delta.Moment_Min.LC3;

AR16_Eta90_Wing_max=AR_16_Eta90_Loads.Loads.Wing_Delta.Moment_Max.LC3;
AR16_Eta90_Wing_min=AR_16_Eta90_Loads.Loads.Wing_Delta.Moment_Min.LC3;

AR16_Eta80_Wing_max=AR_16_Eta80_Loads.Loads.Wing_Delta.Moment_Max.LC3;
AR16_Eta80_Wing_min=AR_16_Eta80_Loads.Loads.Wing_Delta.Moment_Min.LC3;

AR16_Eta70_Wing_max=AR_16_Eta70_Loads.Loads.Wing_Delta.Moment_Max.LC3;
AR16_Eta70_Wing_min=AR_16_Eta70_Loads.Loads.Wing_Delta.Moment_Min.LC3;

AR16_Eta60_Wing_max=AR_16_Eta60_Loads.Loads.Wing_Delta.Moment_Max.LC3;
AR16_Eta60_Wing_min=AR_16_Eta60_Loads.Loads.Wing_Delta.Moment_Min.LC3;

% shear
AR16_Eta100_Wing_shear_max=AR_16_Eta100_Loads.Loads.Wing_Delta.Shear_Max.LC3;
AR16_Eta100_Wing_shear_min=AR_16_Eta100_Loads.Loads.Wing_Delta.Shear_Min.LC3;

AR16_Eta90_Wing_shear_max=AR_16_Eta90_Loads.Loads.Wing_Delta.Shear_Max.LC3;
AR16_Eta90_Wing_shear_min=AR_16_Eta90_Loads.Loads.Wing_Delta.Shear_Min.LC3;

AR16_Eta80_Wing_shear_max=AR_16_Eta80_Loads.Loads.Wing_Delta.Shear_Max.LC3;
AR16_Eta80_Wing_shear_min=AR_16_Eta80_Loads.Loads.Wing_Delta.Shear_Min.LC3;

AR16_Eta70_Wing_shear_max=AR_16_Eta70_Loads.Loads.Wing_Delta.Shear_Max.LC3;
AR16_Eta70_Wing_shear_min=AR_16_Eta70_Loads.Loads.Wing_Delta.Shear_Min.LC3;

AR16_Eta60_Wing_shear_max=AR_16_Eta60_Loads.Loads.Wing_Delta.Shear_Max.LC3;
AR16_Eta60_Wing_shear_min=AR_16_Eta60_Loads.Loads.Wing_Delta.Shear_Min.LC3;


% torque
AR16_Eta100_Wing_torque_max=AR_16_Eta100_Loads.Loads.Wing_Delta.Torque_Max.LC3;
AR16_Eta100_Wing_torque_min=AR_16_Eta100_Loads.Loads.Wing_Delta.Torque_Min.LC3;

AR16_Eta90_Wing_torque_max=AR_16_Eta90_Loads.Loads.Wing_Delta.Torque_Max.LC3;
AR16_Eta90_Wing_torque_min=AR_16_Eta90_Loads.Loads.Wing_Delta.Torque_Min.LC3;

AR16_Eta80_Wing_torque_max=AR_16_Eta80_Loads.Loads.Wing_Delta.Torque_Max.LC3;
AR16_Eta80_Wing_torque_min=AR_16_Eta80_Loads.Loads.Wing_Delta.Torque_Min.LC3;

AR16_Eta70_Wing_torque_max=AR_16_Eta70_Loads.Loads.Wing_Delta.Torque_Max.LC3;
AR16_Eta70_Wing_torque_min=AR_16_Eta70_Loads.Loads.Wing_Delta.Torque_Min.LC3;

AR16_Eta60_Wing_torque_max=AR_16_Eta60_Loads.Loads.Wing_Delta.Torque_Max.LC3;
AR16_Eta60_Wing_torque_min=AR_16_Eta60_Loads.Loads.Wing_Delta.Torque_Min.LC3;



figure % gust load allevation 1**

plot(AR16_Eta100_Y_load(1:end-1)-2,AR16_Eta100_Wing_max,'k-s','MarkerFaceColor','r');
hold on 
plot(AR16_Eta100_Y_load(1:end-1)-2,AR16_Eta100_Wing_min,'k-s','MarkerFaceColor','r')

hold on 

plot(AR16_Eta80_Y_load(1:end-1)-2,AR16_Eta80_Wing_max,'k-s','MarkerFaceColor','g')
hold on 
plot(AR16_Eta80_Y_load(1:end-1)-2,AR16_Eta80_Wing_min,'k-s','MarkerFaceColor','g')

hold on 

plot(AR16_Eta60_Y_load(1:end-1)-2,AR16_Eta60_Wing_max,'k-s','MarkerFaceColor','y')
hold on 
plot(AR16_Eta60_Y_load(1:end-1)-2,AR16_Eta60_Wing_min,'k-s','MarkerFaceColor','y')

set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('$\Delta$ Moment (Nm)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
legend('No FWT','','$\eta$=20$\%$','','$\eta$=40$\%$','Interpreter','latex','FontSize',12)
axis([0 20.05 -1.5e6 1.5e6]);
% ax.PlotBoxAspectRatio = [1 0.5 0.5]

figure % gust load allevation 3**

plot(AR16_Eta100_Y_load(1:end-1)-2,AR16_Eta100_Wing_torque_max,'k-o','MarkerFaceColor','r')
hold on 
plot(AR16_Eta100_Y_load(1:end-1)-2,AR16_Eta100_Wing_torque_min,'k-o','MarkerFaceColor','r')


hold on 

plot(AR16_Eta80_Y_load(1:end-1)-2,AR16_Eta80_Wing_torque_max,'k-o','MarkerFaceColor','g')
hold on 
plot(AR16_Eta80_Y_load(1:end-1)-2,AR16_Eta80_Wing_torque_min,'k-o','MarkerFaceColor','g')


hold on 

plot(AR16_Eta60_Y_load(1:end-1)-2,AR16_Eta60_Wing_torque_max,'k-o','MarkerFaceColor','y')
hold on 
plot(AR16_Eta60_Y_load(1:end-1)-2,AR16_Eta60_Wing_torque_min,'k-o','MarkerFaceColor','y')

set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('$\Delta$ Torque (Nm)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
legend('No FWT','','$\eta$=20$\%$','','$\eta$=40$\%$','Interpreter','latex','FontSize',12)
axis([0 20.05 -1.5e5 1.5e5]);


figure % gust load allevation 2**

plot(AR16_Eta100_Y_load(1:end-1)-2,AR16_Eta100_Wing_shear_max,'k-v','MarkerFaceColor','r')
hold on 
plot(AR16_Eta100_Y_load(1:end-1)-2,AR16_Eta100_Wing_shear_min,'k-v','MarkerFaceColor','r')

hold on 

plot(AR16_Eta80_Y_load(1:end-1)-2,AR16_Eta80_Wing_shear_max,'k-v','MarkerFaceColor','g')
hold on 
plot(AR16_Eta80_Y_load(1:end-1)-2,AR16_Eta80_Wing_shear_min,'k-v','MarkerFaceColor','g')

hold on 

plot(AR16_Eta60_Y_load(1:end-1)-2,AR16_Eta60_Wing_shear_max,'k-v','MarkerFaceColor','y')
hold on 
plot(AR16_Eta60_Y_load(1:end-1)-2,AR16_Eta60_Wing_shear_min,'k-v','MarkerFaceColor','y')

set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('$\Delta$ Vertical shear (N)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
legend('No FWT','','$\eta$=20$\%$','','$\eta$=40$\%$','Interpreter','latex','FontSize',12)
axis([0 20.05 -1.5e5 1.5e5]);



%% Cross section Airfoil 

NACA=load('NACA0012.txt')

figure
plot(NACA(:,1),NACA(:,2),'k-','LineWidth',1)

axis equal
set(gcf,'Color','w')



%% Load distribution 

% add hinge locked results 
% 
% AR16_Eta70_Model=load('D22_results\Res_AR16_Eta_70_Model');
% 
% Param=AR16_Eta70_Model.Param;
% 
% run_folder='C:\Git\A321_sizing\hg_codes\results\cruise_condition';
% 
% [CDi,CD0,CL,k,Distribution,Load_distribution]=Drag_Calculation(Param, run_folder);
% 
% 
% % bending moments 
% 
% figure 
% 
% plot(AR_16_Eta70_Loads.Loads.Y-2,AR_16_Eta70_Loads.Loads.Moment.pullup_cruise,'k-s','MarkerFaceColor','r');
% hold on 
% plot(AR_16_Eta70_Loads.Loads.Y-2,AR_16_Eta70_Loads.Loads.Moment.pullup_sealevel,'k-s','MarkerFaceColor','b')
% hold on 
% plot(AR_16_Eta70_Loads.Loads.Y-2,AR_16_Eta70_Loads.Loads.Moment.dive_sealevel,'k-s','MarkerFaceColor','g')
% hold on 
% plot(AR_16_Eta70_Loads.Loads.Y-2,AR_16_Eta70_Loads.Loads.Moment.g_sealevel_gust1,'k-s','MarkerFaceColor','y')
% 
% hold on 
% plot(Load_distribution.Y-2,Load_distribution.Moment_P2,'k-s','MarkerFaceColor','k')
% 
% axis([0 20.5 0 5e6])
% 
% legend('2.5g 36000ft','2.5g 3000ft','g 3000ft + 1MC')


%%  full load cases

% AR16_Eta70_Model=load('D22_results\Res_AR16_Eta_70_Model');
% 
% Param=AR16_Eta70_Model.Param;
% 
% run_folder='C:\Git\A321_sizing\hg_codes\results\AR16_Eta70_LC_Full';
% 
% [Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis(Param, run_folder);
% 
% AR16_Eta70_Load_All.Y=Y_all;
% AR16_Eta70_Load_All.Loads=Load_distribution;
% AR16_Eta70_Load_All.Wing_Delta=Wing_Delta;
% AR16_Eta70_Load_All.Root_Delta=Root_Delta;
% 
% save('D22_results\AR16_Eta_70_Load_All.mat','AR16_Eta70_Load_All')

AR16_load_all=load('D22_results/AR16_Eta_70_Load_All');

Y_all=AR16_load_all.AR16_Eta70_Load_All.Y;

% Bending moment 
figure 

plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC1,'k-s','MarkerFaceColor','r')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC2,'k-s','MarkerFaceColor','b')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC3,'k-s','MarkerFaceColor','g')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC6,'k-s','MarkerFaceColor','y')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC9,'k-s','MarkerFaceColor','m')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.Cruise,'k-s','MarkerFaceColor','k')

% for sizing 
All_moment=[AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC1,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC2,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC3,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC6,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.LC9,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Moment_P2.Cruise];

sizing_moment=max(All_moment,[],2);

hold on 
plot(Y_all-2,sizing_moment,'r--','LineWidth',2.5)


axis([0 20.5 0 5e6])
legend('2.5g (36000ft)','2.5g (3000ft)','-g (3000ft)','g + 1MC (3000ft)','g - 1MC (3000ft)','Cruise (36000ft)','Sizing','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')

xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Bending moment (Nm)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)



% shear forces
figure 

plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC1,'k-s','MarkerFaceColor','r')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC2,'k-s','MarkerFaceColor','b')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC3,'k-s','MarkerFaceColor','g')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC6,'k-s','MarkerFaceColor','y')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC9,'k-s','MarkerFaceColor','m')
hold on 
plot(Y_all-2,AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.Cruise,'k-s','MarkerFaceColor','k')


% for sizing 
All_shear=[AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC1,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC2,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC3,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC6,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.LC9,...
    AR16_load_all.AR16_Eta70_Load_All.Loads.Shear_P2.Cruise];

sizing_shear=max(All_shear,[],2);

hold on 
plot(Y_all-2,sizing_shear,'r--','LineWidth',2.5)


axis([0 20.5 0 6e5])
legend('2.5g (36000ft)','2.5g (3000ft)','-g (3000ft)','g + 1MC (3000ft)','g - 1MC (3000ft)','Cruise (36000ft)','Sizing','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Shear forces (N)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)



% Torque
figure 

plot(Y_all-2,abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC1),'k-s','MarkerFaceColor','r')
hold on 
plot(Y_all-2,abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC2),'k-s','MarkerFaceColor','b')
hold on 
plot(Y_all-2,abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC3),'k-s','MarkerFaceColor','g')
hold on 
plot(Y_all-2,abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC6),'k-s','MarkerFaceColor','y')
hold on 
plot(Y_all-2,abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC9),'k-s','MarkerFaceColor','m')
hold on 
plot(Y_all-2,abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.Cruise),'k-s','MarkerFaceColor','k')


% for sizing 

All_torque=[abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC1),...
    abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC2),...
    abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC3),...
    abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC6),...
    abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.LC9),...
    abs(AR16_load_all.AR16_Eta70_Load_All.Loads.Torque.Cruise)];

sizing_torque=max(All_torque,[],2);

hold on 
plot(Y_all-2,sizing_torque,'r--','LineWidth',2)

legend('2.5g (36000ft)','2.5g (3000ft)','-g (3000ft)','g + 1MC (3000ft)','g - 1MC (3000ft)','Cruise (36000ft)','Sizing','Interpreter','latex','FontSize',12)
axis([0 20.5 0 2.5e5])
set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Torque (Nm)','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)


%% AR 16 load distribution with different Eta

% AR 16

AR_16_Eta60_Loads=load('D22_results/Res_AR16_Eta_60_Loads');
AR_16_Eta70_Loads=load('D22_results/Res_AR16_Eta_70_Loads');
AR_16_Eta80_Loads=load('D22_results/Res_AR16_Eta_80_Loads');
AR_16_Eta90_Loads=load('D22_results/Res_AR16_Eta_90_Loads');
AR_16_Eta100_Loads=load('D22_results/Res_AR16_Eta_100_Loads');

figure %bending moment 

plot(AR_16_Eta100_Loads.Loads.Y-2,AR_16_Eta100_Loads.Loads.Moment.pullup_cruise,'k-s','MarkerFaceColor','r')
hold on 
plot(AR_16_Eta90_Loads.Loads.Y-2,AR_16_Eta90_Loads.Loads.Moment.pullup_cruise,'k-s','MarkerFaceColor','b')
hold on 
plot(AR_16_Eta80_Loads.Loads.Y-2,AR_16_Eta80_Loads.Loads.Moment.pullup_cruise,'k-s','MarkerFaceColor','g')
hold on 
plot(AR_16_Eta70_Loads.Loads.Y-2,AR_16_Eta70_Loads.Loads.Moment.pullup_cruise,'k-s','MarkerFaceColor','k')
hold on 
plot(AR_16_Eta60_Loads.Loads.Y-2,AR_16_Eta60_Loads.Loads.Moment.pullup_cruise,'k-s','MarkerFaceColor','y')
axis([0 20.5 0 6e6])
set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Bending moment (Nm)','Interpreter','latex','FontSize',12)
legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)

figure %shear 

plot(AR_16_Eta100_Loads.Loads.Y-2,AR_16_Eta100_Loads.Loads.Shear.pullup_cruise,'k-s','MarkerFaceColor','r')
hold on 
plot(AR_16_Eta90_Loads.Loads.Y-2,AR_16_Eta90_Loads.Loads.Shear.pullup_cruise,'k-s','MarkerFaceColor','b')
hold on 
plot(AR_16_Eta80_Loads.Loads.Y-2,AR_16_Eta80_Loads.Loads.Shear.pullup_cruise,'k-s','MarkerFaceColor','g')
hold on 
plot(AR_16_Eta70_Loads.Loads.Y-2,AR_16_Eta70_Loads.Loads.Shear.pullup_cruise,'k-s','MarkerFaceColor','k')
hold on 
plot(AR_16_Eta60_Loads.Loads.Y-2,AR_16_Eta60_Loads.Loads.Shear.pullup_cruise,'k-s','MarkerFaceColor','y')
axis([0 20.5 0 6e5])
set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Vertical shear force (N)','Interpreter','latex','FontSize',12)
legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)


figure %torque

plot(AR_16_Eta100_Loads.Loads.Y-2,AR_16_Eta100_Loads.Loads.Torque.pullup_cruise,'k-s','MarkerFaceColor','r')
hold on 
plot(AR_16_Eta90_Loads.Loads.Y-2,AR_16_Eta90_Loads.Loads.Torque.pullup_cruise,'k-s','MarkerFaceColor','b')
hold on 
plot(AR_16_Eta80_Loads.Loads.Y-2,AR_16_Eta80_Loads.Loads.Torque.pullup_cruise,'k-s','MarkerFaceColor','g')
hold on 
plot(AR_16_Eta70_Loads.Loads.Y-2,AR_16_Eta70_Loads.Loads.Torque.pullup_cruise,'k-s','MarkerFaceColor','k')
hold on 
plot(AR_16_Eta60_Loads.Loads.Y-2,AR_16_Eta60_Loads.Loads.Torque.pullup_cruise,'k-s','MarkerFaceColor','y')
axis([0 20.5 -2.5e5 2.5e5])
set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Torque (Nm)','Interpreter','latex','FontSize',12)
legend('No FWT','$\eta$=10$\%$','$\eta$=20$\%$','$\eta$=30$\%$','$\eta$=40$\%$','Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)

%% range test

Cl_test=0.1:0.1:1.5;

Cd0_test=0.018;
Cdi_test=Cl_test.^2/(10*pi);

LD_test=Cl_test./(Cd0_test + Cdi_test);

figure 

plot(Cl_test,LD_test,'s-')

%% Run trim analysis to work out coast angle 

run_folder='C:\Git\A321_sizing\hg_codes\results\coast_angles';

Model_data=load('D22_results\Res_AR10_Eta_90_Model');

Param=Model_data.Param;

% clear Jig shape
Param.Wing.Jig_Twist=[0,0];
Param.Wing.Jig_Eta=[0,1];

Param.Wing.FWT=[0,0];
Param.Wing.FWT_Eta=[0,1];

[CDi,CD0,CL,k,Distribution,Load_distribution, Displacements]=Static_trim(Param, 2.5, 'AR10_Eta_90_manuever_no_jig', 'lock off', run_folder);

% save('Res_AR16_Eta_70_cruise_Displacement.mat','Displacements')

figure 

plot(Displacements.Y_Data,Displacements.Z_all,'bs')


%% plot results - cruise
% Run_folder=strcat(run_folder,'\hinge_loacked');
model = mni.import_matran(fullfile(run_folder,'AR16_Eta_70_manuvere.dat'));
model.draw

%extract the data
f06 = mni.result.f06(fullfile(run_folder,'AR16_Eta_70_manuvere.f06'));
res_disp =  f06.read_disp;
res_aeroP = f06.read_aeroP;
res_aeroF = f06.read_aeroF;

% apply deformation result
[~,i] = ismember(model.GRID.GID,res_disp.GP);
model.GRID.Deformation = [res_disp.dX(:,i);res_disp.dY(:,i);res_disp.dZ(:,i)];

% apply aero pressure
model.CAERO1.PanelPressure = res_aeroP.Cp;

%apply aero forces
f = [res_aeroF.aeroFx;res_aeroF.aeroFy;res_aeroF.aeroFz;...
    res_aeroF.aeroMx;res_aeroF.aeroMy;res_aeroF.aeroMz];
model.CAERO1.PanelForce = f';

% update the plot to apply deformations and aero pressures + forces
model.update('Scale',1)

view([-90 0])
%% plot results - Manuever
% Run_folder=strcat(run_folder,'\hinge_loacked');
model = mni.import_matran(fullfile(run_folder,'AR16_Eta_70_cruise.dat'));
model.draw

%extract the data
f06 = mni.result.f06(fullfile(run_folder,'AR16_Eta_70_cruise.f06'));
res_disp =  f06.read_disp;
res_aeroP = f06.read_aeroP;
res_aeroF = f06.read_aeroF;

% apply deformation result
[~,i] = ismember(model.GRID.GID,res_disp.GP);
model.GRID.Deformation = [res_disp.dX(:,i);res_disp.dY(:,i);res_disp.dZ(:,i)];

% apply aero pressure
model.CAERO1.PanelPressure = res_aeroP.Cp;

%apply aero forces
f = [res_aeroF.aeroFx;res_aeroF.aeroFy;res_aeroF.aeroFz;...
    res_aeroF.aeroMx;res_aeroF.aeroMy;res_aeroF.aeroMz];
model.CAERO1.PanelForce = f';

% update the plot to apply deformations and aero pressures + forces
model.update('Scale',1)
view([-90 0])

set(gcf,'Color','w')
set(gca,'FontSize',18)
axis([-10 50 -25 25 -2 10])
daspect([1 1 1])

%% Coast angle 

eta_coast=[10,20,30,40];

coast_AR10=[26.31,26.52,24.19,25];
coast_AR13=[34.56,31.14,28.09,24.63];
coast_AR16=[40.99,36.6,33.1,29.1];
coast_AR19=[45.6,42.9,38.4,33.5];
coast_AR22=[50.7,47.95,43.78,38.1];

figure 

% plot(eta_coast,coast_AR10,'k-s','MarkerFaceColor','r')
% hold on 
% plot(eta_coast,coast_AR13,'k-s','MarkerFaceColor','b')
% hold on 
plot(eta_coast,coast_AR16,'b-s','MarkerFaceColor','b')
% hold on 
% plot(eta_coast,coast_AR19,'k-s','MarkerFaceColor','y')
% hold on 
% plot(eta_coast,coast_AR22,'k-s','MarkerFaceColor','k')

set(gcf,'Color','w')
xlabel('$\eta$','Interpreter','latex','FontSize',12)
ylabel('Coast angle ($^{\circ}$)','Interpreter','latex','FontSize',12)
% legend('No FWT','10$\%$ FWT','20$\%$ FWT','30$\%$ FWT','40$\%$ FWT','Interpreter','latex','FontSize',12)

