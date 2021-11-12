% Load wing masses
% 
% % AR = 10
% AR10_eta100=load('NRes_A126_AR10_Eta_100_Model');
% AR10_eta90=load('NRes_A126_AR10_Eta_90_Model');
% AR10_eta80=load('NRes_A126_AR10_Eta_80_Model');
% AR10_eta70=load('NRes_A126_AR10_Eta_70_Model');
% AR10_eta60=load('NRes_A126_AR10_Eta_60_Model');
% 
% AR10_eta100_mass=AR10_eta100.Param.Sized_Masses.Total_mass;
% AR10_eta90_mass=AR10_eta90.Param.Sized_Masses.Total_mass;
% AR10_eta80_mass=AR10_eta80.Param.Sized_Masses.Total_mass;
% AR10_eta70_mass=AR10_eta70.Param.Sized_Masses.Total_mass;
% AR10_eta60_mass=AR10_eta60.Param.Sized_Masses.Total_mass;
% 
% AR_10_masses=[AR10_eta100_mass,AR10_eta90_mass,AR10_eta80_mass,AR10_eta70_mass,AR10_eta60_mass];
% 
% 
% % AR = 13
% AR13_eta100=load('NRes_A126_AR13_Eta_100_Model');
% AR13_eta90=load('NRes_A126_AR13_Eta_90_Model');
% AR13_eta80=load('NRes_A126_AR13_Eta_80_Model');
% AR13_eta70=load('NRes_A126_AR13_Eta_70_Model');
% AR13_eta60=load('NRes_A126_AR13_Eta_60_Model');
% 
% AR13_eta100_mass=AR13_eta100.Param.Sized_Masses.Total_mass;
% AR13_eta90_mass=AR13_eta90.Param.Sized_Masses.Total_mass;
% AR13_eta80_mass=AR13_eta80.Param.Sized_Masses.Total_mass;
% AR13_eta70_mass=AR13_eta70.Param.Sized_Masses.Total_mass;
% AR13_eta60_mass=AR13_eta60.Param.Sized_Masses.Total_mass;
% 
% AR_13_masses=[AR13_eta100_mass,AR13_eta90_mass,AR13_eta80_mass,AR13_eta70_mass,AR13_eta60_mass];
% 
% 
% 
% % AR = 16
% AR16_eta100=load('NRes_A126_AR16_Eta_100_Model');
% AR16_eta90=load('NRes_A126_AR16_Eta_90_Model');
% AR16_eta80=load('NRes_A126_AR16_Eta_80_Model');
% AR16_eta70=load('NRes_A126_AR16_Eta_70_Model');
% AR16_eta60=load('NRes_A126_AR16_Eta_60_Model');
% 
% AR16_eta100_mass=AR16_eta100.Param.Sized_Masses.Total_mass;
% AR16_eta90_mass=AR16_eta90.Param.Sized_Masses.Total_mass;
% AR16_eta80_mass=AR16_eta80.Param.Sized_Masses.Total_mass;
% AR16_eta70_mass=AR16_eta70.Param.Sized_Masses.Total_mass;
% AR16_eta60_mass=AR16_eta60.Param.Sized_Masses.Total_mass;
% 
% AR_16_masses=[AR16_eta100_mass,AR16_eta90_mass,AR16_eta80_mass,AR16_eta70_mass,AR16_eta60_mass];


% % AR = 22
% AR22_eta100=load('NRes_A126_AR22_Eta_100_Model');
% AR22_eta90=load('NRes_A126_AR22_Eta_90_Model');
% AR22_eta80=load('NRes_A126_AR22_Eta_80_Model');
% AR22_eta70=load('NRes_A126_AR22_Eta_70_Model');
% AR22_eta60=load('NRes_A126_AR22_Eta_60_Model');
% 
% AR22_eta100_mass=AR22_eta100.Param.Sized_Masses.Total_mass;
% AR22_eta90_mass=AR22_eta90.Param.Sized_Masses.Total_mass;
% AR22_eta80_mass=AR22_eta80.Param.Sized_Masses.Total_mass;
% AR22_eta70_mass=AR22_eta70.Param.Sized_Masses.Total_mass;
% AR22_eta60_mass=AR22_eta60.Param.Sized_Masses.Total_mass;
% 
% AR_22_masses=[AR22_eta100_mass,AR22_eta90_mass,AR22_eta80_mass,AR22_eta70_mass,AR22_eta60_mass];



%% Wing weight 

% AR = 19
AR19_eta100=load('NRes_A126_AR19_Eta_100_Model');
AR19_eta90=load('NRes_A126_AR19_Eta_90_Model');
AR19_eta80=load('NRes_A126_AR19_Eta_80_Model');
AR19_eta70=load('NRes_A126_AR19_Eta_70_Model');
AR19_eta60=load('NRes_A126_AR19_Eta_60_Model');
AR19_eta50=load('NRes_A126_AR19_Eta_50_Model');

AR19_eta100_mass=AR19_eta100.Param.Sized_Masses.Total_mass;
AR19_eta90_mass=AR19_eta90.Param.Sized_Masses.Total_mass + 400;
AR19_eta80_mass=AR19_eta80.Param.Sized_Masses.Total_mass + 400;
AR19_eta70_mass=AR19_eta70.Param.Sized_Masses.Total_mass + 400;
AR19_eta60_mass=AR19_eta60.Param.Sized_Masses.Total_mass + 400;
AR19_eta50_mass=AR19_eta50.Param.Sized_Masses.Total_mass + 400;

AR_19_masses=[AR19_eta100_mass,AR19_eta90_mass,AR19_eta80_mass,AR19_eta70_mass,AR19_eta60_mass]*2 + 1400*2;

Eta=[0,10,20,30,40];

% baseline 
baseline321=7400*ones(1,5);

figure 

plot(Eta,AR_19_masses,'r-s','MarkerFaceColor','r','MarkerSize',6,'LineWidth',1.2)
hold on 

plot(Eta,baseline321,'k--','MarkerFaceColor','r','MarkerSize',6,'LineWidth',1.2)

xlabel('Size of folding wingtip, $\eta$ ($\%$)','Interpreter','latex','FontSize',12)
ylabel('Overall wing weight (kg)','Interpreter','latex','FontSize',12)

set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Color','w')
grid(gca,'minor')



%% Sizing loads
% plots to show the load alleviation of the FWT device.

AR19_eta100=load('NRes_A126_AR19_Eta_100_Model');
AR19_eta90=load('NRes_A126_AR19_Eta_90_Model');
AR19_eta80=load('NRes_A126_AR19_Eta_80_Model');
AR19_eta70=load('NRes_A126_AR19_Eta_70_Model');
AR19_eta60=load('NRes_A126_AR19_Eta_60_Model');
AR19_eta50=load('NRes_A126_AR19_Eta_50_Model');

AR19_eta100_Param=AR19_eta100.Param;
AR19_eta90_Param=AR19_eta90.Param;
AR19_eta80_Param=AR19_eta80.Param;
AR19_eta70_Param=AR19_eta70.Param;
AR19_eta60_Param=AR19_eta60.Param;
AR19_eta50_Param=AR19_eta50.Param;

[Load_Distribution_eta100, Sizing_Loads_eta100, Box_dimensions_eta100, Box_CrossSec_eta100]=Sizing_Evelope(AR19_eta100_Param,run_folder);
[Load_Distribution_eta90, Sizing_Loads_eta90, Box_dimensions_eta90, Box_CrossSec_eta90]=Sizing_Evelope(AR19_eta90_Param,run_folder);
[Load_Distribution_eta80, Sizing_Loads_eta80, Box_dimensions_eta80, Box_CrossSec_eta80]=Sizing_Evelope(AR19_eta80_Param,run_folder);

[Load_Distribution_eta60, Sizing_Loads_eta60, Box_dimensions_eta60, Box_CrossSec_eta60]=Sizing_Evelope(AR19_eta60_Param,run_folder);
[Load_Distribution_eta50, Sizing_Loads_eta50, Box_dimensions_eta50, Box_CrossSec_eta50]=Sizing_Evelope(AR19_eta50_Param,run_folder);
[Load_Distribution_eta70, Sizing_Loads_eta70, Box_dimensions_eta70, Box_CrossSec_eta70]=Sizing_Evelope(AR19_eta70_Param,run_folder);

figure 
plot(Sizing_Loads_eta100.Y, Sizing_Loads_eta100.Moment_P2,'rs-','MarkerFaceColor','r','MarkerSize',6,'LineWidth',1.2)
hold on 
plot(Sizing_Loads_eta90.Y, Sizing_Loads_eta90.Moment_P2,'gs-','MarkerFaceColor','g','MarkerSize',6,'LineWidth',1.2)
hold on 
plot(Sizing_Loads_eta80.Y, Sizing_Loads_eta80.Moment_P2,'bs-','MarkerFaceColor','g','MarkerSize',6,'LineWidth',1.2)
hold on 
plot(Sizing_Loads_eta70.Y, Sizing_Loads_eta70.Moment_P2,'s-')
hold on 
plot(Sizing_Loads_eta60.Y, Sizing_Loads_eta60.Moment_P2,'v-')
% hold on 
% plot(Sizing_Loads_eta50.Y, Sizing_Loads_eta50.Moment_P2,'o-')



figure % bending load

plot(Load_Distribution_eta100.Y-2, Load_Distribution_eta100.Case1.Moment_P2,'rs-','MarkerFaceColor','r','MarkerSize',6,'LineWidth',0.8)

hold on

plot(Load_Distribution_eta90.Y-2, Load_Distribution_eta90.Case1.Moment_P2,'gs-','MarkerFaceColor','g','MarkerSize',6,'LineWidth',0.8)
% hold on 
% plot(Load_Distribution_eta90.Y, Load_Distribution_eta90.Case3.Moment_P2,'b-')

hold on

plot(Load_Distribution_eta80.Y-2, Load_Distribution_eta80.Case1.Moment_P2,'ks-','MarkerFaceColor','y','MarkerSize',6,'LineWidth',0.8)


hold on 

plot(Load_Distribution_eta70.Y-2, Load_Distribution_eta70.Case1.Moment_P2,'bs-','MarkerFaceColor','b','MarkerSize',6,'LineWidth',0.8)
hold on 
% plot(Load_Distribution_eta70.Y, Load_Distribution_eta70.Case3.Moment_P2,'v-')
% 
% hold on 

plot(Load_Distribution_eta60.Y-2, Load_Distribution_eta60.Case1.Moment_P2,'ks-','MarkerFaceColor','m','MarkerSize',6,'LineWidth',0.8)
hold on 
plot(Load_Distribution_eta80.Y-2, Load_Distribution_eta80.Case3.Moment_P2,'k-.','LineWidth',2.5)


xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Bending moment (Nm)','Interpreter','latex','FontSize',12)

set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Color','w')
grid(gca,'minor')
axis([0 23 0 8e6])

legend('Manoeuvre (hinge locked)','Manoeuvre (free hinge, $\eta$ = 10$\%$)',...
    'Manoeuvre (free hinge, $\eta$ = 20$\%$)','Manoeuvre (free hinge, $\eta$ = 30$\%$)',...
    'Manoeuvre (free hinge, $\eta$ = 40$\%$)','Cruise (hinge locked)','Interpreter','latex','FontSize',12)


% hold on 
% plot(Load_Distribution_eta60.Y, Load_Distribution_eta60.Case3.Moment_P2,'s-')
% 
% hold on

% plot(Sizing_Loads_eta100.Y, Sizing_Loads_eta100.Moment_P2,'s-')
% hold on 
% plot(Sizing_Loads_eta90.Y, Sizing_Loads_eta90.Moment_P2,'s-')
% hold on 
% plot(Sizing_Loads_eta80.Y, Sizing_Loads_eta80.Moment_P2,'s-')
% hold on 
% plot(Sizing_Loads_eta70.Y, Sizing_Loads_eta70.Moment_P2,'s-')
% hold on 
% plot(Sizing_Loads_eta60.Y, Sizing_Loads_eta60.Moment_P2,'v-')





figure % vertical shear

plot(Load_Distribution_eta100.Y-2, Load_Distribution_eta100.Case1.Shear_P2,'rs-','MarkerFaceColor','r','MarkerSize',6,'LineWidth',0.8)

hold on

plot(Load_Distribution_eta90.Y-2, Load_Distribution_eta90.Case1.Shear_P2,'gs-','MarkerFaceColor','g','MarkerSize',6,'LineWidth',0.8)

hold on

plot(Load_Distribution_eta80.Y-2, Load_Distribution_eta80.Case1.Shear_P2,'ks-','MarkerFaceColor','y','MarkerSize',6,'LineWidth',0.8)


hold on 

plot(Load_Distribution_eta70.Y-2, Load_Distribution_eta70.Case1.Shear_P2,'bs-','MarkerFaceColor','b','MarkerSize',6,'LineWidth',0.8)
hold on 
% plot(Load_Distribution_eta70.Y, Load_Distribution_eta70.Case3.Moment_P2,'v-')
% 
% hold on 

plot(Load_Distribution_eta60.Y-2, Load_Distribution_eta60.Case1.Shear_P2,'ks-','MarkerFaceColor','m','MarkerSize',6,'LineWidth',0.8)
hold on 
plot(Load_Distribution_eta80.Y-2, Load_Distribution_eta80.Case3.Shear_P2,'k-.','LineWidth',2.5)


xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Shear force (N)','Interpreter','latex','FontSize',12)

set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Color','w')
grid(gca,'minor')
axis([0 23 0 7e5])

legend('Manoeuvre (hinge locked)','Manoeuvre (free hinge, $\eta$ = 10$\%$)',...
    'Manoeuvre (free hinge, $\eta$ = 20$\%$)','Manoeuvre (free hinge, $\eta$ = 30$\%$)',...
    'Manoeuvre (free hinge, $\eta$ = 40$\%$)','Cruise (hinge locked)','Interpreter','latex','FontSize',12)





figure % torque

Load_Distribution_eta100_=Load_Distribution_eta100;
Load_Distribution_eta90_=Load_Distribution_eta90;
Load_Distribution_eta80_=Load_Distribution_eta80;
Load_Distribution_eta70_=Load_Distribution_eta70;
Load_Distribution_eta60_=Load_Distribution_eta60;


[a1,b1]=find(Load_Distribution_eta100_.Y<7);
[a2,b2]=find(Load_Distribution_eta90_.Y<7);
[a3,b3]=find(Load_Distribution_eta80_.Y<7);
[a4,b4]=find(Load_Distribution_eta70_.Y<7);
[a5,b5]=find(Load_Distribution_eta60_.Y<7);

Load_Distribution_eta100_.Case1.Torque(a1)=-Load_Distribution_eta100_.Case1.Torque(a1);
Load_Distribution_eta90_.Case1.Torque(a2)=-Load_Distribution_eta90_.Case1.Torque(a2);
Load_Distribution_eta80_.Case1.Torque(a3)=-Load_Distribution_eta80_.Case1.Torque(a3);
Load_Distribution_eta70_.Case1.Torque(a4)=-Load_Distribution_eta70_.Case1.Torque(a4);
Load_Distribution_eta60_.Case1.Torque(a5)=-Load_Distribution_eta60_.Case1.Torque(a5);

Load_Distribution_eta80_.Case3.Torque(a3)=-Load_Distribution_eta80_.Case3.Torque(a3);

plot(Load_Distribution_eta100.Y-2, Load_Distribution_eta100_.Case1.Torque,'rs-','MarkerFaceColor','r','MarkerSize',6,'LineWidth',0.8)

hold on

plot(Load_Distribution_eta90.Y-2, Load_Distribution_eta90_.Case1.Torque,'gs-','MarkerFaceColor','g','MarkerSize',6,'LineWidth',0.8)

hold on

plot(Load_Distribution_eta80.Y-2, Load_Distribution_eta80_.Case1.Torque,'ks-','MarkerFaceColor','y','MarkerSize',6,'LineWidth',0.8)


hold on 

plot(Load_Distribution_eta70.Y-2, Load_Distribution_eta70_.Case1.Torque,'bs-','MarkerFaceColor','b','MarkerSize',6,'LineWidth',0.8)
hold on 
% plot(Load_Distribution_eta70.Y, Load_Distribution_eta70.Case3.Moment_P2,'v-')
% 
% hold on 

plot(Load_Distribution_eta60.Y-2, Load_Distribution_eta60_.Case1.Torque,'ks-','MarkerFaceColor','m','MarkerSize',6,'LineWidth',0.8)
hold on 
plot(Load_Distribution_eta80.Y-2, Load_Distribution_eta80_.Case3.Torque,'k-.','LineWidth',2.5)


xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Torque (Nm)','Interpreter','latex','FontSize',12)

set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Color','w')
grid(gca,'minor')
axis([0 23 -2.5e5 2.5e5])

legend('Manoeuvre (hinge locked)','Manoeuvre (free hinge, $\eta$ = 10$\%$)',...
    'Manoeuvre (free hinge, $\eta$ = 20$\%$)','Manoeuvre (free hinge, $\eta$ = 30$\%$)',...
    'Manoeuvre (free hinge, $\eta$ = 40$\%$)','Cruise (hinge locked)','Interpreter','latex','FontSize',12)



 

%% plot results

Run_folder='C:\Git\A321_sizing\hg_codes\results\test_temp2';

% AR19_eta100_Param=AR19_eta100.Param;

AR19_eta70_Param=AR19_eta70.Param;


[~,~,~,~,~,~,TrimLoad_HF1,~,Box_dimensions_A, Box_CrossSec_A]=Static_Trim_v1(AR19_eta70_Param, Run_folder, 'Load_Factor',2.5,'File_Name','Pull_up36000ft','Hinge_Lock','off','Altitude',36000,'Mach_Num',0.78);

[~,~,~,~,~,~,TrimLoad_HL1,~,Box_dimensions_B, Box_CrossSec_B]=Static_Trim_v1(AR19_eta70_Param, Run_folder, 'Load_Factor',1,'File_Name','Cruise36000ft','Hinge_Lock','on','Altitude',36000,'Mach_Num',0.78);


% Run_folder_=strcat(run_folder,'\hinge_locked');
model = mni.import_matran(fullfile(Run_folder,'Pull_up36000ft.dat'));
model.draw

%extract the data
f06 = mni.result.f06(fullfile(Run_folder,'Pull_up36000ft.f06'));
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
model.update('Scale',0.6)

view([-90 0])
set(gcf,'Color','w')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')





Run_folder_=strcat(Run_folder,'\hinge_locked');

model = mni.import_matran(fullfile(Run_folder_,'Cruise36000ft.dat'));
model.draw

%extract the data
f06 = mni.result.f06(fullfile(Run_folder_,'Cruise36000ft.f06'));
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
model.update('Scale',0.6)

view([-90 0])
set(gcf,'Color','w')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')



%% range 

% MTOW
MTOW=93500;

% Mass without wing (kg)
Mass0=40000;

% Payload
Payload=25000*0.5;

% Fuel mass
Fuel_Mass = MTOW-Mass0-Payload-AR_19_masses;

% Fuel_Mass(Fuel_Mass>25000)=12500;

Fuel_Mass(Fuel_Mass>0)=12500;

Range_eta100=Breguet_Range(19, AR_19_masses(1), 126, Payload, Fuel_Mass(1));
Range_eta90=Breguet_Range(19, AR_19_masses(2), 126, Payload, Fuel_Mass(2));
Range_eta80=Breguet_Range(19, AR_19_masses(3), 126, Payload, Fuel_Mass(3));
Range_eta70=Breguet_Range(19, AR_19_masses(4), 126, Payload, Fuel_Mass(4));
Range_eta60=Breguet_Range(19, AR_19_masses(5), 126, Payload, Fuel_Mass(5));

Range=[Range_eta100,Range_eta90,Range_eta80,Range_eta70,Range_eta60];
eta=[0,10,20,30,40];

% Baseline model 
Range_321=Breguet_Range(10, 7400, 126, Payload, 12500);
Baseline=Range_321*ones(1,5);

figure 
plot(eta,Range,'b-s','MarkerFaceColor','b','MarkerSize',6,'LineWidth',1.2)

hold on 

plot(eta,Baseline,'k--','LineWidth',1.2)

xlabel('Size of folding wingtip, $\eta$ ($\%$)','Interpreter','latex','FontSize',12)
ylabel('Range (km)','Interpreter','latex','FontSize',12)

set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Color','w')
grid(gca,'minor')

%% Drag polar 
Cl=0:0.1:1;

AR1=10;
AR2=19;

Cdi1=Cl.^2/(0.4*AR1*pi);
Cdi2=Cl.^2/(0.4*AR2*pi);

Cd0=0.03;

CD1=Cd0+Cdi1;
CD2=Cd0+Cdi2;


figure 

plot(CD1,Cl,'k--','LineWidth',1.2)
hold on 
plot(CD2,Cl,'gs-','MarkerFaceColor','g','MarkerSize',6,'LineWidth',1.2)

xlabel('Drag coefficient','Interpreter','latex','FontSize',12)
ylabel('Lift coefficient','Interpreter','latex','FontSize',12)

set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Color','w')
grid(gca,'minor')

axis([0 0.11 0 1])
