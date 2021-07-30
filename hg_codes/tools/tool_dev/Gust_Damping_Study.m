
% Damp0 0.002
% Damp1 0.02
% Damp2 0.04
% Damp3 0.06
% Damp4 0.08
% Damp5 0.1
% Damp6 0.12

Wing_Model=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR16\Res_AR16_Eta_70_Model.mat');

Param=Wing_Model.Param;

run_folder='C:\Git\A321_sizing\hg_codes\results\Gust_Damping';

Loads_D0=Gust_Analysis_v1(Param,run_folder,'File_Name','Gust_Damp0','March_Number',0.48,'Altitude',3000,'Hinge_Lock','off');


Loads_D1=Gust_Analysis_v1(Param,run_folder,'File_Name','Gust_Damp1','March_Number',0.48,'Altitude',3000,'Hinge_Lock','off');


Loads_D2=Gust_Analysis_v1(Param,run_folder,'File_Name','Gust_Damp2','March_Number',0.48,'Altitude',3000,'Hinge_Lock','off');


Loads_D3=Gust_Analysis_v1(Param,run_folder,'File_Name','Gust_Damp3','March_Number',0.48,'Altitude',3000,'Hinge_Lock','off');


Loads_D4=Gust_Analysis_v1(Param,run_folder,'File_Name','Gust_Damp4','March_Number',0.48,'Altitude',3000,'Hinge_Lock','off');


Loads_D5=Gust_Analysis_v1(Param,run_folder,'File_Name','Gust_Damp5','March_Number',0.48,'Altitude',3000,'Hinge_Lock','off');


Loads_D6=Gust_Analysis_v1(Param,run_folder,'File_Name','Gust_Damp6','March_Number',0.48,'Altitude',3000,'Hinge_Lock','off');





figure % incremental moment 

Gust_Length=linspace(18,214,7);

plot(Gust_Length,Loads_D0.Root_Delta.Max_Moment,'bs-','MarkerFaceColor','b')
hold on 
plot(Gust_Length,Loads_D0.Root_Delta.Min_Moment,'bs-','MarkerFaceColor','b')

hold on 

plot(Gust_Length,Loads_D1.Root_Delta.Max_Moment,'ks-','MarkerFaceColor','r')
hold on 
plot(Gust_Length,Loads_D1.Root_Delta.Min_Moment,'ks-','MarkerFaceColor','r')

hold on 

plot(Gust_Length,Loads_D2.Root_Delta.Max_Moment,'ks-','MarkerFaceColor','k')
hold on 
plot(Gust_Length,Loads_D2.Root_Delta.Min_Moment,'ks-','MarkerFaceColor','k')

hold on 

plot(Gust_Length,Loads_D3.Root_Delta.Max_Moment,'ks-','MarkerFaceColor','g')
hold on 
plot(Gust_Length,Loads_D3.Root_Delta.Min_Moment,'ks-','MarkerFaceColor','g')

hold on 

plot(Gust_Length,Loads_D4.Root_Delta.Max_Moment,'ks-','MarkerFaceColor','y')
hold on 
plot(Gust_Length,Loads_D4.Root_Delta.Min_Moment,'ks-','MarkerFaceColor','y')


hold on 

plot(Gust_Length,Loads_D5.Root_Delta.Max_Moment,'ks-','MarkerFaceColor','m')
hold on 
plot(Gust_Length,Loads_D5.Root_Delta.Min_Moment,'ks-','MarkerFaceColor','m')

set(gcf,'Color','w')

xlabel('Gust Length (m)','interpreter','latex','FontSize',14)

ylabel('Incremental bending moment (Nm)','interpreter','latex','FontSize',14)
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
axis([0 250 -1.5e6 1.5e6])

axP = get(gca,'Position'); 
set(gca, 'Position', axP)

legend('No damping','','Damping = 0.02','','Damping = 0.04','','Damping = 0.06','','Damping = 0.08','','Damping = 0.1','interpreter','latex','Location','eastoutside','FontSize',12);

set('DataAspectRatio',[2,1,1]); 






figure % incremental shear

Gust_Length=linspace(18,214,7);

plot(Gust_Length,Loads_D0.Root_Delta.Max_Shear,'bs-','MarkerFaceColor','b')
hold on 
plot(Gust_Length,Loads_D0.Root_Delta.Min_Shear,'bs-','MarkerFaceColor','b')

hold on 

plot(Gust_Length,Loads_D1.Root_Delta.Max_Shear,'ks-','MarkerFaceColor','r')
hold on 
plot(Gust_Length,Loads_D1.Root_Delta.Min_Shear,'ks-','MarkerFaceColor','r')

hold on 

plot(Gust_Length,Loads_D2.Root_Delta.Max_Shear,'ks-','MarkerFaceColor','k')
hold on 
plot(Gust_Length,Loads_D2.Root_Delta.Min_Shear,'ks-','MarkerFaceColor','k')

hold on 

plot(Gust_Length,Loads_D3.Root_Delta.Max_Shear,'ks-','MarkerFaceColor','g')
hold on 
plot(Gust_Length,Loads_D3.Root_Delta.Min_Shear,'ks-','MarkerFaceColor','g')

hold on 

plot(Gust_Length,Loads_D4.Root_Delta.Max_Shear,'ks-','MarkerFaceColor','y')
hold on 
plot(Gust_Length,Loads_D4.Root_Delta.Min_Shear,'ks-','MarkerFaceColor','y')


hold on 

plot(Gust_Length,Loads_D5.Root_Delta.Max_Shear,'ks-','MarkerFaceColor','m')
hold on 
plot(Gust_Length,Loads_D5.Root_Delta.Min_Shear,'ks-','MarkerFaceColor','m')

set(gcf,'Color','w')

xlabel('Gust Length (m)','interpreter','latex','FontSize',14)

ylabel('Incremental shear force (N)','interpreter','latex','FontSize',14)
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
% axis([0 250 -1.5e6 1.5e6])

axP = get(gca,'Position'); 
set(gca, 'Position', axP)

legend('No damping','','Damping = 0.02','','Damping = 0.04','','Damping = 0.06','','Damping = 0.08','','Damping = 0.1','interpreter','latex','Location','eastoutside','FontSize',12);




figure % incremental torque

Gust_Length=linspace(18,214,7);

plot(Gust_Length,Loads_D0.Root_Delta.Max_Torque,'bs-','MarkerFaceColor','b')
hold on 
plot(Gust_Length,Loads_D0.Root_Delta.Min_Torque,'bs-','MarkerFaceColor','b')

hold on 

plot(Gust_Length,Loads_D1.Root_Delta.Max_Torque,'ks-','MarkerFaceColor','r')
hold on 
plot(Gust_Length,Loads_D1.Root_Delta.Min_Torque,'ks-','MarkerFaceColor','r')

hold on 

plot(Gust_Length,Loads_D2.Root_Delta.Max_Torque,'ks-','MarkerFaceColor','k')
hold on 
plot(Gust_Length,Loads_D2.Root_Delta.Min_Torque,'ks-','MarkerFaceColor','k')

hold on 

plot(Gust_Length,Loads_D3.Root_Delta.Max_Torque,'ks-','MarkerFaceColor','g')
hold on 
plot(Gust_Length,Loads_D3.Root_Delta.Min_Torque,'ks-','MarkerFaceColor','g')

hold on 

plot(Gust_Length,Loads_D4.Root_Delta.Max_Torque,'ks-','MarkerFaceColor','y')
hold on 
plot(Gust_Length,Loads_D4.Root_Delta.Min_Torque,'ks-','MarkerFaceColor','y')


hold on 

plot(Gust_Length,Loads_D5.Root_Delta.Max_Torque,'ks-','MarkerFaceColor','m')
hold on 
plot(Gust_Length,Loads_D5.Root_Delta.Min_Torque,'ks-','MarkerFaceColor','m')

set(gcf,'Color','w')

xlabel('Gust Length (m)','interpreter','latex','FontSize',14)

ylabel('Incremental Torque (Nm)','interpreter','latex','FontSize',14)
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
% axis([0 250 -1.5e6 1.5e6])

axP = get(gca,'Position'); 
set(gca, 'Position', axP)

legend('No damping','','Damping = 0.02','','Damping = 0.04','','Damping = 0.06','','Damping = 0.08','','Damping = 0.1','interpreter','latex','Location','eastoutside','FontSize',12);








figure 

plot(Loads_D0.Y(1:end-1)-2,Loads_D0.Wing_Delta.Max_Moment,'bs-','MarkerFaceColor','b')
hold on 
plot(Loads_D0.Y(1:end-1)-2,Loads_D0.Wing_Delta.Min_Moment,'bs-','MarkerFaceColor','b')

hold on 

plot(Loads_D0.Y(1:end-1)-2,Loads_D2.Wing_Delta.Max_Moment,'rs-','MarkerFaceColor','r')
hold on 
plot(Loads_D0.Y(1:end-1)-2,Loads_D2.Wing_Delta.Min_Moment,'rs-','MarkerFaceColor','r')

hold on 

plot(Loads_D1.Y(1:end-1)-2,Loads_D4.Wing_Delta.Max_Moment,'ks-','MarkerFaceColor','y')
hold on 
plot(Loads_D1.Y(1:end-1)-2,Loads_D4.Wing_Delta.Min_Moment,'ks-','MarkerFaceColor','y')


hold on 

plot(Loads_D1.Y(1:end-1)-2,Loads_D6.Wing_Delta.Max_Moment,'ks-','MarkerFaceColor','g')
hold on 
plot(Loads_D1.Y(1:end-1)-2,Loads_D6.Wing_Delta.Min_Moment,'ks-','MarkerFaceColor','g')



set(gcf,'Color','w')

xlabel('Spanwise distance (m)','interpreter','latex','FontSize',14)

ylabel('Incremental bending moment (Nm)','interpreter','latex','FontSize',14)
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
% axis([0 250 -1.5e6 1.5e6])


legend('No damping','','Damping = 0.04','','Damping = 0.08','','Damping = 0.12','interpreter','latex','Location','northeast','FontSize',12);





figure 

plot(Loads_D0.Y(1:end-1)-2,Loads_D0.Wing_Delta.Max_Shear,'bs-','MarkerFaceColor','b')
hold on 
plot(Loads_D0.Y(1:end-1)-2,Loads_D0.Wing_Delta.Min_Shear,'bs-','MarkerFaceColor','b')

hold on 

plot(Loads_D0.Y(1:end-1)-2,Loads_D2.Wing_Delta.Max_Shear,'rs-','MarkerFaceColor','r')
hold on 
plot(Loads_D0.Y(1:end-1)-2,Loads_D2.Wing_Delta.Min_Shear,'rs-','MarkerFaceColor','r')

hold on 

plot(Loads_D1.Y(1:end-1)-2,Loads_D4.Wing_Delta.Max_Shear,'ks-','MarkerFaceColor','y')
hold on 
plot(Loads_D1.Y(1:end-1)-2,Loads_D4.Wing_Delta.Min_Shear,'ks-','MarkerFaceColor','y')


hold on 

plot(Loads_D1.Y(1:end-1)-2,Loads_D6.Wing_Delta.Max_Shear,'ks-','MarkerFaceColor','g')
hold on 
plot(Loads_D1.Y(1:end-1)-2,Loads_D6.Wing_Delta.Min_Shear,'ks-','MarkerFaceColor','g')



set(gcf,'Color','w')

xlabel('Spanwise distance (m)','interpreter','latex','FontSize',14)

ylabel('Incremental shear force (N)','interpreter','latex','FontSize',14)
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
% axis([0 250 -1.5e6 1.5e6])


legend('No damping','','Damping = 0.04','','Damping = 0.08','','Damping = 0.12','interpreter','latex','Location','northeast','FontSize',12);







figure 

plot(Loads_D0.Y(1:end-1)-2,Loads_D0.Wing_Delta.Max_Torque,'bs-','MarkerFaceColor','b')
hold on 
plot(Loads_D0.Y(1:end-1)-2,Loads_D0.Wing_Delta.Min_Torque,'bs-','MarkerFaceColor','b')

hold on 

plot(Loads_D0.Y(1:end-1)-2,Loads_D2.Wing_Delta.Max_Torque,'rs-','MarkerFaceColor','r')
hold on 
plot(Loads_D0.Y(1:end-1)-2,Loads_D2.Wing_Delta.Min_Torque,'rs-','MarkerFaceColor','r')

hold on 

plot(Loads_D1.Y(1:end-1)-2,Loads_D4.Wing_Delta.Max_Torque,'ks-','MarkerFaceColor','y')
hold on 
plot(Loads_D1.Y(1:end-1)-2,Loads_D4.Wing_Delta.Min_Torque,'ks-','MarkerFaceColor','y')


hold on 

plot(Loads_D1.Y(1:end-1)-2,Loads_D6.Wing_Delta.Max_Torque,'ks-','MarkerFaceColor','g')
hold on 
plot(Loads_D1.Y(1:end-1)-2,Loads_D6.Wing_Delta.Min_Torque,'ks-','MarkerFaceColor','g')



set(gcf,'Color','w')

xlabel('Spanwise distance (m)','interpreter','latex','FontSize',14)

ylabel('Incremental shear force (N)','interpreter','latex','FontSize',14)
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
% axis([0 250 -1.5e6 1.5e6])


legend('No damping','','Damping = 0.04','','Damping = 0.08','','Damping = 0.12','interpreter','latex','Location','northeast','FontSize',12);



figure 
T=linspace(0,2.5,201);
plot(T,Loads_D0.Time_Response.Root_Moment,'k-')
hold on 
% plot(T,Loads_D2.Time_Response.Root_Moment,'r.')
% hold on 
% plot(T,Loads_D4.Time_Response.Root_Moment,'g.')
% hold on 
plot(T,Loads_D6.Time_Response.Root_Moment,'b-')






%% load distributions


AR16_Eta70_Damp0=load('Sizing_with_Damping\Res_AR16_Eta_70_Damp0_Loads');
AR16_Eta70_Damp012=load('Sizing_with_Damping\Res_AR16_Eta_70_Damp012_Loads');



figure % bending moment

plot(AR16_Eta70_Damp0.Loads.Y-2, AR16_Eta70_Damp0.Loads.Moment.g_sealevel_gust1,'r-','LineWidth',1)
 
hold on

plot(AR16_Eta70_Damp0.Loads.Y-2, AR16_Eta70_Damp0.Loads.Moment.hinge_lock,'r--','LineWidth',1)

hold on 

plot(AR16_Eta70_Damp012.Loads.Y-2, AR16_Eta70_Damp012.Loads.Moment.g_sealevel_gust1,'b-','LineWidth',1)
 
hold on

plot(AR16_Eta70_Damp012.Loads.Y-2, AR16_Eta70_Damp012.Loads.Moment.hinge_lock,'b--','LineWidth',1)



axis([0 20.5 0 4e6])
set(gcf,'Color','w')
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Bending moment (Nm)','Interpreter','latex','FontSize',12)

set(gca,'TickLabelInterpreter','latex')
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)

legend('g + 1MC (No damping)', 'Hinge locked case (No damping)','g + 1MC (Damping = 0.12)','Hinge locked case (Damping = 0.12)','Sizing','Interpreter','latex','FontSize',12)





figure % shear forces

plot(AR16_Eta70_Damp0.Loads.Y-2, AR16_Eta70_Damp0.Loads.Shear.g_sealevel_gust1,'r-','LineWidth',1)
 
hold on

plot(AR16_Eta70_Damp0.Loads.Y-2, AR16_Eta70_Damp0.Loads.Shear.hinge_lock,'r--','LineWidth',1)

hold on 

plot(AR16_Eta70_Damp012.Loads.Y-2, AR16_Eta70_Damp012.Loads.Shear.g_sealevel_gust1,'b-','LineWidth',1)
 
hold on

plot(AR16_Eta70_Damp012.Loads.Y-2, AR16_Eta70_Damp012.Loads.Shear.hinge_lock,'b--','LineWidth',1)



% axis([0 20.5 0 4e6])
set(gcf,'Color','w')
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Shear force (N)','Interpreter','latex','FontSize',12)

set(gca,'TickLabelInterpreter','latex')

legend('g + 1MC (No damping)', 'Hinge locked case (No damping)','g + 1MC (Damping = 0.12)','Hinge locked case (Damping = 0.12)','Sizing','Interpreter','latex','FontSize',12)





figure % shear forces

plot(AR16_Eta70_Damp0.Loads.Y-2, abs(AR16_Eta70_Damp0.Loads.Torque.g_sealevel_gust1),'r-','LineWidth',1)
 
hold on

plot(AR16_Eta70_Damp0.Loads.Y-2, abs(AR16_Eta70_Damp0.Loads.Torque.hinge_lock),'r--','LineWidth',1)

hold on 

plot(AR16_Eta70_Damp012.Loads.Y-2, abs(AR16_Eta70_Damp012.Loads.Torque.g_sealevel_gust1),'b-','LineWidth',1)
 
hold on

plot(AR16_Eta70_Damp012.Loads.Y-2, abs(AR16_Eta70_Damp012.Loads.Torque.hinge_lock),'b--','LineWidth',1)



% axis([0 20.5 0 4e6])
set(gcf,'Color','w')
grid on
grid minor
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
ylabel('Torque (Nm)','Interpreter','latex','FontSize',12)

set(gca,'TickLabelInterpreter','latex')

legend('g + 1MC (No damping)', 'Hinge locked case (No damping)','g + 1MC (Damping = 0.12)','Hinge locked case (Damping = 0.12)','Sizing','Interpreter','latex','FontSize',12)
