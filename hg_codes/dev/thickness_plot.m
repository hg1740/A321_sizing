AR10_results=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR10\Res_AR10_Eta_100_Structure');
AR13_results=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR13\Res_AR13_Eta_100_Structure');
AR16_results=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR16\Res_AR16_Eta_100_Structure');
AR19_results=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR19\Res_AR19_Eta_100_Structure');


%Y_data
AR10_Y=AR10_results.Wingbox_property.NodesY;
AR13_Y=AR13_results.Wingbox_property.NodesY;
AR16_Y=AR16_results.Wingbox_property.NodesY;
AR19_Y=AR19_results.Wingbox_property.NodesY;

% Spar thickness
AR10_Spar=AR10_results.Wingbox_property.Spar_Thickness;
AR13_Spar=AR13_results.Wingbox_property.Spar_Thickness;
AR16_Spar=AR16_results.Wingbox_property.Spar_Thickness;
AR19_Spar=AR19_results.Wingbox_property.Spar_Thickness;

% Skin thickness
AR10_Skin=AR10_results.Wingbox_property.Skin_Thickness;
AR13_Skin=AR13_results.Wingbox_property.Skin_Thickness;
AR16_Skin=AR16_results.Wingbox_property.Skin_Thickness;
AR19_Skin=AR19_results.Wingbox_property.Skin_Thickness;

figure 

plot(AR10_Y-2,AR10_Spar*1000,'bs','MarkerFaceColor','b')
hold on 
plot(AR13_Y-2,AR13_Spar*1000,'rs','MarkerFaceColor','r')
hold on 
plot(AR16_Y-2,AR16_Spar*1000,'ks','MarkerFaceColor','g')
% hold on 
% plot(AR16_Y-2,AR19_Spar*1000,'ks','MarkerFaceColor','y')

xlabel('Spanwise distance','interpreter','latex','FontSize',12)

ylabel('Spar thickness (mm)','interpreter','latex','FontSize',12)

set(gcf,'Color','w')
legend('AR=10','AR=13','AR=16','AR=19','interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)

axis([0 20.5 0 70])


figure 

plot(AR10_Y-2,AR10_Skin*1000,'bs','MarkerFaceColor','b')
hold on 
plot(AR13_Y-2,AR13_Skin*1000,'rs','MarkerFaceColor','r')
hold on 
plot(AR16_Y-2,AR16_Skin*1000,'ks','MarkerFaceColor','g')
% hold on 
% plot(AR19_Y-2,AR19_Skin*1000,'ks','MarkerFaceColor','y')

xlabel('Spanwise distance','interpreter','latex','FontSize',12)

ylabel('Skin thickness (mm)','interpreter','latex','FontSize',12)

set(gcf,'Color','w')
legend('AR=10','AR=13','AR=16','AR=19','interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)

axis([0 20.5 0 20])