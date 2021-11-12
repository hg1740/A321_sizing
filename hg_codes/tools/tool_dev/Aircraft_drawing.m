% Top level control parameters --------------------------

Param=eval('A321_v2');

% % Select fold length 
Param.FWT.Fold_eta=1;

% Update Wing properties
Param.Wing.AR=19; % 10.172 for A321

Param.Wing.TotalArea=126;

Param.Wing.AeroPanel_AR=1;

Param.Wing.AeroPanel_Number=12;


%---------------------------------------------------------

if Param.FWT.Fold_eta==1
    
    Param=rmfield(Param,'FWT');
    
end

% Generate wing planform properties 

if isfield(Param,'FWT') 
    [Geo_Wing, Geo_FWT]= Wing_Gen_V1(Param);
else
    [Geo_Wing]= Wing_Gen_V1(Param);
end

Param.Wing.Root_Chord=Geo_Wing.Root_Chord;

Param.Wing.TE_Sweep1=Geo_Wing.TE_Sweep1;

Param.Wing.TE_Sweep2=Geo_Wing.TE_Sweep2;

Param.Wing.Span=Geo_Wing.Span;

Param.Wing.Semi_Span=(Param.Wing.Span-4)/2;

Param.Wing.HalfArea=(Param.Wing.TotalArea-Param.Wing.Root_Chord*4)/2;


[Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v3(Param);
% 
% draw(Aircraft)
% 
draw(FEM_full)

Mid_chord=0.63685*Param.Wing.Root_Chord;
Tip_chord=0.2248*Param.Wing.Root_Chord;
    
origin=[0,2];

pt1=origin;
pt6=origin-[Param.Wing.Root_Chord,0];

pt3=pt1+ [-Param.Wing.Semi_Span*tan(deg2rad(27)), Param.Wing.Semi_Span];
pt2=pt1+(pt3-pt1)*0.27;
pt4=pt3-[Tip_chord,0];
pt5=pt2-[Mid_chord,0];

PT_all.AR19=[pt1;pt2;pt3;pt4;pt5;pt6];

Body_pt1=[0,20];
Body_pt2=[0,-25];
Body_pt=[Body_pt1;Body_pt2];

Body.Eta = [0;0.005268;0.010536;0.015805;0.021073;...
    0.026342;0.03161;0.036879;0.042147;0.047415;0.052684;...
    0.057952;0.063221;0.0684890;0.073758;0.079026;0.084294;0.089563;0.094831;0.1001;0.411022;...
    0.721944;0.736578;0.751213;0.765847;0.780482;0.795117;0.809751;0.824386;0.83902;0.853655;...
    0.868289;0.882924;0.897558;0.912193;0.926827;0.941462;0.956096;0.970731;0.985365;1]';

Body.Radius = [0.01;0.3844030;0.565081;0.707928;0.830682;0.940375;...
    1.04067;1.13377;1.22112;1.30374;1.38237;1.45758;1.52981;1.59941;1.66667;...
    1.73182;1.79508;1.8566;1.91653;1.975;2.11455;2.11455;2.11455;2.11455;2.11455;...
    2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;1.9;...
    1.75;1.6;1.4;1.2;1.0;0.01]';
    
    
figure 

plot(PT_all.AR10(:,2),PT_all.AR10(:,1)+27,'k-','LineWidth',1)

% hold on 
% plot(PT_all.AR13(:,2),PT_all.AR13(:,1)+27,'m-','LineWidth',1)
% 
% hold on 
% plot(PT_all.AR16(:,2),PT_all.AR16(:,1)+27,'b-','LineWidth',1)
% 
hold on 
plot(PT_all.AR19(:,2),PT_all.AR19(:,1)+27,'r-','LineWidth',1)

% hold on 
% plot(PT_all.AR22(:,2),PT_all.AR22(:,1)+27,'r-','LineWidth',1)

hold on 

plot(Body.Radius,(Body.Eta*45-27)+27,'k-','LineWidth',1)

% another half
plot(-PT_all.AR10(:,2),PT_all.AR10(:,1)+27,'k-','LineWidth',1)

hold on 

% plot(-PT_all.AR13(:,2),PT_all.AR13(:,1)+27,'m-','LineWidth',1)
% 
% hold on 
% plot(-PT_all.AR16(:,2),PT_all.AR16(:,1)+27,'b-','LineWidth',1)
% 
% hold on 
plot(-PT_all.AR19(:,2),PT_all.AR19(:,1)+27,'r-','LineWidth',1)

% hold on 
% plot(-PT_all.AR22(:,2),PT_all.AR22(:,1)+27,'r-','LineWidth',1)

hold on 

plot(-Body.Radius,(Body.Eta*45),'k-','LineWidth',1)

axis equal
axis([-30 30 0 45])

set(gcf,'Color','w')
xlabel('Y-axis (m)','Interpreter','latex','FontSize',14)
ylabel('X-axis (m)','Interpreter','latex','FontSize',14)
legend('Baseline model','AR = 19','Interpreter','latex','FontSize',16)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')


% figure 
% 
% plot(Body.Radius,Body.Eta*45)
% axis equal
% AR22_wing.Ptx1=18.3133;
% AR22_wing.Pty1=2;
% 
% AR22_wing.Ptx2=30.7075;
% AR22_wing.Pty2=26.3249;
% 
% AR22_wing.Ptx3=31.6554;
% AR22_wing.Pty3=26.3249;
% 
% AR22_wing.Ptx4=24.3451;
% AR22_wing.Pty4=8.56772;
% 
% AR22_wing.Ptx5=22.53;
% AR22_wing.Pty5=2;


