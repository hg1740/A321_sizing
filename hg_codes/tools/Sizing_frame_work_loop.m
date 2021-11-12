
eta_all=0.5;

for i =1:length(eta_all)
    

% Top level control parameters --------------------------
Param=eval('A321_v2');

% Select fold length 
Param.FWT.Fold_eta=eta_all(i);

% Update Wing properties
Param.Wing.AR=19; % 10.172 for A321

Param.Wing.TotalArea=126; % 126m^2 for A321

% Res_Name=strcat('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A160\AR25\Res_AR25_Eta_70');
Res_Name=strcat('C:\Git\A321_sizing\hg_codes\Data\sizing_v2\AR19\NRes_A126_AR19_Eta_',num2str(eta_all(i)*100));

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

% Update Wing planform 

Param.Wing.Root_Chord=Geo_Wing.Root_Chord;

Param.Wing.TE_Sweep1=Geo_Wing.TE_Sweep1;

Param.Wing.TE_Sweep2=Geo_Wing.TE_Sweep2;

Param.Wing.Span=Geo_Wing.Span;

Param.Wing.Semi_Span=(Param.Wing.Span-4)/2;

Param.Wing.HalfArea=(Param.Wing.TotalArea-Param.Wing.Root_Chord*4)/2;


if isfield(Param,'FWT') 
    
    % update FWT planeform 
    
    Param.FWT.Root_Chord=Geo_FWT.Root_Chord;
    
    Param.FWT.Root_Height=Geo_FWT.Root_Height;
    
    Param.FWT.Tip_Chord=Geo_FWT.Tip_Chord;
    
    Param.FWT.Tip_Height=Geo_FWT.Tip_Height;
    
    
end


%% Model check

if isfield(Param,'FWT')
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v3(Param);
    
    
else
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v3(Param);
    
end

draw(FEM_full);

%% Run analysis

run_folder = 'C:\Git\A321_sizing\hg_codes\results\test_temp'; %[-], folder for exporting the NASTRAN model

%% Material properties : Al7075

Yield_strength = 5.2e8;

%% Saftry factor

Safty_Factor=1.5;

Indicator2=0;

while Indicator2>1.05 || Indicator2<0.95
    
    
    % Jig twist optimisation
    Param_Jig = Jig_Twist(Param,run_folder);
    
     % Sizing load computatipn
    [Load_Distribution, Sizing_Loads, Box_dimensions, Box_CrossSec]=Sizing_Evelope(Param_Jig,run_folder);
    
    Param_Jig.Y=Sizing_Loads.Y;
      
    Skin_thickness0=Param_Jig.Wing.Skin_String.Skin_Thickness;
    
    Indicator1=0;
    
    while Indicator1>1.05 || Indicator1<0.95
        
        % spar sizing
        Param_Spar_Update=Spar_Sizing(Param_Jig, Box_dimensions, Sizing_Loads, Safty_Factor, Yield_strength);
                                            
        % skin stringer panel sizing
        Param_SS_Update=SkinStringer_Sizing(Param_Spar_Update, Box_dimensions, Sizing_Loads, Safty_Factor, Yield_strength);
        

        discrepancy = Param_SS_Update.Wing.Skin_String.Skin_Thickness./Param_Spar_Update.Wing.Skin_String.Skin_Thickness;
        
        Indicator1=max(discrepancy);
        
        disp(['Indicator1 ',num2str(Indicator1)]);
        
        Param_Jig=Param_SS_Update;
        
        
    end
      
    
    discrepancy2 = Param_Jig.Wing.Skin_String.Skin_Thickness./Skin_thickness0;
    
    Indicator2=max(discrepancy2);
    
    
    disp(['Indicator2 ',num2str(Indicator2)]);
    
    Param=Param_Jig;
        
    
end


Mass = Mass_Cal(Param);

[Rib_Thickness,Rib_Weight]=Rib_Sizing(Param,Sizing_Loads,Safty_Factor);

Total_Mass=Mass+Rib_Weight;

    
Param.Rib_Thickness=Rib_Thickness;
    
Param.Sized_Masses.Rib_mass=Rib_Weight;

Param.Sized_Masses.SSS_mass=Mass;

Param.Sized_Masses.Total_mass=Total_Mass;        
        
        
%% Write Results Output Files

% Model data

save(strcat(Res_Name,'_Model.mat'), 'Param')
        
end


% %% temp plot
% 
% AR=[10,12,14,16,18,20];
% 
% Wing_Mass126=([2817,3209,3672,4026.3,4518.1,5002]+1200)*2;
% 
% Wing_Mass140=([3074.4,3548,3934,4463.3,4932.1,5388]+1400)*2;
% 
% figure
% 
% plot(AR, Wing_Mass126,'bs','MarkerFaceColor','b')
% hold on
% plot(AR, Wing_Mass140,'rs','MarkerFaceColor','r')
% 

% axis([10 20 0.5e4 1.5e4]);
% 
% xlabel('Aspect Ratio','Interpreter','latex','FontSize',14)
% ylabel('Overal wing weight (kg)','Interpreter','latex','FontSize',14)
% 
% set(gca,'FontSize',14)
% set(gca,'TickLabelInterpreter','latex')
% 
% set(gcf,'Color','w')
% grid(gca,'minor')
%    
% legend('Wing area = 122 m$^2$','Wing area = 140 m$^2$','Interpreter','latex','FontSize',14,'Location','southeast')        