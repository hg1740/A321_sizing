%% A321 model stretched 

% Top level control parameters --------------------------

    
Param=eval('A321_v2');

% Select fold length 
Param.FWT.Fold_eta=1;

% Update Wing properties
Param.Wing.AR=14; % 10.172 for A321

Param.Wing.TotalArea=112; % 126m^2 for A321


% Update Mass properties
Param.Masses.Payload=20000; % 

% Res_Name=strcat('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A160\AR25\Res_AR25_Eta_70');
Res_Name=strcat('C:\Git\A321_sizing\hg_codes\Data\concept\Res_concept_',num2str(100));


%---------------------------------------------------------

if Param.FWT.Fold_eta==1
    
    Param=rmfield(Param,'FWT');
    
end

% Generate wing planform properties 

if isfield(Param,'FWT') 
    [Geo_Wing, Geo_FWT]= Wing_Gen_temp(Param);
else
    [Geo_Wing]= Wing_Gen_temp(Param);
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
    
    
%     % Jig twist optimisation
%     Param_Jig = Jig_Twist(Param_Sizing,run_folder);
    
     % Sizing load computatipn
    [Load_Distribution, Sizing_Loads, Box_dimensions, Box_CrossSec]=Sizing_Evelope(Param,run_folder);
      
    Skin_thickness0=Param.Wing.Skin_String.Skin_Thickness;
    
    Indicator1=0;
    
    while Indicator1>1.05 || Indicator1<0.95
        
        % spar sizing
        Param_Spar_Update=Spar_Update(Param, Box_dimensions, Sizing_Loads, Safty_Factor, Yield_strength);
                                            
        % skin stringer panel sizing
        Param_SS_Update=SkinStringer_Sizing(Param_Spar_Update, Box_dimensions, Sizing_Loads, Safty_Factor, Yield_strength);
        

        
        discrepancy = Param_SS_Update.Wing.Skin_String.Skin_Thickness./Param_Spar_Update.Wing.Skin_String.Skin_Thickness;
        
        Indicator1=max(discrepancy);
        
        disp(['Indicator1 ',num2str(Indicator1)]);
        
        Param=Param_SS_Update;
        
        
    end
    
    
    
    discrepancy2 = Param.Wing.Skin_String.Skin_Thickness./Skin_thickness0;
    
    Indicator2=max(discrepancy2);
    
    
    disp(['Indicator2 ',num2str(Indicator2)]);
        
    
end



Mass = Mass_Cal(Param);

[Rib_Thickness,Rib_Weight]=Rib_Sizing(Param,Sizing_Loads,Safty_Factor);

Total_Mass=Mass+Rib_Weight;

Param.Wing.Total_mass=Total_Mass;

set(gcf,'Color','w')



