% Top level control parameters --------------------------
Param=eval('A321_v1');

% Select fold length 
Param.FWT.Fold_eta=1;

% Update Wing properties
Param.Wing.AR=10.172; % 10.172 for A321

Param.Wing.TotalArea=126; % 126m^2 for A321

% Res_Name=strcat('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A160\AR25\Res_AR25_Eta_70');
Res_Name=strcat('C:\Git\A321_sizing\hg_codes\Data\sizing_correction\A321\Res_AR10_Eta_',num2str(100));

% Update Initial Guess
% prevous_result=load('D22_results\Res_AR22_Eta_100_Model');

% Param.Wing.Thickness=prevous_result.Param.Wing.Thickness;

% if isfield(prevous_result.Param,'FWT')
%     
%     Param.FWT.Thickness=prevous_result.Param.FWT.Thickness;
%     
% end

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


% %% Model check
% 
% if isfield(Param,'FWT')
%     
%     [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v2(Param);
%     
%     
% else
%     
%     [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v2(Param);
%     
% end
% 
% draw(FEM_full);

%% Run analysis

run_folder = 'C:\Git\A321_sizing\hg_codes\results\test_temp'; %[-], folder for exporting the NASTRAN model

%% Material properties : Al7075

Yield_strength = 5.2e8;

%% Saftry factor

Safty_Factor=1.5;

counter=1;

%initial condition
Constraint_Max=1.1;
Constraint_Min=0.9;

if isfield(Param,'FWT')
    
    Numsec=35;
    record=zeros(100,140);
    
else
    
    Numsec=25;
    record=zeros(100,100);
    
end


while Constraint_Max>1.05 || Constraint_Min<0.95
    
    
    record(counter,:)=[Param.Wing.SparCap_Thickness, Param.Wing.SparWeb_Thickness, Param.Wing.Skin_Thickness, Param.Wing.Stringer_Area];
    
    
    % Jig twist optimisation
    Param_Jig = Jig_Twist(Param,run_folder);
    
    
    % Sizing load computatipn
    [Load_Distribution, Sizing_Loads, Box_dimensions, Box_CrossSec]=Sizing_Evelope(Param_Jig,run_folder);
    
    
    % Internal stresses computation
    Internal_Stresses=Internal_Stress_Calc(Param_Jig, Box_CrossSec, Box_dimensions, Sizing_Loads);
    
    
    % Constraint checks
    [Constraints, Constraint_uppers, Constraint_Max, Constraint_Min]  = Constraint_check(Internal_Stresses, Safty_Factor);
    

    % Thickness update
    Param_Thickness_Update=Thickness_Update(Param_Jig,Constraint_uppers,Internal_Stresses, Safty_Factor, Yield_strength);
    


    Param = Param_Thickness_Update;
    
    counter=counter+1;
    
    disp(['Max_constraint ',num2str(Constraint_Max)]);
    
    disp(['Min_constraint ',num2str(Constraint_Min)]);
    
    
end
        
        
        
%% Write Results Output Files

% Model data

save(strcat(Res_Name,'_Model.mat'), 'Param')
        
        
        
        
        
        