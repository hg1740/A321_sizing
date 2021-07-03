
% Aircraft Model selection 

Param=eval('Hodges');


% FWT update

Param.FWT.Fold_angle=10;

Param.FWT.Flare_angle=25;

Param.FWT.Fold_eta=0.8;

Param.FWT.Root_Chord=5;

Param.FWT.Root_Height=5*0.15;

Param.FWT.Tip_Chord=5;

Param.FWT.Tip_Height=5*0.15;

Param.FWT.Thickness=0.002;

Param.FWT.Hinge_Stiffness=1e-4;


% Span update

Param.Wing.Span=(34-4)/Param.FWT.Fold_eta + 4;

Param.Wing.Semi_Span=(Param.Wing.Span-4)/2;

% Thickness guess update

Param.Wing.Thickness=[0.0223433000000000,0.0221239000000000,0.0221090520000000,0.0223257510000000,0.0221091280000000,0.0222707790000000,0.0221770220000000,0.0225898920000000,0.0229393770000000,0.0222648460000000,0.0226110480000000,0.0218236060000000,0.0220503620000000,0.0210794660000000,0.0211146280000000,0.0199168260000000,0.0203244920000000,0.0189565460000000,0.0185290220000000,0.0167415890000000,0.0158305210000000,0.0140566400000000,0.0112756260000000,0.00856440100000000,0.00795457000000000,0.00862650400000000,0.00844402300000000,0.00825119400000000,0.00804985600000000,0.00784702700000000,0.00763136700000000,0.00741284200000000,0.00717615100000000,0.00692504600000000,0.00667816000000000,0.00641480600000000,0.00614569300000000,0.00586020500000000,0.00556658800000000,0.00526062100000000,0.00494440600000000,0.00460896900000000,0.00428528500000000,0.00397395200000000,0.00367685200000000,0.00332539500000000,0.00304101800000000,0.00293599900000000,0.00287679900000000,0.00269164800000000,9.51000000000000e-05,8.95000000000000e-05,8.43000000000000e-05,7.96000000000000e-05,7.54000000000000e-05,7.14000000000000e-05,6.73000000000000e-05,6.30000000000000e-05,5.86000000000000e-05,5.46000000000000e-05,5.03000000000000e-05,4.62000000000000e-05,4.20000000000000e-05,3.80000000000000e-05,3.39000000000000e-05,3.01000000000000e-05,2.60000000000000e-05,2.25000000000000e-05,1.92000000000000e-05,1.64000000000000e-05,1.32000000000000e-05,1.05000000000000e-05,9.63000000000000e-06,8.75000000000000e-06,7.05000000000000e-06];

% Model check 
[Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models(Param);

draw(Aircraft)

%% Run folder 

run_folder = ['C:\Git\A321_sizing\hg_codes\results\Hodges_ETA20']; %[-], folder for exporting the NASTRAN model

%% Saftry Factor

SF=1;

%% Yield strength 

Yield_strength=5.2e8;


%% Sizing loop

cond_set1=[1.1,1.1,1.1,1.1];
cond_set2=[1.1,1.1,1.1,1.1];
counter=1;
record=zeros(100,75);

while max(cond_set1)>1.05 || min(cond_set2)<0.95
    
    record(counter,:)=Param.Wing.Thickness;
    
    
    [Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis(Param, run_folder);
    
    % collect internal stresses
    RVon_skn=Internal_Stresses.VonMise_skin;
    RVon_spr=Internal_Stresses.VonMise_spar;
    Rsigmab_skn=Internal_Stresses.Skin_BucklingCritical;
    Rtaub_skn=Internal_Stresses.Skin_buckling_shear;
    Rsigma_pp=Internal_Stresses.Skin_buckling_principal;
    Rsigma_strg=Internal_Stresses.Stringer_compressive;
    Rsigma_crip=Internal_Stresses.Stringer_crippling;
    Rsigma_col=Internal_Stresses.Stringer_columnbuckling;
    
    % skin check
    Check_von_skin1=SF*RVon_skn/Yield_strength;
    Check_von_skin2=SF*Rsigma_pp./Rsigmab_skn;
    
    V_skin_change=step_size(Check_von_skin1);
    B_skin_change=step_size(Check_von_skin2);
    
    skin_coefficient=max([V_skin_change;B_skin_change]);
    
    Param.Wing.Thickness(26:50)=Param.Wing.Thickness(26:50).*skin_coefficient;
    
    % spar check
    Check_von_spar=SF*RVon_spr/Yield_strength;
    
    Spar_coefficient=step_size(Check_von_spar);
    
    Param.Wing.Thickness(1:25)=Param.Wing.Thickness(1:25).*Spar_coefficient;
    
    % stringers check
    Check_strg=SF*Rsigma_strg./Rsigma_col;
    
    Stringer_coefficient=step_size(Check_strg);
    
    Param.Wing.Thickness(51:75)=Param.Wing.Thickness(51:75).*Stringer_coefficient;
    
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

% Calculate Wingf mass
[Wingbox_mass, Secondary_mass, Wing_total_mass]=Wing_masses_v1(Param.Wing.Thickness,Box_dimensions.Inboard.Height,Box_dimensions.Inboard.Width, Y_eta, Param.Masses.MTOW, Param.Wing.HalfArea, Param.Wing.LE_Sweep, Param.Wing.Semi_Span);
        
% Drag calculation 
[CDi,CD0,CL,k,Distribution]=Drag_Calculation(Param, run_folder);


% Y_all=Y_all';
% 
% % AMoment_P2.LC1=Load_distribution.Moment_P2.LC1';
% % AShear_P2.LC1=Load_distribution.Shear_P2.LC1';
% % ATorque.LC1=Load_distribution.Torque.LC1';
% 
% 
% 
figure
plot(Y_all, Load_distribution.Moment_P2.LC1,'s-')

figure
plot(Y_all, Load_distribution.Shear_P2.LC1,'s-')


figure
plot(Y_all(1:end-1), Wing_Delta.Moment_Max.LC1,'s-')
% 
% 
figure
gust_length=linspace(18,214,7);
plot(gust_length, Root_Delta.Moment_Max.LC1,'s-')

hold on 

plot(gust_length, Root_Delta.Moment_Min.LC1,'s-')
% 
% 
% figure 
% plot(Distribution.Y,Distribution.Lift_var,'s-')
% 
% figure 
% plot(Distribution.Y,Distribution.WJ,'s-')


%% Long example 

Param=eval('Hodges');

% Update thickness

Param.Wing.Thickness=[0.01*ones(1,25),0.003*ones(1,25),1e-10*ones(1,25)];
Param.Wing.Root_Chord=6;

% Update masses

Param.Masses.Engine=0.01*7362/2;

Param.Masses.Pylon=0.01*1239/2;

Param.Masses.Secondary_Mass=0.01*1000;

Param.Masses.Fuel_Fraction=0.723*0.01;

Param.Masses.Fuel_Density=840; %g/L

Param.Masses.Fuel_Mass=Param.Masses.Fuel_Capacity * Param.Masses.Fuel_Fraction * Param.Masses.Fuel_Density/1000; 

run_folder = ['C:\Git\A321_sizing\hg_codes\results\Long_compare']; %[-], folder for exporting the NASTRAN model

% Drag calculation 
[CDi,CD0,CL,k,Distribution]=Drag_Calculation(Param, run_folder);

figure 
plot(Distribution.Y,Distribution.Alphai,'b.')


figure 
plot(Distribution.Y,Distribution.Lift_var,'b.')


model = mni.import_matran(fullfile(run_folder,'A321_36000ft_1g_drag.dat'));
model.draw

%extract the data
f06 = mni.result.f06(fullfile(run_folder,'A321_36000ft_1g_drag.f06'));
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


