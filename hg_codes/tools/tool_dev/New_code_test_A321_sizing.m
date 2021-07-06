%A321

x=[0.0191135910000000,0.0193953730000000,0.0201137600000000,0.0213823650000000,0.0229999860000000,0.0255251240000000,0.0285624380000000,0.0372113470000000,0.0373207110000000,0.0372799900000000,0.0368848460000000,0.0378039760000000,0.0379687300000000,0.0379303110000000,0.0378464170000000,0.0368959070000000,0.0365928770000000,0.0358741000000000,0.0351157990000000,0.0325204110000000,0.0321366360000000,0.0317223020000000,0.0280889430000000,0.0226448800000000,0.0118821390000000,0.00707914700000000,0.00718361400000000,0.00730636500000000,0.00745189900000000,0.00763438400000000,0.00783787200000000,0.00810313600000000,0.00823100800000000,0.00803824900000000,0.00781619300000000,0.00760075000000000,0.00736034800000000,0.00709415700000000,0.00680991600000000,0.00650766900000000,0.00618441800000000,0.00597655500000000,0.00574729800000000,0.00545697300000000,0.00516124700000000,0.00473319700000000,0.00418339100000000,0.00343190400000000,0.00256919700000000,0.000896999000000000,6.20000000000000e-05,6.38000000000000e-05,6.58000000000000e-05,6.82000000000000e-05,7.14000000000000e-05,7.49000000000000e-05,7.99000000000000e-05,8.21000000000000e-05,7.78000000000000e-05,7.33000000000000e-05,6.92000000000000e-05,6.48000000000000e-05,6.03000000000000e-05,5.55000000000000e-05,5.06000000000000e-05,4.55000000000000e-05,4.24000000000000e-05,3.92000000000000e-05,3.52000000000000e-05,3.13000000000000e-05,2.63000000000000e-05,2.05000000000000e-05,1.40000000000000e-05,8.02000000000000e-06,9.15000000000000e-07];

Param.Wing.Root_Chord=6;

Param.Wing.Thickness=x;

Param.Wing.StringerPitch=0.24;

Param.Wing.LE_Sweep=27; 

Param.Wing.TE_Sweep1=0;

Param.Wing.TE_Sweep2=16.5;

Param.Wing.BeamLoc=0.4;

Param.Wing.Kink=0.27;

Param.Wing.AeroPanel_Number=10;

Param.Wing.AeroPanel_AR=2.5;

Param.Wing.Span=35.8;

Param.Wing.Semi_Span=(Param.Wing.Span-4)/2;

Param.Wing.Dihedral=5;

Param.Wing.TotalArea=126;

Param.Wing.HalfArea=(Param.Wing.TotalArea-Param.Wing.Root_Chord*4)/2;




Param.FWT.Fold_angle=10;

Param.FWT.Flare_angle=25;

Param.FWT.Fold_eta=0.75;

Param.FWT.Root_Chord=0.5;

Param.FWT.Root_Height=0.5;

Param.FWT.Tip_Chord=0.3;

Param.FWT.Tip_Height=0.1;

Param.FWT.Thickness=0.002;

Param.FWT.Hinge_Stiffness=1e-4;



Param.Layout.Fuselage_Length=45;

Param.Layout.Wing_Position=20;

Param.Layout.Engine_Position=4.29;

Param.Layout.Horizontal_Tail_Position=42;

Param.Layout.Vertical_Tail_Position=41;


Param.Masses.Secondary_Mass=1000;

Param.Masses.Fuselage_Structure=25000;

Param.Masses.Payload=25000;

Param.Masses.Engine=7362/2;

Param.Masses.Pylon=1239/2;

Param.Masses.MTOW=93500;

Param.Masses.OWE=48500; 

Param.Masses.Fuselage_shell_mass = 2*pi*2*0.004*2800*44.5;

Param.Masses.Horizontal_tail=682; 

Param.Masses.Vertical_tail=522;

Param.Masses.Fuel_Capacity=32940; % L

Param.Masses.Fuel_Fraction=0.2;

Param.Masses.Fuel_Density=840; %g/L

Param.Masses.Fuel_Mass=Param.Masses.Fuel_Capacity * Param.Masses.Fuel_Fraction * Param.Masses.Fuel_Density/1000; 

%% check model 

%  [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models(Param)
 
 

%% run_folder

    run_folder = [
        'C:\Git\A321_sizing\hg_codes\results\test_temp1']; %[-], folder for exporting the NASTRAN model
    
    %% Material properties : Al7075
    
    Yield_strength = 5.2e8;
    
    %% Saftry factor
    
    SF=1;
    
    %% Initiate the loop
 
    Wing_total_mass=1;
    Guess_wing_mass=2;
    
        
    %% start the loop
    
    
    counter=1;
    
    % mass loop
    while abs(Wing_total_mass/Guess_wing_mass - 1)>=0.05
        
        
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models(Param);
        

        [Wingbox_mass, Secondary_mass, Wing_total_mass]=Wing_masses_v1(Param.Wing.Thickness,Box_dimensions.Inboard.Height,Box_dimensions.Inboard.Width, Y_eta, Param.Masses.MTOW, Param.Wing.HalfArea, Param.Wing.LE_Sweep, Param.Wing.Semi_Span);
        
        % update secondary mass
        
        Param.Masses.Secondary_Mass=Secondary_mass;
        
        Guess_wing_mass=Wing_total_mass;
        
        % update fuselage structual mass
        
        Param.Masses.Fuselage_Structure = Param.Masses.OWE - 2*Guess_wing_mass - Param.Masses.Horizontal_tail - Param.Masses.Vertical_tail - Param.Masses.Engine*2 - Param.Masses.Pylon*2 - Param.Masses.Fuselage_shell_mass;
        
        % size loop
        %initial condition
        cond_set1=[1.1,1.1,1.1,1.1];
        cond_set2=[1.1,1.1,1.1,1.1];
        
        Numsec=25;
        
        record=zeros(100,75);
        
        
        while max(cond_set1)>1.02 || min(cond_set2)<0.98
            
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
        
        
        [Wingbox_mass, Secondary_mass, Wing_total_mass]=Wing_masses_v1(Param.Wing.Thickness,Box_dimensions.Inboard.Height,Box_dimensions.Inboard.Width, Y_eta, Param.Masses.MTOW, Param.Wing.HalfArea, Param.Wing.LE_Sweep, Param.Wing.Semi_Span);
         
        Aircraft_total_mass=Param.Masses.Fuselage_Structure + Param.Masses.Fuselage_shell_mass+Wing_total_mass*2 + Param.Masses.Horizontal_tail + Param.Masses.Vertical_tail + Param.Masses.Engine*2 + Param.Masses.Pylon*2;
        
        progress_indicator=abs(Wing_total_mass/Guess_wing_mass - 1);
        
        disp(progress_indicator);
        
        
        
    end
    
    
    %% test
    Param=eval('A321');
    
    run_folder = ['C:\Git\A321_sizing\hg_codes\results\test_temp1']; %[-], folder for exporting the NASTRAN model
        
    [CDi,CD0,CL,k,Distribution]=Drag_Calculation(Param, run_folder);

   [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models(Param);
   
   
   figure 
   plot(Distribution.Y,Distribution.WJ,'b.')
   
      figure 
   plot(Distribution.Y,Distribution.Alphai,'b.')
   
   figure 
   plot(Distribution.Y,Distribution.Chord_Length,'b.')
   
   figure 
   plot(Distribution.Y,Distribution.Lift_var,'b.')
   
   figure
   plot(Distribution.Y,Distribution.CDi_var,'b.')
   
      figure
   plot(Distribution.Y,Distribution.CL_var,'b.')
   
   
    
    figure
    plot(Y_all,Load_distribution.Moment_P2.LC1,'s-')
    hold on
    plot(Y_all,Load_distribution.Moment_P2.LC4,'s-')
    
    
    figure 
    gust_length=linspace(18,214,7);
    plot(gust_length,Root_Delta.Moment_Max.LC1,'s-')
    hold on 
    plot(gust_length,Root_Delta.Moment_Max.LC2,'v-')
    hold on 
    plot(gust_length,Root_Delta.Moment_Max.LC3,'o-')
    
    hold on
    
    plot(gust_length,Root_Delta.Moment_Min.LC1,'s-')
    hold on
    plot(gust_length,Root_Delta.Moment_Min.LC2,'v-')
    hold on
    plot(gust_length,Root_Delta.Moment_Min.LC3,'o-')
    
    
    
    figure
    plot(Y_all(1:end-1),Wing_Delta.Moment_Max.LC1,'s-')
    hold on
    plot(Y_all(1:end-1),Wing_Delta.Moment_Max.LC2,'v-')
    hold on
    plot(Y_all(1:end-1),Wing_Delta.Moment_Max.LC3,'o-')
    
    
    figure
    plot(Y_all,Param.Wing.Thickness(1:25),'s-')

    figure
    plot(Y_all,Param.Wing.Thickness(26:50),'s-')

    figure
    plot(Y_all,Param.Wing.Thickness(51:75),'s-')
    
    
    
