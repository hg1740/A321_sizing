% %% A321 model stretched 

% % Top level control parameters --------------------------

Param=eval('A321');

% Select fold length 
Param.FWT.Fold_eta=0.7;

% Update Wing properties
Param.Wing.AR=16; % 10.172 for A321

Param.Wing.TotalArea=126; % 126m^2 for A321

Res_Name=strcat('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_Damping\Res_AR16_Eta_70_Damp012');

% Update Initial Guess
prevous_result=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR16\Res_AR16_Eta_70_Model');

Param.Wing.Thickness=prevous_result.Param.Wing.Thickness;

if isfield(prevous_result.Param,'FWT')
    
    Param.FWT.Thickness=prevous_result.Param.FWT.Thickness;
    
end

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
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v1(Param);
    
    
else
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
    
end

draw(FEM_full);




%% Run analysis

% Jig twist optimisation 

run_folder = ['C:\Git\A321_sizing\hg_codes\results\AR_Eta_Study']; %[-], folder for exporting the NASTRAN model


 %% Material properties : Al7075
    
    Yield_strength = 5.2e8;
    
    %% Saftry factor
    
    SF=1.5;
        
   %% Sizing Analysis 

    counter=1;

    % size loop
    %initial condition
    cond_set1=[1.1,1.1,1.1,1.1];
    cond_set2=[1.1,1.1,1.1,1.1];
    
    if isfield(Param,'FWT')
        
        Numsec=35;
        record=zeros(100,105);
        
    else
        
        Numsec=25;
        record=zeros(100,75);
        
    end

    % Initialising Jig twist array
    
    [~,CDi,CD0,CL,k,Distribution,Load_distribution, Displacements_Res]=Static_Trim_v1(Param, run_folder, 'Load_Factor',1,'File_Name','Jig_twist_1g','Hinge_Lock','on');

    
    Jig_Twist=zeros(1,length(Distribution.Y));
      
    % Cruising air condition
    
    FlightPoint=awi.model.FlightPoint;
    FlightPoint.Mach=0.78;
    FlightPoint.AcVelocity=FlightPoint.Mach*340;
    FlightPoint.Altitude = 36000;
    Air_condition=getFlightPointData(FlightPoint,'ISA');
        
        
        
    while max(cond_set1)>1.05 || min(cond_set2)<0.95
        
        % Trigger Jig twist 
        
         Discrepancy_L=1;
        
        % Jig twist optimisation

        while max(Discrepancy_L)>0.06
            
            [~,CDi,CD0,CL,k,Distribution,~, Displacements_Res]=Static_Trim_v1(Param, run_folder, 'Load_Factor',1,'File_Name','Jig_twist_1g','Hinge_Lock','on');
            
            Area = trapz(Distribution.Y,Distribution.Lift_var);
            
            PanelY=Distribution.Y;
            
            Panel_Eta=Distribution.Y/(0.5*Param.Wing.Span);
            
            Target_L0=2*Area/(pi*Param.Wing.Span*0.5);
            
            Target_L=Target_L0*sqrt((1-(Distribution.Y/(Param.Wing.Span*0.5)).^2));
            
            Delta_L=Distribution.Lift_var-Target_L;
            
            Discrepancy_L=abs(Delta_L./Target_L);
            
            % Find Eta values for connector + wing part 1 +  wing part 2 + FWT (if exist)
            
            Y_conn=Distribution.Y_conn-Distribution.Y_conn(1);
            
            Y_conn_eta=Y_conn/Y_conn(end);
            
            
            Y_wing=[Distribution.Y_wing1, Distribution.Y_wing2]-Distribution.Y_wing1(1);
            
            Y_wing_eta=Y_wing/Y_wing(end);
            
            
            if isfield(Param,'FWT')
                
                Y_fwt=Distribution.Y_fwt-Distribution.Y_fwt(1);
                
                Y_fwt_eta=Y_fwt/Y_fwt(end);
                
            end
            
            % Update Jig Twist
            
            Delta_Jig=1.2*Delta_L./(Air_condition.AirDensity .* Air_condition.AcVelocity.^2 .* pi.* Distribution.Chord_Length);
            
            Jig_correction = rad2deg(Delta_Jig);
            
            Jig_Twist=Jig_Twist-Jig_correction';
            
            % assign/allocate new twist values to different part of the wing
            
            Jig_Twist_R=Jig_Twist(end/2+1:end);
            
            if isfield(Param,'FWT')
                
                Num_Y=[numel(Y_conn),numel(Y_wing),numel(Y_fwt)];
                
                Cum_Num=cumsum(Num_Y);
                
                Jig_Twist_Connector=Jig_Twist_R(1:Cum_Num(1));
                
                Jig_Twist_Wing=Jig_Twist_R(1+Cum_Num(1):Cum_Num(2));
                
                Jig_Twist_FWT=Jig_Twist_R(1+Cum_Num(2):Cum_Num(3));
                
                
                
                Param.Connector.Jig_Twist=deg2rad([Jig_Twist_Connector]);
                
                Param.Connector.Jig_Eta=Y_conn_eta;
                
                Param.Wing.Jig_Twist=deg2rad([Jig_Twist_Wing]);
                
                Param.Wing.Jig_Eta=Y_wing_eta;
                
                Param.FWT.Jig_Twist=deg2rad([Jig_Twist_FWT]);
                
                Param.FWT.Jig_Eta=Y_fwt_eta;
                
            else
                
                Num_Y=[numel(Y_conn),numel(Y_wing)];
                
                Cum_Num=cumsum(Num_Y);
                
                Jig_Twist_Connector=Jig_Twist_R(1:Cum_Num(1));
                
                Jig_Twist_Wing=Jig_Twist_R(1+Cum_Num(1):Cum_Num(2));
                
                
                
                Param.Connector.Jig_Twist=deg2rad([Jig_Twist_Connector]);
                
                Param.Connector.Jig_Eta=Y_conn_eta;
                
                Param.Wing.Jig_Twist=deg2rad([Jig_Twist_Wing]);
                
                Param.Wing.Jig_Eta=Y_wing_eta;
                
                
            end
            
            disp(max(Discrepancy_L))
            
            
        end
        
        % Run sizing loop with updated Param after Jig twist loop
        
        if isfield(Param,'FWT')
            
            record(counter,:)=[Param.Wing.Thickness,Param.FWT.Thickness(2:11),Param.FWT.Thickness(13:22),Param.FWT.Thickness(24:33)];
            
            [Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis_Fast_v1(Param, run_folder);
            
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
            
            Param.Wing.Thickness(26:50)=Param.Wing.Thickness(26:50).*skin_coefficient(1:25);
            Param.FWT.Thickness(13:22)=Param.FWT.Thickness(13:22).*skin_coefficient(26:35);
            
            % spar check
            Check_von_spar=SF*RVon_spr/Yield_strength;
            
            Spar_coefficient=step_size(Check_von_spar);
            
            Param.Wing.Thickness(1:25)=Param.Wing.Thickness(1:25).*Spar_coefficient(1:25);
            Param.FWT.Thickness(2:11)=Param.FWT.Thickness(2:11).*Spar_coefficient(26:35);
            
            % stringers check
            Check_strg=SF*Rsigma_strg./Rsigma_col;
            
            Stringer_coefficient=step_size(Check_strg);
            
            Param.Wing.Thickness(51:75)=Param.Wing.Thickness(51:75).*Stringer_coefficient(1:25);
            Param.FWT.Thickness(24:33)=Param.FWT.Thickness(24:33).*Stringer_coefficient(26:35);
            
            % end condition check
            cond1=max(Check_von_skin1(1:end-3));
            cond2=max(Check_von_skin2(1:end-3));
            cond3=max(Check_von_spar(1:end-3));
            cond4=max(Check_strg(1:end-3));
            
            Skin_conds=[Check_von_skin1;Check_von_skin2];
            Skin_ub=max(Skin_conds);
            
            cond_skin=min(Skin_ub(1:end-3));
            
            cond5=min(Check_von_skin1(1:end-3));
            cond6=min(Check_von_skin2(1:end-3));
            cond7=min(Check_von_spar(1:end-3));
            cond8=min(Check_strg(1:end-3));
            
            cond_set1=[cond1,cond2,cond3,cond4];
            cond_set2=[cond_skin,cond7,cond8];
            
            disp(cond_set1);
            disp(cond_set2);
            
            counter=counter+1;
            
        else
            
            record(counter,:)=Param.Wing.Thickness;
            
            
            
            [Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis_Fast_v1(Param, run_folder);
            
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
            
            Skin_conds=[Check_von_skin1;Check_von_skin2];
            Skin_ub=max(Skin_conds);
            cond_skin=min(Skin_ub);
            
            cond5=min(Check_von_skin1);
            cond6=min(Check_von_skin2);

            cond7=min(Check_von_spar);
            cond8=min(Check_strg);
            
            cond_set1=[cond1,cond2,cond3,cond4];
            cond_set2=[cond_skin,cond7,cond8];
            
            disp(cond_set1);
            disp(cond_set2);
            
            counter=counter+1;
            
        end
        
        
    end
    
    
    %% Write Results Output Files 
    
    % Model data 
    
    save(strcat(Res_Name,'_Model.mat'), 'Param')
    
    % Aerodynamics data L/D
    
    Aerodynamics.Y=Distribution.Y;
    
    Aerodynamics.CL=CL;
    
    Aerodynamics.CDi=CDi;
    
    Aerodynamics.CD0=CD0;
    
    Aerodynamics.k=k;
    
    Aerodynamics.WJ=Distribution.WJ;
    
    Aerodynamics.Alphai=Distribution.Alphai;
    
    Aerodynamics.Chord_Length=Distribution.Chord_Length;
    
    Aerodynamics.CL_var=Distribution.CL_var;
    
    Aerodynamics.Cd_var=Distribution.CDi_var;
    
    Aerodynamics.Lift=Distribution.Lift_var;
    
    Aerodynamics.Panel_width=Distribution.Panel_width;
    
    save(strcat(Res_Name,'_Aero.mat'), 'Aerodynamics')
   
   
    % Wingbox structural data
    
    Wingbox_property.NodesY=Y_all;
    
    Wingbox_property.Bwidth=Box_dimensions.Inboard.Width';
    
    Wingbox_property.Bheight=Box_dimensions.Inboard.Height';
    
    Wingbox_property.Izz=Box_CrossSec.Izz';
    
    Wingbox_property.Ixx=Box_CrossSec.Ixx';
    
    Wingbox_property.A=Box_CrossSec.A';
    
    
    if isfield(Param,'FWT')
        
        Wingbox_property.Spar_Thickness=[Param.Wing.Thickness(1:25),Param.FWT.Thickness(2:11)]';
        
        Wingbox_property.Skin_Thickness=[Param.Wing.Thickness(26:50),Param.FWT.Thickness(13:22)]';
        
        Wingbox_property.Strg_Area=[Param.Wing.Thickness(51:75),Param.FWT.Thickness(24:33)]';
        
    else
        
        Wingbox_property.Spar_Thickness=Param.Wing.Thickness(1:25)';
        
        Wingbox_property.Skin_Thickness=Param.Wing.Thickness(26:50)';
        
        Wingbox_property.Strg_Area=Param.Wing.Thickness(51:75)';
        
    end
    

%     a1 = Wingbox_property.A(1:end-1);
%     a2 = Wingbox_property.A(2:end);
%     mean_area=(a1+a2)/2;
%     
%     l2=Wingbox_property.NodesY(2:end);
%     l1=Wingbox_property.NodesY(1:end-1);
%     seg_length=l2-l1;
%     
%     Vol_wing=sum(mean_area.*seg_length);
%     
%     rho_al=2800;
%     
%     Wingbox_mass=Vol_wing*rho_al;

    [Wingbox_mass,Secondary_mass,Total_mass]=Wing_Mass_Calc_v1(Wingbox_property, Param);

    Wingbox_property.Wingbox_Mass=Wingbox_mass;
    Wingbox_property.Secondary_Mass=Secondary_mass;
    Wingbox_property.Total_Mass=Total_mass;
    
    
    save(strcat(Res_Name,'_Structure.mat'),'Wingbox_property');
    
    % Load distribution 
    
    Loads.Y=Y_all;
    
    Loads.Moment.pullup_cruise=Load_distribution.Moment_P2.LC1;
    
    Loads.Moment.pullup_sealevel=Load_distribution.Moment_P2.LC2;
    
    Loads.Moment.dive_sealevel=Load_distribution.Moment_P2.LC3;
    
    Loads.Moment.g_sealevel_gust1=Load_distribution.Moment_P2.LC6;
    
    Loads.Moment.g_sealevel_gust2=Load_distribution.Moment_P2.LC9;
    
    Loads.Moment.hinge_lock=Load_distribution.Moment_P2.HL;
    
    Loads.Moment.zero_fuel=Load_distribution.Moment_P2.ZF;
    
    
    Loads.Shear.pullup_cruise=Load_distribution.Shear_P2.LC1;
    
    Loads.Shear.pullup_sealevel=Load_distribution.Shear_P2.LC2;
    
    Loads.Shear.dive_sealevel=Load_distribution.Shear_P2.LC3;
    
    Loads.Shear.g_sealevel_gust1=Load_distribution.Shear_P2.LC6;
    
    Loads.Shear.g_sealevel_gust2=Load_distribution.Shear_P2.LC9;
    
    Loads.Shear.hinge_lock=Load_distribution.Shear_P2.HL;
    
    Loads.Shear.zero_fuel=Load_distribution.Shear_P2.ZF;
    
    
    Loads.Torque.pullup_cruise=Load_distribution.Torque.LC1;
    
    Loads.Torque.pullup_sealevel=Load_distribution.Torque.LC2;
    
    Loads.Torque.dive_sealevel=Load_distribution.Torque.LC3;
    
    Loads.Torque.g_sealevel_gust1=Load_distribution.Torque.LC6;
    
    Loads.Torque.g_sealevel_gust2=Load_distribution.Torque.LC9;
    
    Loads.Torque.hinge_lock=Load_distribution.Torque.HL;
    
    Loads.Torque.zero_fuel=Load_distribution.Torque.ZF;
    
    if isfield(Param,'FWT')
        Loads.HL=Load_distribution.HL;
    end
    
    Loads.ZF=Load_distribution.ZF;
    
    
    Loads.Wing_Delta=Wing_Delta;
    
    Loads.Root_Delta.Gust_Length=linspace(18,214,7);
    
    Loads.Root_Delta=Root_Delta;
    
    Loads.Displacement=Displacements_Res;
   
    save(strcat(Res_Name,'_Loads.mat'),'Loads');
    
    
%     %% Wing box weight record
%     
%     Wingbox_Weight_Recard(i)=Wingbox_mass;
    
% 
% end
%     
%     toc
%     
    
    
    
%     figure 
%     
%     for i=6:12
%         plot(record(i,1:25))
%         hold on
%     end
    
%     figure 
%     
%     plot(Loads.Y(1:end-1),Loads.Wing_Delta.Moment_Max.LC3,'bs')
%     hold on 
%     plot(Loads.Y(1:end-1),Loads.Wing_Delta.Moment_Min.LC3,'bs')
    
% figure    
% for i=25:60
% plot(Y_all(1:25),record(i,1:25))
% hold on 
% end

% figure
% 
% plot(Param.FWT.Thickness(1:11))
% 
% figure 
% plot(Y_all(1:25),record(50,1:25),'r-')
% hold on 
% plot(Y_all(1:25),record(60,1:25),'b-')
% hold on 
% plot(Y_all(1:25),record(70,1:25),'k-')

%% plot results
% % Run_folder=strcat(run_folder,'\hinge_loacked');
% model = mni.import_matran(fullfile(run_folder,'a321_cruise_2p5g.dat'));
% model.draw
% 
% %extract the data
% f06 = mni.result.f06(fullfile(run_folder,'a321_cruise_2p5g.f06'));
% res_disp =  f06.read_disp;
% res_aeroP = f06.read_aeroP;
% res_aeroF = f06.read_aeroF;
% 
% % apply deformation result
% [~,i] = ismember(model.GRID.GID,res_disp.GP);
% model.GRID.Deformation = [res_disp.dX(:,i);res_disp.dY(:,i);res_disp.dZ(:,i)];
% 
% % apply aero pressure
% model.CAERO1.PanelPressure = res_aeroP.Cp;
% 
% %apply aero forces
% f = [res_aeroF.aeroFx;res_aeroF.aeroFy;res_aeroF.aeroFz;...
%     res_aeroF.aeroMx;res_aeroF.aeroMy;res_aeroF.aeroMz];
% model.CAERO1.PanelForce = f';
% 
% % update the plot to apply deformations and aero pressures + forces
% model.update('Scale',1)
% 
% 
% 
% figure
% 
% plot(Distribution.Y,Distribution.Lift_var,'b.')
% 
% figure
% 
% plot(Y_all,Load_distribution.Moment_P2.LC1,'b.')
% 
% hold on
% plot(Y_all,Load_distribution.Moment_P2.Cruise,'r.')
% 
% %% mass calculation 

% [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v1(Param);
% 
% Thickness_sparall=[Param.Wing.Thickness(1:25),Param.FWT.Thickness(2:11)];
% Thickness_skinall=[Param.Wing.Thickness(26:50),Param.FWT.Thickness(13:22)];
% Area_strgall=[Param.Wing.Thickness(51:75),Param.FWT.Thickness(24:33)];
% 
% 
% FWT_y_eta=Y_eta(end)+linspace(FWT.RData(1),FWT.RData(2),11);
% 
% Eta_all=[Y_eta,FWT_y_eta(2:end)];
% 
% Box_height=Box_dimensions.Inboard.Height;
% 
% Box_width=Box_dimensions.Inboard.Width;
% 
% Num_stringers=Box_width/0.24;
% 
% cs_=(Thickness_sparall.* Box_height + Thickness_skinall.* Box_width)*2 + Area_strgall.*Num_stringers*2;
% a1=cs_(1:end-1);
% a2=cs_(2:end);
% mean_area=(a1+a2)/2;
% 
% l2=Eta_all(2:end);
% l1=Eta_all(1:end-1);
% seg_length=l2-l1;
% 
% Vol_wing=sum(mean_area.*seg_length);
% rho_al=2700;
% 
% Wingbox_mass=Vol_wing*rho_al;
% 
% 
% figure 
% plot(Eta_all,Thickness_sparall,'rs','MarkerFaceColor','r')
% 
% xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
% 
% ylabel('Spar thickness (m)','Interpreter','latex','FontSize',12)
% 
% set(gcf,'Color','w');
% 
% figure 
% plot(Eta_all,Area_strgall,'rs','MarkerFaceColor','r')
% 
% xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
% 
% ylabel('Stringer area (m$^2$)','Interpreter','latex','FontSize',12)
% 
% set(gcf,'Color','w');
% 
% figure 
% 
% plot(Eta_all,Thickness_skinall,'rs','MarkerFaceColor','r')
% 
% xlabel('Spanwise distance (m)','Interpreter','latex','FontSize',12)
% 
% ylabel('Skin thickness (m)','Interpreter','latex','FontSize',12)
% 
% set(gcf,'Color','w');
% 
% %% Drag compare
% 
% polar_Cl=0:0.05:0.9;
% 
% Reference_Cd=0.017+0.037*polar_Cl.^2;
% 
% New_Cd=0.017+0.0154*polar_Cl.^2;
% 
% figure 
% plot(Reference_Cd,polar_Cl,'bs','MarkerFaceColor','b')
% hold on 
% plot(New_Cd,polar_Cl,'rs','MarkerFaceColor','r')
% 
% xlabel('Drag coefficient $C_D$','Interpreter','latex','FontSize',12)
% 
% ylabel('Lift coefficient $C_L$','Interpreter','latex','FontSize',12)
% 
% set(gcf,'Color','w');
% legend('Reference Aircraft','AR 19 Model','Interpreter','latex','FontSize',12)
% 
% 
% figure 
% plot(Y_all,Load_distribution.Shear_P2.LC1,'bs')
% hold on 
% plot(Y_all,Load_distribution.Shear_P2.LC2,'bs')
% hold on 
% plot(Y_all,Load_distribution.Shear_P2.LC3,'rs')
% hold on 
% plot(Y_all,Load_distribution.Shear_P2.LC6,'bs')
% % hold on 
% % plot(Y_all,Load_distribution.Shear_P2.LC9,'bs')
% hold on 
% plot(Y_all,Load_distribution.Shear_P2.Cruise,'bs')



