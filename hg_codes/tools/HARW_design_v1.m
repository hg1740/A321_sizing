%% A321 model stretched 


% Top level control parameters --------------------------

Param=eval('A321');

% % Select fold length 
Param.FWT.Fold_eta=0.8;

% Update Wing properties
Param.Wing.AR=10.172;

Param.Wing.TotalArea=126;

%---------------------------------------------------------



if Param.FWT.Fold_eta==1
    
    Param=rmfield(Param,'FWT')
    
end


% Generate wing planform properties 

if isfield(Param,'FWT') 
    [Geo_Wing, Geo_FWT]= Wing_Gen_V1(Param)
else
    [Geo_Wing]= Wing_Gen_V1(Param)
end

% Update Wing planform 

Param.Wing.Root_Chord=Geo_Wing.Root_Chord;

Param.Wing.TE_Sweep1=Geo_Wing.TE_Sweep1;

Param.Wing.TE_Sweep2=Geo_Wing.TE_Sweep2;

Param.Wing.Span=Geo_Wing.Span;

Param.Wing.Semi_Span=(Param.Wing.Span-4)/2;

Param.Wing.HalfArea=(Param.Wing.TotalArea-Param.Wing.Root_Chord*4)/2;

% Initial Jig Twist 

Param.Connector.Jig_Twist=deg2rad([0,0]);

Param.Connector.Jig_Eta=[0,1];

Param.Wing.Jig_Twist=deg2rad([0,0]);

Param.Wing.Jig_Eta=[0,1];


if isfield(Param,'FWT') 
    
    % FWT update geometry
    
    Param.FWT.Root_Chord=Geo_FWT.Root_Chord;
    
    Param.FWT.Root_Height=Geo_FWT.Root_Height;
    
    Param.FWT.Tip_Chord=Geo_FWT.Tip_Chord;
    
    Param.FWT.Tip_Height=Geo_FWT.Tip_Height;

    % update initial jig twist
    
    Param.FWT.Jig_Twist=deg2rad([0,0]);
    
    Param.FWT.Jig_Eta=[0,1];
    
    % Build Aircraft FEM 
    
%     [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v1(Param);
%     
% else
%     
%     [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
    
end




%% Run analysis

% Jig twist optimisation 


run_folder = ['C:\Git\A321_sizing\hg_codes\results\test_temp']; %[-], folder for exporting the NASTRAN model

[Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis_Fast(Param, run_folder)



[CDi,CD0,CL,k,Distribution,Load_distribution]=Drag_Calculation(Param, run_folder);

% Initialising Jig twist array

Jig_Twist=zeros(1,length(Distribution.Y));

Discrepancy_L=1;

% Cruising air condition 

FlightPoint=awi.model.FlightPoint;
FlightPoint.Mach=0.78;
FlightPoint.AcVelocity=FlightPoint.Mach*340;
FlightPoint.Altitude = 36000;
Air_condition=getFlightPointData(FlightPoint,'ISA');


while max(Discrepancy_L)>0.05
    
    [CDi,CD0,CL,k,Distribution]=Drag_Calculation(Param, run_folder);
    
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



 %% Material properties : Al7075
    
    Yield_strength = 5.2e8;
    
    %% Saftry factor
    
    SF=1.5;
        
   %% start the loop

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
    
    [CDi,CD0,CL,k,Distribution,Load_distribution]=Drag_Calculation(Param, run_folder);
    
    
    Jig_Twist=zeros(1,length(Distribution.Y));

    
    % Cruising air condition
    
    FlightPoint=awi.model.FlightPoint;
    FlightPoint.Mach=0.78;
    FlightPoint.AcVelocity=FlightPoint.Mach*340;
    FlightPoint.Altitude = 36000;
    Air_condition=getFlightPointData(FlightPoint,'ISA');
        
        
        
    while max(cond_set1)>1.05 || min(cond_set2)<0.95
        
         Discrepancy_L=1;
        
        % Jig twist optimisation

        while max(Discrepancy_L)>0.05
            
            [CDi,CD0,CL,k,Distribution]=Drag_Calculation(Param, run_folder);
            
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
        
        
        
        record(counter,:)=Param.Wing.Thickness;
        
        [Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis_Fast(Param, run_folder);
        
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
        
        













% figure
% 
% plot(Y_all-2, Load_distribution.Moment_P2.LC1,'b-s','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Moment_P2.LC2,'r-o','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Moment_P2.LC3,'k-v','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Moment_P2.LC6,'b-v','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Moment_P2.LC9,'r-o','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Moment_P2.Cruise,'k-','LineWidth',1)
% 
% 
% legend('2.5g 36000 ft','2.5g 3000 ft','-g 3000 ft','g +1MC 3000 ft','g - 1MC 3000 ft','Cruise','Interpreter','latex','FontSize',12)
% 
% set(gcf,'Color','w')
% xlabel('Span distance (m)','Interpreter','latex','FontSize',12)
% ylabel('Bending Moment (Nm)','Interpreter','latex','FontSize',12)
% 
% 
% figure
% 
% plot(Y_all-2, Load_distribution.Shear_P2.LC1,'b-s','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Shear_P2.LC2,'r-o','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Shear_P2.LC3,'k-v','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Shear_P2.LC6,'b-v','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Shear_P2.LC9,'r-o','LineWidth',0.8)
% hold on
% plot(Y_all-2, Load_distribution.Shear_P2.Cruise,'k-','LineWidth',1)
% 
% 
% legend('2.5g 36000 ft','2.5g 3000 ft','-g 3000 ft','g +1MC 3000 ft','g - 1MC 3000 ft','Cruise','Interpreter','latex','FontSize',12)
% 
% set(gcf,'Color','w')
% xlabel('Span distance (m)','Interpreter','latex','FontSize',12)
% ylabel('Shear Force (N)','Interpreter','latex','FontSize',12)








% 
% figure 
% plot(Load_distribution.Y,Load_distribution.Moment_P2,'b.')
% 
% 
% figure 
% plot(Load_distribution.Y,Load_distribution.Shear_P2,'b.')








% figure 
% 
% % plot(Distribution.Y,Distribution.Lift_var,'r.','MarkerSize',8)
% % 
% % hold on 
% 
% plot(Distribution.Y,Target_L,'k--','LineWidth',1)
% 
% legend('Lift distribution on the wing','Targeted lift distribution','Interpreter','latex','FontSize',12)
% 
% set(gcf,'Color','w')
% xlabel('Span distance (m)','Interpreter','latex','FontSize',12)
% ylabel('Lift per unit span (N/m)','Interpreter','latex','FontSize',12)
% axis([-20 20 0 4e4])
% 
% 
% figure 
% plot(Distribution.Y,-Distribution.Alphai,'b.')
% hold on 
% plot([-20,20],[0.023,0.023],'k--','LineWidth',1)
% set(gcf,'Color','w')
% xlabel('Span distance (m)','Interpreter','latex','FontSize',12)
% ylabel('Induced AoA (Radius)','Interpreter','latex','FontSize',12)
% axis([-20 20 0 0.04])
% legend('Calculated result','Analytical prediction','Interpreter','latex','FontSize',12)
% 
% 
% 
% figure 
% plot(Distribution.Y,-Distribution.CDi_var,'b.')
% 
% set(gcf,'Color','w')
% xlabel('Span distance (m)','Interpreter','latex','FontSize',12)
% ylabel('Induced drag coefficient','Interpreter','latex','FontSize',12)
% % axis([-20 20 0 0.04])
% legend('Calculated result','Analytical prediction','Interpreter','latex','FontSize',12)

% %% plot results
% 
% model = mni.import_matran(fullfile(run_folder,'a321_36000ft_1g_drag.dat'));
% model.draw
% 
% %extract the data
% f06 = mni.result.f06(fullfile(run_folder,'a321_36000ft_1g_drag.f06'));
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



% %% test drag calculation - elliptical
% 
% % calculate vorticity
%     % Air conditions
%     FlightPoint=awi.model.FlightPoint;
%     FlightPoint.Mach=0.78;
%     FlightPoint.AcVelocity=0.78*340;
%     FlightPoint.Altitude = 36000;
%     getFlightPointData(FlightPoint,'ISA');
%     DynPressure = FlightPoint.DynPressure;
%     
%     
%     
%     Lift_all=Target_L;
%     
%     Y_all=Distribution.Y;
%     
%     panel_width_all=Distribution.Panel_width;
%     
%     Lift_all_abs=Lift_all.*panel_width_all';
%     
%     Gamma_all=Lift_all/(FlightPoint.AirDensity*FlightPoint.AcVelocity);
%     
%     % calculate the derivative
%     dGdy= gradient(Gamma_all(:)) ./ gradient(Y_all(:));
%     
%     % downwash
%     wj=zeros(1,length(Y_all));
%     alphai=zeros(1,length(Y_all));
%     
%     for i=1:length(Y_all)
%         
%         w=-(1/(4*pi))* dGdy.* panel_width_all'./(Y_all(i)-Y_all);
%         
%         w=w( ~any( isnan( w ) | isinf( w ), 2 ),: );
%         
%         wj(i)=sum(w);
%         
%         alphai(i)=wj(i)/FlightPoint.AcVelocity;
%         
%     end
%     
%     Dragi_var=Lift_all_abs'.*sin(alphai);
%     Drag_force=sum([Dragi_var]);
%     
%     %     if ~isempty(FWT_)
%     %         Total_area=(Wingbox_right.SurfaceArea+FWT_R.SurfaceArea)*2 + Param.Wing.Root_Chord*Param.Layout.Fuselage_Width;
%     %     else
%     %         Total_area=(Wingbox_right.SurfaceArea)*2 + Param.Wing.Root_Chord*Param.Layout.Fuselage_Width;
%     %     end
%     
%     Total_area=126;
%     
%     Cdi=Drag_force/(DynPressure*Total_area);
%     
%     figure
%     plot(Y_all,wj,'b.')
%     
%     figure
%     plot(Y_all,dGdy,'b.')
%     
%         figure
%     plot(Y_all,Gamma_all,'b.')
    

% Param.Connector.Jig_Twist=deg2rad([x(1),x(2)]);
% 
% Param.Connector.Jig_Eta=[0,1];
% 
% Param.Wing.Jig_Twist=deg2rad([x(3),x(4),x(5)]);
% 
% Param.Wing.Jig_Eta=[0,0.27,1];

% old code------------------------------------------------------------------
% AllFEM = flatlist(FEM_full);
% Wing_part = AllFEM(ismember({AllFEM.Name}, {'Connector_Right','A320Wing_right','FWT_R'}));
% 
% 
% Num_Panel=zeros(numel(Wing_part),2);
% 
% for i=1:numel(Wing_part)
%     
%     for j=1:numel(Wing_part(i).AeroPanels)
%         
%         Num_Panel(i,j)=Wing_part(i).AeroPanels(j).NumPanels;
%         
%     end
%     
% end
% 
% Total_PanelNumber=sum(Num_Panel,'all');
% 
% Num_Y=Total_PanelNumber/Param.Wing.AeroPanel_Number;


%% fmincon example

% % Initial guess of the twist
% x0=[0,0,0,0,0];
% 
% % Define variable bound
% lb=[-10, -10, -10, -10, -10];
% ub=[10,15,15,15, 15];
% 
% % % % Settings
% % options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','CheckGradients',true);
% 
% options = optimset('FinDiffRelStep',1e-3);
% % 
% % 
% % % Run optimisation 
% [xopt, fva,ef,output,lamda] = fmincon(@objective,x0,[],[],[],[],lb,ub,[],options);
% 
% % [x,fval] = fminunc(@objective,x0);
% 
% 
% % handle=@objective;
% % 
% % obj=handle(x0);
% 
% 
% % Define objective function
% 
% 
% function delta =Find_Diff(x)
% 
% Param=eval('A321');
% 
% run_folder = ['C:\Git\A321_sizing\hg_codes\results\test_temp']; %[-], folder for exporting the NASTRAN model
% 
% Param.Connector.Jig_Twist=deg2rad([x(1),x(2)]);
% 
% Param.Connector.Jig_Eta=[0,1];
% 
% Param.Wing.Jig_Twist=deg2rad([x(3),x(4),x(5)]);
% 
% Param.Wing.Jig_Eta=[0,0.27,1];
% 
% [CDi,CD0,CL,k,Distribution]=Drag_Calculation(Param, run_folder);
% 
% Area = trapz(Distribution.Y,Distribution.Lift_var);
% 
% Target_L0=2*Area/(pi*Param.Wing.Span*0.5);
% 
% Target_L=Target_L0*sqrt(1-(Distribution.Y/(Param.Wing.Span*0.5)).^2);
% 
% delta_L=abs(Distribution.Lift_var-Target_L);
% 
% delta=sum(delta_L);
% 
% end 
% 
% 
% function obj = objective(x)
% 
% obj=1e5*Find_Diff(x);
% 
% disp(obj)
% 
% disp(x)
% 
% end
% 
% 
% % Constraint
% function [c, ceq]=constraint(x)
% 
% c = x(1)+x(2)+x(3)+x(4)-100;
% 
% ceq=[];
% 
% end
% 
% 
% 
% x=[-0.761,0.2518,-1.54,2.05,2.176];
% 
% Param=eval('A321');
% 
% run_folder = ['C:\Git\A321_sizing\hg_codes\results\test_temp']; %[-], folder for exporting the NASTRAN model
% 
% Param.Connector.Jig_Twist=deg2rad([x(1),x(2)]);
% 
% Param.Connector.Jig_Eta=[0,1];
% 
% Param.Wing.Jig_Twist=deg2rad([x(3),x(4),x(5)]);
% 
% Param.Wing.Jig_Eta=[0,0.27,1];
% 
% [CDi,CD0,CL,k,Distribution]=Drag_Calculation(Param, run_folder);
% 
% Area = trapz(Distribution.Y,Distribution.Lift_var);
% 
% Target_L0=2*Area/(pi*Param.Wing.Span*0.5);
% 
% Target_L=Target_L0*sqrt(1-(Distribution.Y/(Param.Wing.Span*0.5)).^2);
% 
% figure 
% 
% plot(Distribution.Y,Distribution.Lift_var,'b.','MarkerSize',8)
% 
% hold on 
% 
% plot(Distribution.Y,Target_L,'rs')
% 
% set(gcf,'Color','w')
% xlabel('Span distance (m)','Interpreter','latex','FontSize',12)
% ylabel('Lift per unit span (N/m)','Interpreter','latex','FontSize',12)





    
% %    a=area(Distribution.Y,Distribution.Lift_var)
%    
%    A = trapz(Distribution.Y,Distribution.Lift_var);
%    
%    Target_L0=2*A/(pi*Param.Wing.Span*0.5);
%    
% 
%    
%    Target_L=Target_L0*sqrt(1-(Distribution.Y/(Param.Wing.Span*0.5)).^2);
%    
%    delta_L=abs(Distribution.Lift_var-Target_L);
%    
%    figure 
%    plot(Distribution.Y,Distribution.Lift_var,'b.')
%    
%    hold on 
%    
%    plot(Distribution.Y,Target_L,'rs')
%    
%    figure 
%    plot(Distribution.Y,delta_L,'rs')
   

   
   

% %% plot results
% 
% model = mni.import_matran(fullfile(run_folder,'a321_36000ft_1g_drag.dat'));
% model.draw
% 
% %extract the data
% f06 = mni.result.f06(fullfile(run_folder,'a321_36000ft_1g_drag.f06'));
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


