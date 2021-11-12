function Param_opt = Jig_Twist(Param,run_folder) 



    % Initialising Jig twist array
    
    [~,~,~,~,~,Distribution,~, ~]=Static_Trim_v1(Param, run_folder, 'Load_Factor',1,'File_Name','Jig_twist_1g','Hinge_Lock','on');

    
    Jig_Twist=zeros(1,length(Distribution.Y));
      
    % Cruising air condition
    
    FlightPoint=awi.model.FlightPoint;
    FlightPoint.Mach=0.78;
    FlightPoint.AcVelocity=FlightPoint.Mach*340;
    FlightPoint.Altitude = 36000;
    Air_condition=getFlightPointData(FlightPoint,'ISA');
    
    

% Trigger Jig twist 
        
         Discrepancy_L=1;
        
        % Jig twist optimisation

        while max(Discrepancy_L)>0.1
            
            [~,~,~,~,~,Distribution,~,Displacements_Res,~, ~]=Static_Trim_v1(Param, run_folder, 'Load_Factor',1,'File_Name','Jig_twist_1g','Hinge_Lock','on');
            
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
        
        % optimised param 
        
        Param_opt=Param;




end