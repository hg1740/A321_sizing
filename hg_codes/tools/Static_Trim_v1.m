function [FEM_full,CDi,CD0,CL,k,Aerodynamic_distribution,Load_distribution,Displacements_Res,Box_dimensions, Box_CrossSec]=Static_Trim_v1(Param, run_folder, varargin)


    p=inputParser;

    addParameter(p,'Load_Factor',1,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))
    addParameter(p,'File_Name','file_name',@(x)validateattributes(x,{'char'},{'nonempty'}))
    addParameter(p,'Hinge_Lock','on',@(x)any(validatestring(x,{'on','off'})))
    addParameter(p,'Altitude',36000,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))
    addParameter(p,'Mach_Num',0.78,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))

    parse(p,varargin{:})

    % update Panel aspect ratio 
    
    Param.Wing.AeroPanel_AR=1;
    
    if ~isempty(p.Results.Hinge_Lock) && isfield(Param,'FWT')
        
        if contains(p.Results.Hinge_Lock,'on')
            
            % $$$---  Hinge lock   ---$$$
            
            Param.FWT.Hinge_Stiffness=1e14;
            
            % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            % create a sub-directory for hinge-locked cases
            
            if ~exist(strcat(run_folder,'\hinge_locked'),'dir')
                
                mkdir(run_folder,'hinge_locked')
                
                run_folder=strcat(run_folder,'\hinge_locked');
                
            else
                run_folder=strcat(run_folder,'\hinge_locked');
            end
            
            
        end
        
    end
    
    % create Aircraft model
    
    if isfield(Param,'FWT')
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT_R]=Aircraft_Models_v3(Param);
        
    else
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v3(Param);
    end
    
    
    export(FEM_full, run_folder);
       
     %% steady state for cruising, 36000 ft, g
    TrimLoadcase = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = p.Results.Altitude;
    mach_number       = p.Results.Mach_Num;
    aircraft_velocity = mach_number*340;
    flap_angle=0;
    
    TrimLoadcase.Name = p.Results.File_Name;
    
    TrimLoadcase.Altitude   = altitude;
    TrimLoadcase.Mach       = mach_number;
    TrimLoadcase.AcVelocity = aircraft_velocity;
    TrimLoadcase.AcMass = acMass;
    
    TrimLoadcase.PitchAngle=0;
    TrimLoadcase.RollAngle =0;
    TrimLoadcase.ID = 1030;
    TrimLoadcase.LoadFactor = p.Results.Load_Factor;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase.CsDeflection=flap_angle*pi/180;
    
    
    %% NASTRAN method - RUN SOL 144
    
    NastranMethods = awi.methods.Nastran;
    NastranMethods.AnalysisModel = FEM_full;
    MassCases=awi.model.MassCases.empty;
    
    % AELINK Cards
    % get control surfaces
    controlsurfs = awi.methods.Nastran.getControlSurfs(FEM_full);
    elev_surfs = controlsurfs(cellfun(@(x)contains(x,'elev'),controlsurfs));
    if length(elev_surfs) == 2
        label_d = elev_surfs{1};
        aelink = mni.printing.cards.AELINK(elev_surfs{1},{{elev_surfs{2},1}});
    else
        error('Expecting only two elevator control surfaces')
    end
    
    %Write trim file
    trimFile = NastranMethods.writeTrimFile(Aircraft, TrimLoadcase,...
        MassCases,run_folder,'DatFilename',p.Results.File_Name,...
        'ExtraCards',aelink);
    
    % clear directory
    
    delete(strcat(run_folder, '\*.xdb'));
    delete(strcat(run_folder, '\*.h5'));
    delete(strcat(run_folder, '\*.log'));
    delete(strcat(run_folder, '\*.f06'));
    delete(strcat(run_folder, '\*.f04'));
    delete(strcat(run_folder, '\*.op4'));
    
    
    NastranMethods.runNastran(trimFile);
    
    %% Load 
    
    Trim1_res=NastranMethods.extractNastranResults(strcat(run_folder,'\',p.Results.File_Name,'.f06'),'ReadF06',true,'ReadHDF5',false);
    
    WingNodes=NastranMethods.WingNode;
    FWTNodes=NastranMethods.FWTNode;
    
    if ~isempty(FWTNodes)
        
        All_nodes=[WingNodes(1:end-1),FWTNodes(1:end)];
        
    else
        
        All_nodes=[WingNodes(1:end)];
        
    end
    
    X_all=[All_nodes.X];
    Y_all=X_all(2,:)';
    
    index1=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[All_nodes(1:end-1).GID]);
    index2=ismember([Trim1_res.f06data.Bendingmoment.LGRID],[All_nodes(end).GID]);
    
    Trim1.M_P1=[Trim1_res.f06data.Bendingmoment.UMPLN1(index1),Trim1_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Trim1.M_P2=[Trim1_res.f06data.Bendingmoment.UMPLN2(index1),Trim1_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Trim1.T=[Trim1_res.f06data.Bendingmoment.UTORQUE1(index1),Trim1_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
    
    Trim1.S_P1=[Trim1_res.f06data.Bendingmoment.USPLN1(index1),Trim1_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Trim1.S_P2=[Trim1_res.f06data.Bendingmoment.USPLN2(index1),Trim1_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    

    Moment_P1_Loadcase1=Trim1.M_P1; % internal load M1 for 2.5 g pull up
    Moment_P2_Loadcase1=Trim1.M_P2; % internal load M2 for 2.5 g pull up
    Torque_Loadcase1=Trim1.T;       % internal load T for 2.5 g pull up
    Shear_P1_Loadcase1=Trim1.S_P1;
    Shear_P2_Loadcase1=Trim1.S_P2;
    
    % correction stresses at hinge
    Moment_P2_Loadcase1(end-10)=(Moment_P2_Loadcase1(end-9)+Moment_P2_Loadcase1(end-11))/2;
    Shear_P2_Loadcase1(end-10)=(Shear_P2_Loadcase1(end-9)+Shear_P2_Loadcase1(end-11))/2;
    Torque_Loadcase1(end-10)=(Torque_Loadcase1(end-9)+Torque_Loadcase1(end-11))/2;
    
    % Result output
    Load_distribution.Y=Y_all;
    Load_distribution.Moment_P2=Moment_P2_Loadcase1';
    Load_distribution.Shear_P2=Shear_P2_Loadcase1';
    Load_distribution.Torque=Torque_Loadcase1';
    
    
    
    %% Displacements
    
    
    structural_nodes = h5read(strcat(run_folder,'\',p.Results.File_Name,'.h5'),'/NASTRAN/INPUT/NODE/GRID');
    Displacements=h5read(strcat(run_folder,'\',p.Results.File_Name,'.h5'),'/NASTRAN/RESULT/NODAL//DISPLACEMENT');
    
    indexi=ismember([Displacements.ID],[WingNodes.GID]);
    
    Disp_Z_Wing=Displacements.Z(indexi);
    
    Disp_Z_all=Disp_Z_Wing;
    
    
    if ~isempty(FWTNodes)
        
        indexo=ismember([Displacements.ID],[FWTNodes.GID]);
        
        Disp_Z_FWT=Displacements.Z(indexo);
        
        Y_FWT=structural_nodes.X(2,indexo);
        
%         Y_FWT=X_FWT.X(2,:)';
        
        Disp_Z_all=[Disp_Z_Wing(1:end);Disp_Z_FWT(2:end)];
        
        Coast_Angle=rad2deg(atan((Disp_Z_FWT(end)-Disp_Z_FWT(1))/(Y_FWT(end)-Y_FWT(1))));
        
    end
    
    % Result out put
    Displacements_Res.Y_Data=Y_all;
    
    
    Displacements_Res.Z_all=Disp_Z_all;
    Displacements_Res.Z_wing=Disp_Z_Wing;
    
    
    if ~isempty(FWTNodes)
        
        Displacements_Res.FWT_Y_Data=Y_FWT;
        Displacements_Res.Coast_Angle=Coast_Angle;
        Displacements_Res.Z_FWT=Disp_Z_FWT;
    end
    

    %% Drag 
    
    % create an object array 
    All_FEM=flatlist(FEM_full);
    
    Wing_conn=All_FEM(ismember({All_FEM.Name}, 'Connector_Right'));
    Wing_inboard = All_FEM(ismember({All_FEM.Name}, 'A320Wing_right'));
    FWT_ = All_FEM(ismember({All_FEM.Name}, 'A320Wing_right_FWT'));
    
    Wing_conn_ID=Wing_conn.AeroSets.E;
    Wing_inboard_ID=Wing_inboard.AeroSets(1).E;
    Wing_outboard_ID=Wing_inboard.AeroSets(2).E;
    
    % aeroforces at each segment
    res_aeroF = mni.result.f06(strcat(run_folder,'/',p.Results.File_Name,'.f06')).read_aeroF;
    
    Index_conn=ismember([res_aeroF.PanelID(1,:)],Wing_conn_ID);
    Index_wing_inboard=ismember([res_aeroF.PanelID(1,:)],Wing_inboard_ID);
    Index_wing_outboard=ismember([res_aeroF.PanelID(1,:)],Wing_outboard_ID);
    
    
    Lift_conn=res_aeroF.aeroFz(Index_conn);
    Lift_wing_inboard=res_aeroF.aeroFz(Index_wing_inboard);
    Lift_wing_outboard=res_aeroF.aeroFz(Index_wing_outboard);
    
    % calculate panel width for each segment: conn + wing_bef_kink +
    % wing_aft_kink
    num_panel_conn=numel(Wing_conn_ID)/Param.Wing.AeroPanel_Number;
    num_panel_inboard=numel(Wing_inboard_ID)/Param.Wing.AeroPanel_Number;
    num_panel_outboard=numel(Wing_outboard_ID)/Param.Wing.AeroPanel_Number;
     
    kink_pos=Wingbox_right.YData(2);
    
    width_conn=0.5 * Param.Layout.Fuselage_Width / num_panel_conn;
    width_wing1=kink_pos/num_panel_inboard;
    width_wing2=(Wingbox_right.YData(3)-Wingbox_right.YData(2))/num_panel_outboard;
    
    
    if ~isempty(FWT_)
        FWT_ID=FWT_.AeroSets.E;
        Index_FWT=ismember([res_aeroF.PanelID(1,:)],FWT_ID);
        Lift_FWT=res_aeroF.aeroFz(Index_FWT);
        num_panel_FWT=numel(FWT_ID)/Param.Wing.AeroPanel_Number;
        width_fwt=FWT_R.Span/num_panel_FWT;
        
        panel_width =[width_conn*ones(1,num_panel_conn),width_wing1*ones(1,num_panel_inboard),width_wing2*ones(1,num_panel_outboard),width_fwt*ones(1,num_panel_FWT)]; 
        
    else
        
        panel_width =[width_conn*ones(1,num_panel_conn),width_wing1*ones(1,num_panel_inboard),width_wing2*ones(1,num_panel_outboard)];

    end
    

    % Y-Data
    % find corresponding Y positions
    Y_conn=0.5*width_conn:width_conn*0.999:2;
    
    Y_wing1=2+0.5*width_wing1:width_wing1*0.999:2+Wingbox_right.YData(2);
    
    Y_wing2=2+Wingbox_right.YData(2)+ 0.5*width_wing2:width_wing2*0.999:2+Wingbox_right.YData(3);
    
    if ~isempty(FWT_)
        
        Y_fwt=2+Wingbox_right.YData(3)+0.5*width_fwt:width_fwt*0.999:2+Wingbox_right.Span+FWT_R.Span;
        
        Y=[Y_conn,Y_wing1,Y_wing2,Y_fwt];
        
    else
        Y=[Y_conn,Y_wing1,Y_wing2];
        
    end
    

    % normalise lift by panel width to obtain lift per unit span
    
    if ~isempty(FWT_)
        lift_wing=[Lift_conn';Lift_wing_inboard';Lift_wing_outboard'; Lift_FWT'];
        
        lift_wing_matrix=reshape(lift_wing, Param.Wing.AeroPanel_Number, numel(lift_wing)/Param.Wing.AeroPanel_Number);
        
        % abs. value of the lift
        lift_wing_abs=sum(lift_wing_matrix);
        
        Cum_num=cumsum([num_panel_conn,num_panel_inboard,num_panel_outboard,num_panel_FWT]);
        
        lift_wing_matrix(:,1:Cum_num(1))=lift_wing_matrix(:,1:Cum_num(1))/width_conn;
        lift_wing_matrix(:,Cum_num(1)+1:Cum_num(2))=lift_wing_matrix(:,Cum_num(1)+1:Cum_num(2))/width_wing1;
        lift_wing_matrix(:,Cum_num(2)+1:Cum_num(3))=lift_wing_matrix(:,Cum_num(2)+1:Cum_num(3))/width_wing2;
        lift_wing_matrix(:,Cum_num(3)+1:Cum_num(4))=lift_wing_matrix(:,Cum_num(3)+1:Cum_num(4))/width_fwt;
        
    else
        lift_wing=[Lift_conn';Lift_wing_inboard';Lift_wing_outboard'];
        
        lift_wing_matrix=reshape(lift_wing, Param.Wing.AeroPanel_Number, numel(lift_wing)/Param.Wing.AeroPanel_Number);
        
        % abs. value of the lift
        lift_wing_abs=sum(lift_wing_matrix);
        
        Cum_num=cumsum([num_panel_conn,num_panel_inboard,num_panel_outboard]);
        
        lift_wing_matrix(:,1:Cum_num(1))=lift_wing_matrix(:,1:Cum_num(1))/width_conn;
        lift_wing_matrix(:,Cum_num(1)+1:Cum_num(2))=lift_wing_matrix(:,Cum_num(1)+1:Cum_num(2))/width_wing1;
        lift_wing_matrix(:,Cum_num(2)+1:Cum_num(3))=lift_wing_matrix(:,Cum_num(2)+1:Cum_num(3))/width_wing2;    
        
    end

    
    % lift per unit span
    lift_wing_var=sum(lift_wing_matrix);
    
    % whole wing span lift distribution
    Y_left=sort(-Y);
    lift_wing_var_left=flip(lift_wing_var);
    lift_wing_abs_left=flip(lift_wing_abs);
    panel_width_left=flip(panel_width);
    
    Y_all=[Y_left,Y]';
    Lift_all=[lift_wing_var_left,lift_wing_var];
    Lift_all_abs=[lift_wing_abs_left,lift_wing_abs];
    
    panel_width_all=[panel_width_left,panel_width];
    
    % calculate vorticity
    % Air conditions
    FlightPoint=awi.model.FlightPoint;
    FlightPoint.Mach=0.78;
    FlightPoint.AcVelocity=0.78*340;
    FlightPoint.Altitude = 36000;
    getFlightPointData(FlightPoint,'ISA');
    DynPressure = FlightPoint.DynPressure;
    
    Gamma_all=Lift_all/(FlightPoint.AirDensity*FlightPoint.AcVelocity);
    
    % calculate the derivative
    dGdy= gradient(Gamma_all(:)) ./ gradient(Y_all(:));
    
    % downwash
    wj=zeros(1,length(Y_all));
    alphai=zeros(1,length(Y_all));
    
    for i=1:length(Y_all)
        
        w=-(1/(4*pi))*dGdy.* panel_width_all'./(Y_all(i)-Y_all);
        
        w=w( ~any( isnan( w ) | isinf( w ), 2 ),: );
        
        wj(i)=sum(w);
        
        alphai(i)=wj(i)/FlightPoint.AcVelocity;
        
    end
    
    Dragi_var=Lift_all_abs.*sin(alphai);
    Drag_force=sum([Dragi_var]);
    
    if ~isempty(FWT_)
        Total_area=(Wingbox_right.SurfaceArea+FWT_R.SurfaceArea)*2 + Param.Wing.Root_Chord*Param.Layout.Fuselage_Width;
    else
        Total_area=(Wingbox_right.SurfaceArea)*2 + Param.Wing.Root_Chord*Param.Layout.Fuselage_Width;
    end
    
    CDi=Drag_force/(DynPressure*Total_area);
    CL=sum(Lift_all_abs)/(DynPressure*Total_area);
    k=CDi/(CL^2);
    
    % Lift & Drag distrbutions 
    
    if ~isempty(FWT_)
        
        Wing_Chords=[Wingbox_right.Chord(1),Wingbox_right.Chord(1),Wingbox_right.Chord(2),FWT_R.Chord(2)];
        
        Eta_Chord=[0, Wingbox_right.YData(1:2) + 0.5*Param.Layout.Fuselage_Width, Param.Wing.Semi_Span + 0.5*Param.Layout.Fuselage_Width];
        
        Chord_Length_R=interp1(Eta_Chord,Wing_Chords,Y);
        
        Chord_Length_L=sort(Chord_Length_R);
        
        Chord_Lengths=[Chord_Length_L, Chord_Length_R];
        
        CL_var=2*Gamma_all./(FlightPoint.AcVelocity * Chord_Lengths);
        
        CDi_var=Dragi_var./(DynPressure.*Chord_Lengths.* panel_width_all);
        
    else
        
        Wing_Chords=[Wingbox_right.Chord(1),Wingbox_right.Chord(1),Wingbox_right.Chord(2),Wingbox_right.Chord(3)];
        
        Eta_Chord=[0, Wingbox_right.YData(1:2) + 0.5*Param.Layout.Fuselage_Width, Param.Wing.Semi_Span + 0.5*Param.Layout.Fuselage_Width];
        
        Chord_Length_R=interp1(Eta_Chord,Wing_Chords,Y);
        
        Chord_Length_L=sort(Chord_Length_R);
        
        Chord_Lengths=[Chord_Length_L, Chord_Length_R];
        
        CL_var=2*Gamma_all./(FlightPoint.AcVelocity * Chord_Lengths);
        
        CDi_var=Dragi_var./(DynPressure.*Chord_Lengths.* panel_width_all);
        
    end
    
    Aerodynamic_distribution.Y=Y_all;
    Aerodynamic_distribution.Y_conn=Y_conn;
    Aerodynamic_distribution.Y_wing1=Y_wing1;
    Aerodynamic_distribution.Y_wing2=Y_wing2;
    
    if ~isempty(FWT_)
        Aerodynamic_distribution.Y_fwt=Y_fwt;
    end
         
    Aerodynamic_distribution.WJ=wj';
    Aerodynamic_distribution.Alphai=alphai';
    Aerodynamic_distribution.Chord_Length=Chord_Lengths';
    Aerodynamic_distribution.CL_var=CL_var';
    Aerodynamic_distribution.CDi_var=CDi_var';
    Aerodynamic_distribution.Lift_var=Lift_all';
    
    Aerodynamic_distribution.Panel_width=panel_width_all;
    
    %% calculate zero-lift drag
    
    CD0=0.017;
    
   
    
    
    
    

end
