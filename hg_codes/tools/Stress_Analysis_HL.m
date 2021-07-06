function [Static_Loads,Delta,Upper_Bound]=Stress_Analysis_HL(Param, run_folder)

    
    if  isfield(Param,'FWT')
     
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
    
    % create Aircraft model
    
    if isfield(Param,'FWT')
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT_R]=Aircraft_Models_v1(Param);
        
    else
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
    end
    
    
    export(FEM_full, run_folder);
       
     %% cruising, 36000 ft, g
    TrimLoadcase = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;
    
    TrimLoadcase.Name = 'Cruise_3000ft_1g_hinge_locked';
    
    TrimLoadcase.Altitude   = altitude;
    TrimLoadcase.Mach       = mach_number;
    TrimLoadcase.AcVelocity = aircraft_velocity;
    TrimLoadcase.AcMass = acMass;
    
    TrimLoadcase.PitchAngle=0;
    TrimLoadcase.RollAngle =0;
    TrimLoadcase.ID = 1030;
    TrimLoadcase.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase.CsDeflection=flap_angle*pi/180;
    
    
    %% Gust threshold 25% of gust velocity 

    GustLoadcase = awi.model.LoadCase;
    GustLoadcase.Altitude   = 3000;
    GustLoadcase3.AcVelocity = 0.48*340;
    GustLoadcase3.AcMass = 500;
    GustLoadcase3.Mach = 0.48;
    GustLoadcase3.GustLength = linspace(18,214,7);
    
    % Gust direction: positive or negative hit
    GustLoadcase.GustDirection=0.25;
    
    FlightPoint=awi.model.FlightPoint;
    FlightPoint.Mach=0.48;
    FlightPoint.AcVelocity=FlightPoint.Mach*340;
    FlightPoint.Altitude = 3000;
    getFlightPointData(FlightPoint,'ISA');
    
      
    %% Write input
    
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
        MassCases,run_folder,'DatFilename','Cruise_3000ft_1g_hinge_locked',...
        'ExtraCards',aelink);

    
    % Write gust
    gustfile=NastranMethods1.writeGustFile(Aircraft, GustLoadcase, MassCases, FlightPoint, run_folder,'DatFilename','gust_threshold_3000ft_hinge_locked');
    
    
    %% Run analysis
    
    % clear directory
    
    delete(strcat(run_folder, '\*.xdb'));
    delete(strcat(run_folder, '\*.h5'));
    delete(strcat(run_folder, '\*.log'));
    delete(strcat(run_folder, '\*.f06'));
    delete(strcat(run_folder, '\*.f04'));
    delete(strcat(run_folder, '\*.op4'));
    
    NastranMethods.runNastran(trimFile);
    NastranMethods1.runNastran(gustfile);
    
    %% Load 
    
    Trim1_res=NastranMethods.extractNastranResults(strcat(run_folder,'\','Cruise_3000ft_1g_hinge_locked','.f06'),'ReadF06',true,'ReadHDF5',false);
    
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
    
    Moment_P1_Loadcase1=Trim1.M_P1; % internal load M1 for 1.5 g pull up
    Moment_P2_Loadcase1=Trim1.M_P2; % internal load M2 for 1.5 g pull up
    Torque_Loadcase1=Trim1.T;       % internal load T for 1.5 g pull up
    Shear_P1_Loadcase1=Trim1.S_P1;
    Shear_P2_Loadcase1=Trim1.S_P2;
    
    % correction stresses at hinge
    Moment_P2_Loadcase1(end-10)=(Moment_P2_Loadcase1(end-9)+Moment_P2_Loadcase1(end-11))/2;
    Shear_P2_Loadcase1(end-10)=(Shear_P2_Loadcase1(end-9)+Shear_P2_Loadcase1(end-11))/2;
    Torque_Loadcase1(end-10)=(Torque_Loadcase1(end-9)+Torque_Loadcase1(end-11))/2;
    

    % extract gust result: Max-Min
    NumSteps=201;
    [Root_Delta, Wing_Delta]=Gust_peaks(All_nodes,GustLoadcase,run_folder,'\gust_threshold_3000ft_hinge_locked.h5',NumSteps);
    
    % correction incremental stresses at hinge
    Wing_Delta.Max_Moment(end-9)=(Wing_Delta.Max_Moment(end-8) + Wing_Delta.Max_Moment(end-10))/2;
    Wing_Delta.Min_Moment(end-9)=(Wing_Delta.Min_Moment(end-8) + Wing_Delta.Min_Moment(end-10))/2;
    
    Wing_Delta.Max_Shear(end-9)=(Wing_Delta.Max_Shear(end-8) + Wing_Delta.Max_Shear(end-10))/2;
    Wing_Delta.Min_Shear(end-9)=(Wing_Delta.Min_Shear(end-8) + Wing_Delta.Min_Shear(end-10))/2;
    
    Wing_Delta.Max_Torque(end-9)=(Wing_Delta.Max_Torque(end-8) + Wing_Delta.Max_Torque(end-10))/2;
    Wing_Delta.Min_Torque(end-9)=(Wing_Delta.Min_Torque(end-8) + Wing_Delta.Min_Torque(end-10))/2;
    
 
    % for sizing - upper bound
    Moment_P2_pull_up= 1.5*Moment_P2_Loadcase1;
    Moment_P2_gust=Moment_P2_Loadcase1 + [Wing_Delta.Max_Moment; 0]';
    Sizing_Moment_P2=max([Moment_P2_pull_up;Moment_P2_gust]);
    
    Shear_P2_pull_up = 1.5*Shear_P2_Loadcase1;
    Shear_P2_gust=Shear_P2_Loadcase1 + [Wing_Delta.Max_Shear; 0]';
    Sizing_Shear_P2=max([Shear_P2_pull_up;Shear_P2_gust]);
    
    Torque_pull_up = 1.5*Torque_Loadcase1;
    Torque_gust=Torque_Loadcase1 + [Wing_Delta.Max_Torque; 0]';
    Sizing_Torque=max([Torque_pull_up;Torque_gust]);
    
    % Result output - static loads
    Static_Loads.Y=Y_all;
    Static_Loads.Moment_P2=Moment_P2_Loadcase1';
    Static_Loads.Shear_P2=Shear_P2_Loadcase1';
    Static_Loads.Torque=Torque_Loadcase1';
    
    % Result output - delta
    Delta.Wing= Wing_Delta;
    Delta.Root= Root_Delta;
    
    % Result output - upper bound
    Upper_Bound.Moment_P2=Sizing_Moment_P2';
    Upper_Bound.Shear_P2=Sizing_Shear_P2';
    Upper_Bound.Torque=Sizing_Torque';
    

end
