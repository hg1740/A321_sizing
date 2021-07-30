function [Static_Loads,Upper_Bound]=Stress_Analysis_ZF(Param, run_folder)

    
 
     
    % $$$---  Fuel removal   ---$$$

    Param.Masses.Fuel_Mass=1;

    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    % create a sub-directory for hinge-locked cases

    if ~exist(strcat(run_folder,'\zero_fuel'),'dir')

        mkdir(run_folder,'zero_fuel')

        run_folder=strcat(run_folder,'\zero_fuel');

    else
        run_folder=strcat(run_folder,'\zero_fuel');
    end

    
    % create Aircraft model
    
    if isfield(Param,'FWT')
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT_R]=Aircraft_Models_v1(Param);
        
    else
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
    end
    
    
    export(FEM_full, run_folder);
       
     %% cruising, 36000 ft, 2.5 g
     
    TrimLoadcase1 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;
    
    TrimLoadcase1.Name = 'Cruise_36000ft_2p5g_zero_fuel';
    
    TrimLoadcase1.Altitude   = altitude;
    TrimLoadcase1.Mach       = mach_number;
    TrimLoadcase1.AcVelocity = aircraft_velocity;
    TrimLoadcase1.AcMass = acMass;
    
    TrimLoadcase1.PitchAngle=0;
    TrimLoadcase1.RollAngle =0;
    TrimLoadcase1.ID = 1030;
    TrimLoadcase1.LoadFactor = 2.5;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase1.CsDeflection=flap_angle*pi/180;
    
    
    %% dive, 3000 ft, - g
    
    TrimLoadcase2 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 3000;
    mach_number       = 0.48;
    aircraft_velocity = 0.48*340;
    flap_angle=0;
    
    TrimLoadcase2.Name = 'Dive_3000ft_g_zero_fuel';
    
    TrimLoadcase2.Altitude   = altitude;
    TrimLoadcase2.Mach       = mach_number;
    TrimLoadcase2.AcVelocity = aircraft_velocity;
    TrimLoadcase2.AcMass = acMass;
    
    TrimLoadcase2.PitchAngle=0;
    TrimLoadcase2.RollAngle =0;
    TrimLoadcase2.ID = 1030;
    TrimLoadcase2.LoadFactor = -1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase2.CsDeflection=flap_angle*pi/180;
    

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
    
    %Write trim file 1
    trimFile1 = NastranMethods.writeTrimFile(Aircraft, TrimLoadcase1,...
        MassCases,run_folder,'DatFilename','Cruise_36000ft_2p5g_zero_fuel',...
        'ExtraCards',aelink);
    
    %Write trim file 2
    trimFile2 = NastranMethods.writeTrimFile(Aircraft, TrimLoadcase2,...
        MassCases,run_folder,'DatFilename','Dive_3000ft_g_zero_fuel',...
        'ExtraCards',aelink);
 
    
    %% Run analysis
    
    % clear directory
    
    delete(strcat(run_folder, '\*.xdb'));
    delete(strcat(run_folder, '\*.h5'));
    delete(strcat(run_folder, '\*.log'));
    delete(strcat(run_folder, '\*.f06'));
    delete(strcat(run_folder, '\*.f04'));
    delete(strcat(run_folder, '\*.op4'));
    
    NastranMethods.runNastran(trimFile1);
    NastranMethods.runNastran(trimFile2);

    
    %% Load 
    
    Trim1_res=NastranMethods.extractNastranResults(strcat(run_folder,'\','Cruise_36000ft_2p5g_zero_fuel','.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim2_res=NastranMethods.extractNastranResults(strcat(run_folder,'\','Dive_3000ft_g_zero_fuel','.f06'),'ReadF06',true,'ReadHDF5',false);
    
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
    
    % result data Trim 1: ZF pull up
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
    
    
    % result data Trim 2: ZF dive
    Trim2.M_P1=[Trim2_res.f06data.Bendingmoment.UMPLN1(index1),Trim2_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Trim2.M_P2=[Trim2_res.f06data.Bendingmoment.UMPLN2(index1),Trim2_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Trim2.T=[Trim2_res.f06data.Bendingmoment.UTORQUE1(index1),Trim2_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
    
    Trim2.S_P1=[Trim2_res.f06data.Bendingmoment.USPLN1(index1),Trim2_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Trim2.S_P2=[Trim2_res.f06data.Bendingmoment.USPLN2(index1),Trim2_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    Moment_P1_Loadcase2=Trim2.M_P1; % internal load M1 for 1.5 g pull up
    Moment_P2_Loadcase2=Trim2.M_P2; % internal load M2 for 1.5 g pull up
    Torque_Loadcase2=Trim2.T;       % internal load T for 1.5 g pull up
    Shear_P1_Loadcase2=Trim2.S_P1;
    Shear_P2_Loadcase2=Trim2.S_P2;
    
   
    % find upper bound for sizing
    
    % Moment
    Sizing_Moment_P2=max([Moment_P2_Loadcase1;Moment_P2_Loadcase2]);
    
    % Shear
    Sizing_Shear_P2=max([Shear_P2_Loadcase1;Shear_P2_Loadcase2]);
    
    % Torque
    Sizing_Torque=max([Torque_Loadcase1;Torque_Loadcase2]);
    

    %%  Result output - static loads
    Static_Loads.Y=Y_all;
    Static_Loads.ZF_pullup.Moment_P2=Moment_P2_Loadcase1';
    Static_Loads.ZF_pullup.Shear_P2=Shear_P2_Loadcase1';
    Static_Loads.ZF_pullup.Torque=Torque_Loadcase1';
    
    Static_Loads.ZF_dive.Moment_P2=Moment_P2_Loadcase2';
    Static_Loads.ZF_dive.Shear_P2=Shear_P2_Loadcase2';
    Static_Loads.ZF_dive.Torque=Torque_Loadcase2';
    

    % Result output - upper bound
    Upper_Bound.Moment_P2=Sizing_Moment_P2';
    Upper_Bound.Shear_P2=Sizing_Shear_P2';
    Upper_Bound.Torque=Sizing_Torque';
    

end
