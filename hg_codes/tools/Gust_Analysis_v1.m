function Loads=Gust_Analysis_v1(Param,run_folder,varargin)


    p=inputParser;

    addParameter(p,'File_Name','file_name',@(x)validateattributes(x,{'char'},{'nonempty'}))
    addParameter(p,'Hinge_Lock','on',@(x)any(validatestring(x,{'on','off'})))
    addParameter(p,'Altitude',36000,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))
    addParameter(p,'March_Number',0.78,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))
    addParameter(p,'Gust_Eta',1,@(x)validateattributes(x,{'numeric'},{'scalar', 'nonnan', 'finite', 'real','nonempty'}))

    parse(p,varargin{:})

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
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT_R]=Aircraft_Models_v1(Param);

    else
        [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
    end

    export(FEM_full, run_folder);
    
    
    %% set-up loadcase 
    
    GustLoadcase = awi.model.LoadCase;
    GustLoadcase.Altitude   = p.Results.Altitude;
    GustLoadcase.Mach = p.Results.March_Number;
    GustLoadcase.AcVelocity = GustLoadcase.Mach*340;
    GustLoadcase.AcMass = 500;
    GustLoadcase.GustLength = linspace(18,214,7);%[18:20:214];
    
    % Gust direction: positive or negative hit
    GustLoadcase.GustDirection=p.Results.Gust_Eta;
    
    FlightPoint=awi.model.FlightPoint;
    FlightPoint.Mach=p.Results.March_Number;
    FlightPoint.AcVelocity=FlightPoint.Mach*340;
    FlightPoint.Altitude = p.Results.Altitude;
    getFlightPointData(FlightPoint,'ISA');
    
    
    %% Run analysis
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_full;
    MassCases=awi.model.MassCases.empty;
    
    gustfile=NastranMethods1.writeGustFile(Aircraft, GustLoadcase, MassCases, FlightPoint, run_folder,'DatFilename',p.Results.File_Name);
    
    delete(strcat(run_folder, '\*.xdb'));
    delete(strcat(run_folder, '\*.h5'));
    delete(strcat(run_folder, '\*.log'));
    delete(strcat(run_folder, '\*.f06'));
    delete(strcat(run_folder, '\*.f04'));
    delete(strcat(run_folder, '\*.op4'));
    
    delete(strcat(run_folder, '\*.xdb.*'));
    delete(strcat(run_folder, '\*.h5.*'));
    delete(strcat(run_folder, '\*.log.*'));
    delete(strcat(run_folder, '\*.f06.*'));
    delete(strcat(run_folder, '\*.f04.*'));
    
    NastranMethods1.runNastran(gustfile);
    
    
    %% Result extract
    
    WingNodes=NastranMethods1.WingNode;
    FWTNodes=NastranMethods1.FWTNode;
    
    if ~isempty(FWTNodes)
        
        All_nodes=[WingNodes(1:end-1),FWTNodes(1:end)];
        
    else
        
        All_nodes=[WingNodes(1:end)];
        
    end
    
    X_all=[All_nodes.X];
    Y_all=X_all(2,:)';
    
    NumSteps=201;
    
    Result_file=strcat('\',p.Results.File_Name,'.h5');
    
    % extract peaks
    [Root_Delta, Wing_Delta]=Gust_peaks(All_nodes,GustLoadcase,run_folder,Result_file,NumSteps);
    
    % time history at root
    [Moment_Root, Torque_Root, Shear_Root]=Gust_Thistory(WingNodes,GustLoadcase,run_folder,Result_file,NumSteps);
    
   % write output
    Loads.Y=Y_all;
    Loads.Wing_Delta=Wing_Delta;
    Loads.Root_Delta=Root_Delta;
    
    Loads.Time_Response.Root_Moment=Moment_Root;
    Loads.Time_Response.Root_Torque=Torque_Root;
    Loads.Time_Response.Root_Shear=Shear_Root;
    
    


end

