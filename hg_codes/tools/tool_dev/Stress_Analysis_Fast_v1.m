function [Aircraft, Wingbox_right, FEM_full, Y_all, Box_dimensions, Box_CrossSec, Load_distribution, Wing_Delta, Root_Delta, Internal_Stresses]=Stress_Analysis_Fast_v1(Param, run_folder)
    

    % update Panel aspect ratio 
    
    Param.Wing.AeroPanel_AR=2.5;


    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec]=Aircraft_Models_v1(Param);

    export(FEM_full, run_folder);

    %% LoadCase 1: cruising altitude 2.5g pull up
    TrimLoadcase1 = awi.model.LoadCase;

    acMass = 94000;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;
    

    TrimLoadcase1.Name = 'A321_cruise_2p5g';
    TrimLoadcase1.Altitude   = altitude;
    TrimLoadcase1.Mach       = mach_number;
    TrimLoadcase1.AcVelocity = aircraft_velocity;
    TrimLoadcase1.AcMass = acMass;
    

    TrimLoadcase1.PitchAngle=0;
    TrimLoadcase1.RollAngle =0;
    TrimLoadcase1.ID = 1020;
    TrimLoadcase1.LoadFactor = 2.5;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase1.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase1)
    
    
    %% Loadcase 2: sea level 2.5g pull up
    TrimLoadcase2 = awi.model.LoadCase;

    acMass = 500;
    altitude          = 3000;
    mach_number       = 0.48;
    aircraft_velocity = 0.48*340;
    flap_angle=0;

    TrimLoadcase2.Name = 'A321_sea_level_2.5g';

    TrimLoadcase2.Altitude   = altitude;
    TrimLoadcase2.Mach       = mach_number;
    TrimLoadcase2.AcVelocity = aircraft_velocity;
    TrimLoadcase2.AcMass = acMass;

    TrimLoadcase2.PitchAngle=0;
    TrimLoadcase2.RollAngle =0;
    TrimLoadcase2.ID = 1030;
    TrimLoadcase2.LoadFactor = 2.5;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase2.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase2)
    
    %% Loadcase 3: sea level - g sea level dive
    TrimLoadcase3 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 3000;
    mach_number       = 0.48;
    aircraft_velocity = 0.48*340;
    flap_angle=0;

    TrimLoadcase3.Name = 'A321_sea_level_neg_g';

    TrimLoadcase3.Altitude   = altitude;
    TrimLoadcase3.Mach       = mach_number;
    TrimLoadcase3.AcVelocity = aircraft_velocity;
    TrimLoadcase3.AcMass = acMass;

    TrimLoadcase3.PitchAngle=0;
    TrimLoadcase3.RollAngle =0;
    TrimLoadcase3.ID = 1030;
    TrimLoadcase3.LoadFactor = -1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase3.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase3)
    
    
    %% steady state for gust analysis, 3000 ft, g
    TrimLoadcase6 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 3000;
    mach_number       = 0.48;
    aircraft_velocity = 0.48*340;
    flap_angle=0;
    
    TrimLoadcase6.Name = 'A321_3000ft_g';
    
    TrimLoadcase6.Altitude   = altitude;
    TrimLoadcase6.Mach       = mach_number;
    TrimLoadcase6.AcVelocity = aircraft_velocity;
    TrimLoadcase6.AcMass = acMass;
    
    TrimLoadcase6.PitchAngle=0;
    TrimLoadcase6.RollAngle =0;
    TrimLoadcase6.ID = 1030;
    TrimLoadcase6.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase6.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase6)


    %% Gust analysis Case 3 : sea level,  + 1MC
    
    GustLoadcase3 = awi.model.LoadCase;
    %     GustLoadcase.LoadCaseType = 'Pratt Gust';
    GustLoadcase3.Altitude   = 3000;
    GustLoadcase3.AcVelocity = 0.48*340;
    GustLoadcase3.AcMass = 500;
    GustLoadcase3.Mach = 0.48;
    GustLoadcase3.GustLength = linspace(18,214,7);%[18:20:214];
    
    % Gust direction: positive or negative hit
    GustLoadcase3.GustDirection=1;
    
    FlightPoint3=awi.model.FlightPoint;
    FlightPoint3.Mach=0.48;
    FlightPoint3.AcVelocity=FlightPoint3.Mach*340;
    FlightPoint3.Altitude = 3000;
    getFlightPointData(FlightPoint3,'ISA');
    

    %% write inputs 
    
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_full;
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

    
    % sizing case 1
    trimFile1 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase1, MassCases,run_folder,'DatFilename','A321_cruise_2p5g','ExtraCards',aelink);
    
    % sizing case 2
    trimFile2 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase2, MassCases,run_folder,'DatFilename','A321_sealevel_2p5g','ExtraCards',aelink);
    
    % sizing case 3
    trimFile3 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase3, MassCases,run_folder,'DatFilename','A321_sealevel_neg_1g','ExtraCards',aelink);
    
    % sizing case 6
    trimFile6 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase6, MassCases,run_folder,'DatFilename','A321_3000ft_1g','ExtraCards',aelink);

    % gust 3
    gustfile3=NastranMethods1.writeGustFile(Aircraft, GustLoadcase3, MassCases, FlightPoint3, run_folder,'DatFilename','gust_analysis_g_3000ft_pos_1MC');
    

    %%  run analysis
    
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

    
    NastranMethods1.runNastran(trimFile1);
    NastranMethods1.runNastran(trimFile2);
    NastranMethods1.runNastran(trimFile3);

    NastranMethods1.runNastran(trimFile6);
    
    
    NastranMethods1.runNastran(gustfile3);
    
    
    %%  ***  Hinge lock cases  ***
    

    if isfield(Param,'FWT')
        
        [HL_Static,HL_Delta,HL_Upper_Bound]=Stress_Analysis_HL(Param, run_folder);

    else
        HL_Upper_Bound.Moment_P2=zeros(25,1);
        HL_Upper_Bound.Shear_P2=zeros(25,1);
        HL_Upper_Bound.Torque=zeros(25,1);

    end
        

    %%  ***  Zero fuel cases  ***
    
    
    
    
    [ZF_Static_Loads,ZF_Upper_Bound]=Stress_Analysis_ZF(Param, run_folder);
    
    
    
    

    %% Result extraction 
    

    Trim1_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_cruise_2p5g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim2_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_sealevel_2p5g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim3_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_sealevel_neg_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    
    Trim6_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_3000ft_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    
    
    WingNodes=NastranMethods1.WingNode;
    FWTNodes=NastranMethods1.FWTNode;
    
    if ~isempty(FWTNodes)
        
        All_nodes=[WingNodes(1:end-1),FWTNodes(1:end)];
        
    else
        
        All_nodes=[WingNodes(1:end)];
        
    end
    
    X_all=[All_nodes.X];
    Y_all=X_all(2,:)';
    
    index1=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[All_nodes(1:end-1).GID]);
    index2=ismember([Trim1_res.f06data.Bendingmoment.LGRID],[All_nodes(end).GID]);
    
    % extract static results:
    
    % forces from trim 1: 2.5 g 36000ft
    Trim1.M_P1=[Trim1_res.f06data.Bendingmoment.UMPLN1(index1),Trim1_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Trim1.M_P2=[Trim1_res.f06data.Bendingmoment.UMPLN2(index1),Trim1_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Trim1.T=[Trim1_res.f06data.Bendingmoment.UTORQUE1(index1),Trim1_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque

    Trim1.S_P1=[Trim1_res.f06data.Bendingmoment.USPLN1(index1),Trim1_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Trim1.S_P2=[Trim1_res.f06data.Bendingmoment.USPLN2(index1),Trim1_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    % forces from trim 2: 2.5g Sea level
    Trim2.M_P1=[Trim2_res.f06data.Bendingmoment.UMPLN1(index1),Trim2_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Trim2.M_P2=[Trim2_res.f06data.Bendingmoment.UMPLN2(index1),Trim2_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Trim2.T=[Trim2_res.f06data.Bendingmoment.UTORQUE1(index1),Trim2_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque

    Trim2.S_P1=[Trim2_res.f06data.Bendingmoment.USPLN1(index1),Trim2_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Trim2.S_P2=[Trim2_res.f06data.Bendingmoment.USPLN2(index1),Trim2_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    % forces from trim 3: -g sea level
    Trim3.M_P1=[Trim3_res.f06data.Bendingmoment.UMPLN1(index1),Trim3_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Trim3.M_P2=[Trim3_res.f06data.Bendingmoment.UMPLN2(index1),Trim3_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Trim3.T=[Trim3_res.f06data.Bendingmoment.UTORQUE1(index1),Trim3_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
    
    Trim3.S_P1=[Trim3_res.f06data.Bendingmoment.USPLN1(index1),Trim3_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Trim3.S_P2=[Trim3_res.f06data.Bendingmoment.USPLN2(index1),Trim3_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    
    
    % forces from trim 6: g 3000ft
    Trim6.M_P1=[Trim6_res.f06data.Bendingmoment.UMPLN1(index1),Trim6_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Trim6.M_P2=[Trim6_res.f06data.Bendingmoment.UMPLN2(index1),Trim6_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Trim6.T=[Trim6_res.f06data.Bendingmoment.UTORQUE1(index1),Trim6_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
    
    Trim6.S_P1=[Trim6_res.f06data.Bendingmoment.USPLN1(index1),Trim6_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Trim6.S_P2=[Trim6_res.f06data.Bendingmoment.USPLN2(index1),Trim6_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
       
    
    % extract gust result: Max-Min
    NumSteps=201;
    
    [Root_Delta_LC3, Wing_Delta_LC3]=Gust_peaks(All_nodes,GustLoadcase3,run_folder,'\gust_analysis_g_3000ft_pos_1MC.h5',NumSteps);
    
    
    %% Internal Load Summery 
    
    % Load case 1: 2.5 g pull up crusing altutude 
    Moment_P1_Loadcase1=abs(Trim1.M_P1); % internal load M1 for 2.5 g pull up
    Moment_P2_Loadcase1=abs(Trim1.M_P2); % internal load M2 for 2.5 g pull up
    Torque_Loadcase1=abs(Trim1.T);       % internal load T for 2.5 g pull up
    Shear_P1_Loadcase1=abs(Trim1.S_P1);
    Shear_P2_Loadcase1=abs(Trim1.S_P2);
    
    % Load case 2: 2.5 g pull up sea level 
    Moment_P1_Loadcase2=abs(Trim2.M_P1); % internal load M1 for 2.5 g pull up
    Moment_P2_Loadcase2=abs(Trim2.M_P2); % internal load M2 for 2.5 g pull up
    Torque_Loadcase2=abs(Trim2.T);       % internal load T for 2.5 g pull up
    Shear_P1_Loadcase2=abs(Trim2.S_P1);
    Shear_P2_Loadcase2=abs(Trim2.S_P2);
    
    % Load case 3: - g dive at sea level
    Moment_P1_Loadcase3=abs(Trim3.M_P1); % internal load M1 for 2.5 g pull up
    Moment_P2_Loadcase3=abs(Trim3.M_P2); % internal load M2 for 2.5 g pull up
    Torque_Loadcase3=abs(Trim3.T);       % internal load T for 2.5 g pull up
    Shear_P1_Loadcase3=abs(Trim3.S_P1);
    Shear_P2_Loadcase3=abs(Trim3.S_P2);
    
    
    
    % Load case 6: 3000ft g and + 1MC
    Moment_P1_Loadcase6=abs(Trim6.M_P1);
    Moment_P2_Loadcase6=abs(Trim6.M_P2) + [Wing_Delta_LC3.Max_Moment; 0]';
    Torque_Loadcase6=abs(Trim6.T) + [Wing_Delta_LC3.Max_Torque; 0]';
    Shear_P1_Loadcase6=abs(Trim6.S_P1);
    Shear_P2_Loadcase6=abs(Trim6.S_P2) + [Wing_Delta_LC3.Max_Shear; 0]';
    
    
    
    % Load case 9: 3000ft g and - 1MC
    Moment_P1_Loadcase9=abs(Trim6.M_P1);
    Moment_P2_Loadcase9=abs(Trim6.M_P2) + abs([Wing_Delta_LC3.Min_Moment; 0]');
    Torque_Loadcase9=abs(Trim6.T) + abs([Wing_Delta_LC3.Min_Torque; 0]');
    Shear_P1_Loadcase9=abs(Trim6.S_P1);
    Shear_P2_Loadcase9=abs(Trim6.S_P2) + abs([Wing_Delta_LC3.Min_Shear; 0]');
    
    
    % hinge locked cases
    
    Moment_P2_HingeLock=HL_Upper_Bound.Moment_P2;
    
    Shear_P2_HingeLock=HL_Upper_Bound.Shear_P2;
    
    Torque_HingeLock=HL_Upper_Bound.Torque;
    
    % Maximum forces used for sizing
    M_P1=max([Moment_P1_Loadcase1; Moment_P1_Loadcase2; Moment_P1_Loadcase3;  Moment_P1_Loadcase6; Moment_P1_Loadcase9]);
    M_P2=max([Moment_P2_Loadcase1; Moment_P2_Loadcase2; Moment_P2_Loadcase3;  Moment_P2_Loadcase6; Moment_P2_Loadcase9; abs(Moment_P2_HingeLock)'; ZF_Upper_Bound.Moment_P2']);
    T=max([Torque_Loadcase1; Torque_Loadcase2; Torque_Loadcase3; Torque_Loadcase6; Torque_Loadcase9; abs(Torque_HingeLock'); abs(ZF_Upper_Bound.Torque)']);
    
    S_P1=max([Shear_P1_Loadcase1; Shear_P1_Loadcase2; Shear_P1_Loadcase3; Shear_P1_Loadcase6; Shear_P1_Loadcase9]);
    S_P2=max([Shear_P2_Loadcase1; Shear_P2_Loadcase2; Shear_P2_Loadcase3; Shear_P2_Loadcase6; Shear_P2_Loadcase9; abs(Shear_P2_HingeLock'); ZF_Upper_Bound.Shear_P2']);
        
  
    %% Result output
    
    % Out of plane bending moment distributions for all load cases
    Load_distribution.Moment_P2.LC1=Moment_P2_Loadcase1';
    Load_distribution.Moment_P2.LC2=Moment_P2_Loadcase2';
    Load_distribution.Moment_P2.LC3=Moment_P2_Loadcase3';
    
    Load_distribution.Moment_P2.LC6=Moment_P2_Loadcase6';

    Load_distribution.Moment_P2.LC9=Moment_P2_Loadcase9';
    
    Load_distribution.Moment_P2.HL=Moment_P2_HingeLock;
    
    Load_distribution.Moment_P2.ZF=ZF_Upper_Bound.Moment_P2;
    
    % Vertical shear force distributions for all load cases
    Load_distribution.Shear_P2.LC1=Shear_P2_Loadcase1';
    Load_distribution.Shear_P2.LC2=Shear_P2_Loadcase2';
    Load_distribution.Shear_P2.LC3=Shear_P2_Loadcase3';
   
    Load_distribution.Shear_P2.LC6=Shear_P2_Loadcase6';
    
    Load_distribution.Shear_P2.LC9=Shear_P2_Loadcase9';
    
    Load_distribution.Shear_P2.HL=Shear_P2_HingeLock;
    
    Load_distribution.Shear_P2.ZF=ZF_Upper_Bound.Shear_P2;
    
    % Torque distribution for all load cases 
    Load_distribution.Torque.LC1=(Trim1.T)';
    Load_distribution.Torque.LC2=(Trim2.T)';
    Load_distribution.Torque.LC3=(Trim3.T)';
    
    Load_distribution.Torque.LC6=(Trim6.T + [Wing_Delta_LC3.Max_Torque; 0]')';
   
    Load_distribution.Torque.LC9=(Trim6.T + abs([Wing_Delta_LC3.Min_Torque; 0]'))';
    
    Load_distribution.Torque.HL=Torque_HingeLock;
    
    Load_distribution.Torque.ZF=ZF_Upper_Bound.Torque;
    
    
    % HL & ZF
    
    if isfield(Param,'FWT')
        
        Load_distribution.HL=HL_Static;
        
        Wing_Delta.HL=HL_Delta.Wing_Delta;
        
        Root_Delta.HL=HL_Delta.Root_Delta;
    
    end
    
    Load_distribution.ZF=ZF_Static_Loads;
    
    
    
    % Delta moment
    
    Wing_Delta.Moment_Max.LC3 = Wing_Delta_LC3.Max_Moment;
    

    Wing_Delta.Moment_Min.LC3 = Wing_Delta_LC3.Min_Moment;
    
    % Delta shear
    
    Wing_Delta.Shear_Max.LC3=Wing_Delta_LC3.Max_Shear;
    
    Wing_Delta.Shear_Min.LC3=Wing_Delta_LC3.Min_Shear;
    
    % Delta torque
    
    Wing_Delta.Torque_Max.LC3=Wing_Delta_LC3.Max_Torque;
    

    Wing_Delta.Torque_Min.LC3=Wing_Delta_LC3.Min_Torque;
   
    
   
    
    
    % Root Delta moment
    Root_Delta.Moment_Max.LC3 = Root_Delta_LC3.Max_Moment';
    

    Root_Delta.Moment_Min.LC3 = Root_Delta_LC3.Min_Moment';
    
    % Root Delta shear

    Root_Delta.Shear_Max.LC3=Root_Delta_LC3.Max_Shear';
    

    Root_Delta.Shear_Min.LC3=Root_Delta_LC3.Min_Shear';
    
    % Root Delta torque

    Root_Delta.Torque_Max.LC3=Root_Delta_LC3.Max_Torque';
    

    Root_Delta.Torque_Min.LC3=Root_Delta_LC3.Min_Torque';
    
   
    
    
    %% calculate stresses
    
    
    if isfield(Param,'FWT')
        
        Spar_Thickness=Param.Wing.Thickness(1:25);
        Skin_Thickness=Param.Wing.Thickness(26:50);
        Stringer_Area=Param.Wing.Thickness(51:75);
        
        FWT_Spar=Param.FWT.Thickness(1:11);
        FWT_Skin=Param.FWT.Thickness(12:22);
        FWT_Stringer=Param.FWT.Thickness(23:33);
        
        % update
        Spar_Thickness=[Spar_Thickness, FWT_Spar(2:end)];
        Skin_Thickness=[Skin_Thickness, FWT_Skin(2:end)];
        Stringer_Area=[Stringer_Area,FWT_Stringer(2:end)];
        
        Stringer_length=sqrt(Stringer_Area/0.36);
        Stringer_thickness=0.12*Stringer_length;
        
        NumSec=35;
        
    else
        
        Spar_Thickness=Param.Wing.Thickness(1:25);
        Skin_Thickness=Param.Wing.Thickness(26:50);
        Stringer_Area=Param.Wing.Thickness(51:75);
        
        Stringer_length=sqrt(Stringer_Area/0.36);
        Stringer_thickness=0.12*Stringer_length;
        
        NumSec=25;
        
    end
    
    E=70e9;

    sigma_skn=ones(1,NumSec);
    tau_skin=ones(1,NumSec);

    sigma_spr=ones(1,NumSec); 
    tau_spr=ones(1,NumSec);

    Von_skn=ones(1,NumSec);
    Von_spr=ones(1,NumSec);

    sigmab_skn=ones(1,NumSec);
    taub_skn=ones(1,NumSec);

    sigma_strg=ones(1,NumSec);
    sigma_crip=ones(1,NumSec);
    sigma_col=ones(1,NumSec);
    sigma_pp=ones(1,NumSec);

    for jj=1:NumSec

        % geometry info.
        ts=Stringer_thickness(jj);
        ds=Stringer_length(jj);
        strg_n=Param.Wing.StringerPitch;

        hs=Box_dimensions.Inboard.Height(jj);
        ws=Box_dimensions.Inboard.Width(jj);
        
        Ixx=Box_CrossSec.Ixx(jj);
        Izz=Box_CrossSec.Izz(jj);
        
        t1=Spar_Thickness(jj); % spar
        t2=Skin_Thickness(jj); %skin

        % load info
        m1=abs(M_P1(jj));
        m2=abs(M_P2(jj));
        tc=abs(T(jj));
        sp1=abs(S_P1(jj));
        sp2=abs(S_P2(jj));

        % Calc.
        sigma_skn(jj)=0.5*m2*hs/Ixx;
        tau_skin(jj)= tc/(2*hs*ws*t2)+ sp2/(2*hs*t1+2*ws*t2);
%         tau_skin(jj)=0;

        sigma_spr(jj)=0.5*m1*ws/Izz;
        % old shear
%         tau_spr(jj)=sp2/(2*hs*t1) + abs(tc)/(2*hs*ws*t1);
        Q=ws*t2*hs/2 + hs^2*t1/4;
       
        V=sp2/(2*hs*t1+2*ws*t2);
        % Corrected vertical shear!!!!
%         V=sp2;
       
        tau_spr(jj)=V*Q/(2*Ixx*t1) + tc/(2*hs*ws*t1);

        % Von mises: skin and spar
        Von_skn(jj)=sqrt(sigma_skn(jj)^2);
        Von_spr(jj)=sqrt(sigma_spr(jj)^2+3*tau_spr(jj)^2);

        % stringers compressive stress

        sigma_strg(jj)=m2*(0.5*hs-0.5*ds)/Ixx;

        % critical buckling stress - for skin--------------------------
        sigmab_skn(jj)=4*pi^2*E*(t2/strg_n)^2/(12*(1-0.33^2));
        taub_skn(jj)=5.6*pi^2*E*(t2/strg_n)^2/(12*(1-0.33^2));      
        sigma_pp(jj) = 0.5*sigma_skn(jj)+sqrt((0.5*sigma_skn(jj))^2+tau_skin(jj)^2);

        % yield stress - skin and spar tensile failure
        sigma_Y=5.2E8;

        % striner crippling stress--------------------------------------
        A1=-0.7885; B1=0.6194; A2=-0.8046; B2=1.2117;

        Rc1_=B1*(sqrt(sigma_Y/E)*(ds/ts))^A1;
        Rc2_=B2*(sqrt(sigma_Y/E)*(ds/ts))^A2;

        sigmac1_=min(Rc1_,1.45)*sigma_Y;
        sigmac2_=min(Rc2_,1.45)*sigma_Y;

        sigma_crip(jj)=(ds*ts*sigmac1_*2 + ds*ts*sigmac2_)/(ts*ds*3);
        %-------------------------------------------------------------

        % column buckling

%         w_eff=1.7*t2*sqrt(E/abs(sigma_strg(jj)));
%         A_eff=w_eff*t2;
        L=0.65; %take as the rib pitch 0.65 m 
        A=3*ds*ts;

        I=(ts*ds^3/12)+(ds*ts^3/12 + ds*ts*(ds/2)^2)*2;
        rg=sqrt(I/A);
%         rg=sqrt(I/(ds*ts*3+A_eff));
%         K=L/(rg*1);
        
        K=L/(rg*2);

        Kcr=pi*sqrt(2*E/sigma_crip(1));

        if K < Kcr      
            sigma_col(jj)=sigma_crip(jj)*(1-sigma_crip(jj)*K^2/(4*pi^2*E));       
        elseif K > Kcr
            sigma_col(jj)=pi^2*E/K^2;
            
        end
    end
    
    % summery of internal loads calculation 
    
    Internal_Stresses.VonMise_skin=Von_skn;
    
    Internal_Stresses.VonMise_spar=Von_spr;
    
    Internal_Stresses.Skin_BucklingCritical=sigmab_skn;
    
    Internal_Stresses.Skin_buckling_shear=taub_skn;
    
    Internal_Stresses.Skin_buckling_principal=sigma_pp;
    
    Internal_Stresses.Stringer_compressive=sigma_strg;
    
    Internal_Stresses.Stringer_crippling=sigma_crip;
    
    Internal_Stresses.Stringer_columnbuckling=sigma_col;
    



end