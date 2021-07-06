function [Res_TrimLoadcase1,Res_TrimLoadcase2,Res_TrimLoadcase3,Res_TrimLoadcase4,Res_TrimLoadcase5,Res_TrimLoadcase6]=Static_Load_Envelop(run_folder,Aircraft,FEM_full)

 %% LoadCase 1: cruising altitude 2.5g pull up
    TrimLoadcase1 = awi.model.LoadCase;
    
    acMass = 94000;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;
    
%     FlightPoint=awi.model.FlightPoint;
%     FlightPoint.Mach=0.78;
%     FlightPoint.AcVelocity=FlightPoint.Mach*340;
%     FlightPoint.Altitude = 36000;
%     getFlightPointData(FlightPoint,'ISA');

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
    
    
    %% steady state for gust analysis, 36000 ft, g
    TrimLoadcase4 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 36000;
    mach_number       = 0.78;
    aircraft_velocity = 0.78*340;
    flap_angle=0;

    TrimLoadcase4.Name = 'A321_36000ft_g';

    TrimLoadcase4.Altitude   = altitude;
    TrimLoadcase4.Mach       = mach_number;
    TrimLoadcase4.AcVelocity = aircraft_velocity;
    TrimLoadcase4.AcMass = acMass;

    TrimLoadcase4.PitchAngle=0;
    TrimLoadcase4.RollAngle =0;
    TrimLoadcase4.ID = 1030;
    TrimLoadcase4.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase4.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase4)
    
    %% steady state for gust analysis, 20000 ft, g
    TrimLoadcase5 = awi.model.LoadCase;
    
    acMass = 500;
    altitude          = 20000;
    mach_number       = 0.6;
    aircraft_velocity = 0.6*340;
    flap_angle=0;
    
    TrimLoadcase5.Name = 'A321_20000ft_g';
    
    TrimLoadcase5.Altitude   = altitude;
    TrimLoadcase5.Mach       = mach_number;
    TrimLoadcase5.AcVelocity = aircraft_velocity;
    TrimLoadcase5.AcMass = acMass;
    
    TrimLoadcase5.PitchAngle=0;
    TrimLoadcase5.RollAngle =0;
    TrimLoadcase5.ID = 1030;
    TrimLoadcase5.LoadFactor = 1;
    
    % CS deflection - flap deflection angle in radian
    TrimLoadcase5.CsDeflection=flap_angle*pi/180;
    
    build(TrimLoadcase5)
    
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
    
    %% Run analysis
    
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_full;
    MassCases=awi.model.MassCases.empty;
    ID0=200;
    
    % sizing case 1
    trimFile1 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase1, MassCases,run_folder,'DatFilename','A321_cruise_2p5g');
    
    % sizing case 2
    trimFile2 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase2, MassCases,run_folder,'DatFilename','A321_sealevel_2p5g');
    
    % sizing case 3
    trimFile3 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase3, MassCases,run_folder,'DatFilename','A321_sealevel_neg_1g');
    
    % sizing case 4
    trimFile4 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase4, MassCases,run_folder,'DatFilename','A321_36000ft_1g');
    
    % sizing case 5
    trimFile5 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase5, MassCases,run_folder,'DatFilename','A321_20000ft_1g');
    
    % sizing case 6
    trimFile6 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase6, MassCases,run_folder,'DatFilename','A321_3000ft_1g');
    
    
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
    NastranMethods1.runNastran(trimFile4);
    NastranMethods1.runNastran(trimFile5);
    NastranMethods1.runNastran(trimFile6);
    
    %% Result extraction 
    
    Trim1_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_cruise_2p5g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim2_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_sealevel_2p5g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim3_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_sealevel_neg_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim4_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_36000ft_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim5_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_20000ft_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    Trim6_res=NastranMethods1.extractNastranResults(strcat(run_folder,'\A321_3000ft_1g.f06'),'ReadF06',true,'ReadHDF5',false);
    
    WingNodes=NastranMethods1.WingNode;
    
    index1=ismember([Trim1_res.f06data.Bendingmoment.UGRID],[WingNodes(1:end).GID]);
    index2=ismember([Trim1_res.f06data.Bendingmoment.LGRID],[WingNodes(end).GID]);
    
    % extract forces:
    
    % forces from trim 1: 2.5 g 36000ft
    Res_TrimLoadcase1.M_P1=[Trim1_res.f06data.Bendingmoment.UMPLN1(index1),Trim1_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Res_TrimLoadcase1.M_P2=[Trim1_res.f06data.Bendingmoment.UMPLN2(index1),Trim1_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Res_TrimLoadcase1.T=[Trim1_res.f06data.Bendingmoment.UTORQUE1(index1),Trim1_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque

    Res_TrimLoadcase1.S_P1=[Trim1_res.f06data.Bendingmoment.USPLN1(index1),Trim1_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Res_TrimLoadcase1.S_P2=[Trim1_res.f06data.Bendingmoment.USPLN2(index1),Trim1_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    % forces from trim 2: 2.5 g Sea level
    Res_TrimLoadcase2.M_P1=[Trim2_res.f06data.Bendingmoment.UMPLN1(index1),Trim2_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Res_TrimLoadcase2.M_P2=[Trim2_res.f06data.Bendingmoment.UMPLN2(index1),Trim2_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Res_TrimLoadcase2.T=[Trim2_res.f06data.Bendingmoment.UTORQUE1(index1),Trim2_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque

    Res_TrimLoadcase2.S_P1=[Trim2_res.f06data.Bendingmoment.USPLN1(index1),Trim2_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Res_TrimLoadcase2.S_P2=[Trim2_res.f06data.Bendingmoment.USPLN2(index1),Trim2_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    % forces from trim 3: -g sea level
    Res_TrimLoadcase3.M_P1=[Trim3_res.f06data.Bendingmoment.UMPLN1(index1),Trim3_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Res_TrimLoadcase3.M_P2=[Trim3_res.f06data.Bendingmoment.UMPLN2(index1),Trim3_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Res_TrimLoadcase3.T=[Trim3_res.f06data.Bendingmoment.UTORQUE1(index1),Trim3_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
    
    Res_TrimLoadcase3.S_P1=[Trim3_res.f06data.Bendingmoment.USPLN1(index1),Trim3_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Res_TrimLoadcase3.S_P2=[Trim3_res.f06data.Bendingmoment.USPLN2(index1),Trim3_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    % forces from trim 4: g 36000ft
    Res_TrimLoadcase4.M_P1=[Trim4_res.f06data.Bendingmoment.UMPLN1(index1),Trim4_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Res_TrimLoadcase4.M_P2=[Trim4_res.f06data.Bendingmoment.UMPLN2(index1),Trim4_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Res_TrimLoadcase4.T=[Trim4_res.f06data.Bendingmoment.UTORQUE1(index1),Trim4_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
    
    Res_TrimLoadcase4.S_P1=[Trim4_res.f06data.Bendingmoment.USPLN1(index1),Trim4_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Res_TrimLoadcase4.S_P2=[Trim4_res.f06data.Bendingmoment.USPLN2(index1),Trim4_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    % forces from trim 5: g 20000ft
    Res_TrimLoadcase5.M_P1=[Trim5_res.f06data.Bendingmoment.UMPLN1(index1),Trim5_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Res_TrimLoadcase5.M_P2=[Trim5_res.f06data.Bendingmoment.UMPLN2(index1),Trim5_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Res_TrimLoadcase5.T=[Trim5_res.f06data.Bendingmoment.UTORQUE1(index1),Trim5_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
    
    Res_TrimLoadcase5.S_P1=[Trim5_res.f06data.Bendingmoment.USPLN1(index1),Trim5_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Res_TrimLoadcase5.S_P2=[Trim5_res.f06data.Bendingmoment.USPLN2(index1),Trim5_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    % forces from trim 6: g 3000ft
    Res_TrimLoadcase6.M_P1=[Trim6_res.f06data.Bendingmoment.UMPLN1(index1),Trim6_res.f06data.Bendingmoment.LMPLN1(index2)]; % in-plane moment
    Res_TrimLoadcase6.M_P2=[Trim6_res.f06data.Bendingmoment.UMPLN2(index1),Trim6_res.f06data.Bendingmoment.LMPLN2(index2)]; % out of plane moment
    Res_TrimLoadcase6.T=[Trim6_res.f06data.Bendingmoment.UTORQUE1(index1),Trim6_res.f06data.Bendingmoment.LTORQUE1(index2)];% torque
    
    Res_TrimLoadcase6.S_P1=[Trim6_res.f06data.Bendingmoment.USPLN1(index1),Trim6_res.f06data.Bendingmoment.LSPLN1(index2)]; % in plane shear
    Res_TrimLoadcase6.S_P2=[Trim6_res.f06data.Bendingmoment.USPLN2(index1),Trim6_res.f06data.Bendingmoment.LSPLN2(index2)]; % out of plane shear
    
    
    
    
    
    





end