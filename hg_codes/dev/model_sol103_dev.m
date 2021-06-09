%% sizing parameters thickness1 = spar, thickness2 = skin

    thickness1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
        x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
        x(21),x(22)];

    thickness2=[x(23),x(24),x(25),x(26),x(27),x(28),x(29),x(30),x(31),x(32)...
        x(33),x(34),x(35),x(36),x(37),x(38),x(39),x(40),x(41),x(42)...
        x(43),x(44)];

    Astrg=[x(45),x(46),x(47),x(48),x(49),x(50),x(51),x(52),x(53),x(54)...
        x(55),x(56),x(57),x(58),x(59),x(60),x(61),x(62),x(63),x(64)...
        x(65),x(66)];
    
   
  %% Wingbox 1 - right and control surf.

    Wingbox_right = awi.model.LiftingSurface;
    Wingbox_right.Name = 'A320Wing_right';
    Wingbox_right.Origin=[15,0,0];
    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox_right.ActiveSet = 'sSet';

    %Wing dimensions
    Wingbox_right.SpanVector  = 'Y';
    Wingbox_right.Span        = 15;   %34.1/2;
    Wingbox_right.LESweep     = [27, 27];
    Wingbox_right.LESweep_eta = [0, 1];
    Wingbox_right.TESweep     = [0, 16.59, 16.59];
    Wingbox_right.TESweep_eta = [0, 0.27, 1];
    Wingbox_right.RootChord   = 6;
    
    % testing non-sweep    
%     Wingbox_right.SpanVector  = 'Y';
%     Wingbox_right.Span        = 15;   %34.1/2;
%     Wingbox_right.LESweep     = [0, 0];
%     Wingbox_right.LESweep_eta = [0, 1];
%     Wingbox_right.TESweep     = [0, 0, 0];
%     Wingbox_right.TESweep_eta = [0, 0.27, 1];
%     Wingbox_right.RootChord   = 4;
    
    %Dihedral 
    Wingbox_right.Dihedral=[5,5];
    Wingbox_right.Dihedral_eta=[0,1];
    

    %Make sure the beam is at the midchord
    all_eta           = Wingbox_right.Eta_;
    Wingbox_right.BeamLoc     = repmat(0.4, size(all_eta));
%     Wingbox_right.BeamLoc     = [0.34,0.4,0.4];
    Wingbox_right.BeamLoc_eta = all_eta;

    %Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
    FrontSpar_right = awi.model.Spar;
    FrontSpar_right.XLoc = [0.15, 0.15];
    FrontSpar_right.Eta  = [0   , 1];
    RearSpar_right = awi.model.Spar;
    RearSpar_right.XLoc = [0.65, 0.65];
    RearSpar_right.Eta  = [0   , 1];

    Wingbox_right.add([FrontSpar_right, RearSpar_right]);

    %Define internal layout
    Wingbox_right.RibPitch      = 0.65;
    Wingbox_right.StringerPitch = 0.15;

    %Make the material
    E  = 70e9; %[N/m^2], typical YM of aluminium
    nu = 0.333;
    rho=2810;
    Mat = awi.model.Material;
    Mat.E  = E;
    Mat.Nu = nu;
    Mat.G  = E / (2 * (1 + nu));
    Mat.Rho=rho;
    Wingbox_right.Material_eta = [0, 1];
    Wingbox_right.Material     = [Mat, Mat];
    
    build(Wingbox_right)
    
    
    %% Create discretised boxbeam with varied cross section prperties along the span 

    NumSec=22;
    d_strg=sqrt(Astrg/0.36);
    t_strg=0.12*d_strg;

    % etaS=linspace(0,Wingbox_right.Span,NumSec);

    % set width and height array 
    YData=Wingbox_right.YData;
    SparWidth=Wingbox_right.Chord*0.5;

    RootH=Wingbox_right.Chord(1)*0.15;
    MidH=Wingbox_right.Chord(2)*0.12;
    TipH=Wingbox_right.Chord(end)*0.11;


    % set up eta values
    eta1_=linspace(0,0.27,7);
    eta2_=linspace(0.27,1,16);
    etaS=[eta1_(1:end-1),eta2_(1:end)];

    RData=Wingbox_right.RData;
    eta_R=RData/RData(end);
    eta_Y=YData/YData(end);
    etaRS=interp1(eta_Y,eta_R,etaS);

    Bwidth=interp1(RData/RData(end),SparWidth,etaRS);
%     Bheight=interp1(RData/RData(end),0.79*[RootH,MidH,TipH],etaRS);
    
    Bheight=interp1(RData/RData(end),0.79*[RootH,RootH,RootH],etaRS);

    % stringer pitch 
    strg_n=0.24;

    %intialise data array
    A_val=zeros(1,NumSec);
    Ixx_val=zeros(1,NumSec);
    Izz_val=zeros(1,NumSec);
    J_val=zeros(1,NumSec);
    
    %NSM
    NSM_val=zeros(1,NumSec);
    NSI_val=zeros(1,NumSec);
    

    for ii=1:NumSec

        boxname=strcat('Box',string(ii));
        boxname=awi.model.BoxBeam;
        boxname.BoxType='SymmetricBox';
        boxname.Height=Bheight(ii);
        boxname.Width=Bwidth(ii);
        boxname.CoverThickness=thickness2(ii);
        boxname.SparThickness=thickness1(ii);

        NumStrg=floor(Bwidth(ii)/strg_n);

        ts=t_strg(ii);
        ds=d_strg(ii);
        hs=Bheight(ii);
        ws=Bwidth(ii);

        Istrg_xx=((ts*ds^3/12)+(ds*ts^3/12 + ds*ts*(ds/2)^2)*2 + 3*ds*ts*(hs/2-ds/2)^2)*NumStrg*2;
        Istrg_zz_=(ds*ts^3/12 + (ts*ds^3/12 + ts*ds*(ds/2)^2)*2);

        if mod(NumStrg,2)==0
            offset=0.12:strg_n:ws/2;
            Istrg_zz=(NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;

        elseif mod(NumStrg,2)==1
            offset=0:strg_n:ws/2;
            Istrg_zz=(NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;

        end

        getGeometricProps(boxname)
        A_val(ii)=boxname.Abb+0;
        Ixx_val(ii)=boxname.Ixx+Istrg_xx;
        Izz_val(ii)=boxname.Izz+Istrg_zz;
        J_val(ii)=boxname.Jbb;
        
        % NSM       
        NSM_val(ii)=boxname.NSM;
        NSI_val(ii)=boxname.NSI;

    end


    eta_=etaRS;
    Wingbox_right.A   =  A_val;
    Wingbox_right.A_eta=eta_;

    Wingbox_right.I11 = Izz_val;
    Wingbox_right.I11_eta=eta_;

    Wingbox_right.I22 = Ixx_val;
    Wingbox_right.I22_eta = eta_;

    Wingbox_right.J   = J_val;
    Wingbox_right.J_eta= eta_;
    
    % NSM and NSI
    Wingbox_right.NSM=NSM_val;
    Wingbox_right.NSM_eta= eta_;
    
    Wingbox_right.NSI=NSI_val;
    Wingbox_right.NSI_eta= eta_;

  
    % Aeropanel definition
    
    % AeroPanelLength
    %     NumAeroPanel
    Wingbox_right.NumAeroPanel=20;
   

    %% Mass definition
    
    % total wing mass
    
    [wing_mass,total_mass]=Mass_calc_v1(x);
    
    wingmass_eta=0.04:0.04:1;
    Mwidth=interp1(RData/RData(end),SparWidth,wingmass_eta);
    
    mass_set=(total_mass-wing_mass)*(Mwidth)/sum(Mwidth);
      
%     m=total_mass/19;
    
    for i=1:1:25
        handle=strcat('PM_right','i');
        handle=awi.model.PointMass;
        handle.SOffset=0+i*0.04;
        handle.Mass=mass_set(i);
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
        handle.MassGroup='Group1';
        Wingbox_right.add(handle);

    end
    
%     % Landing Gear
%     Landing_Gear=awi.model.PointMass;   
%     Landing_Gear.SOffset=0.23;
%     Landing_Gear.Mass=2491/2;
%     Wingbox_right.add(Landing_Gear);
%     
%     % engine
%     engine=awi.model.PointMass;   
%     engine.SOffset=0.25;
% %     engine.XOffset=-1.5;
% %     engine.YOffset=4;
% %     engine.ZOffset=-1.8;
%     engine.Mass=3681+1239/2;
%     Wingbox_right.add(engine);
    
 
    %% attachments - engine
    
    Engine=awi.model.BluffBody;
    Engine.Name='Engine';
    
    % cylinder body
    Engine.Radius=[1.4, 1.4, 1];
    Engine.Eta =  [0, 0.6, 1];
    
%     Engine.Origin = [16.238147-3.5, 4.05 , -1.8];
    Engine.Origin = [16.238147-3.5, 4.05 , 0.35432909];
  
    Engine.Length = 3.5;
%     Engine.XOffset=16.471008-3.5;
%     Engine.YOffset=5.8170588;
%     Engine.ZOffset=-2;
    
    %Make the material
    E1  = 76e9; %[N/m^2],set as a rigid body
    nu = 0.333;
    Mat1 = awi.model.Material;
    Mat1.E  = E1;
    Mat1.Nu = nu;
    Mat1.G  = E1 / (2 * (1 + nu));
    Mat1.Rho=2800;
    
    
    % use the strong material
    Engine.Material_eta = [0, 1];
    Engine.Material     = [Mat1, Mat1];

    Engine.A   = 0.04432;
    Engine.I11 = 0.002;
    Engine.I22 = 0.002;
    Engine.J   = 0.001636;
    
    %Aeropanel althoufh it is useless now
    Engine.AeroPanelLength=0.5;
    
    % add engine mass
    engine_mass=awi.model.PointMass;   
    engine_mass.SOffset=0.4;
    engine_mass.Mass=3681;
    Engine.add(engine_mass);
    
    % add pylon
    pylon_mass=awi.model.PointMass;   
    pylon_mass.SOffset=0.9;
    pylon_mass.Mass=1239/2;
    Engine.add(pylon_mass);

    build(Engine)
      
    Wingbox_right.add(Engine)
    build(Wingbox_right);
%     
%     draw(Wingbox_right)

    %% Build aircraft model
    Aircraft = awi.model.Aircraft;

%     Aircraft.add(Body);
        
    Aircraft.add(Wingbox_right)

    %The analysis methods require an 'awi.model.Aircraft' object
    % This is because some information is only known at the aircraft level,
    % e.g. all-up mass, reference span, area, etc.
    % Aircraft = awi.model.Aircraft;
    % Aircraft.add(LS);

    Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea]);
    Aircraft.RefSpan  = Wingbox_right.Span;
    Aircraft.RefChord = Wingbox_right.RootChord;

    %% Generate the FEM - SOL 103 Modal analysis

    run_folder = [
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\half_ac_dev']; %[-], folder for exporting the NASTRAN model

    % Convert to a finite element model and draw it
    FEM_wing = convertToFE(Aircraft);

    % %Export it to a file
    export(FEM_wing, run_folder);
    
    
    
    %% insert SPC
    
    NastranMethods1 = awi.methods.Nastran;
    NastranMethods1.AnalysisModel = FEM_half;
    MassCases=awi.model.MassCases.empty;
    RefGrid=NastranMethods1.RefNode;

    fid = fopen(strcat(run_folder,'\NastranHeaderFile.dat'),'r');
    
    cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true );
    cac=cac{1};
    fclose( fid );
    
    mid_line0=22;
    mid_line=length(cac)-1; 
    fid = fopen( strcat(run_folder,'\half_model_103.dat'), 'w' );
    
    % write existing head file
    for jj = 1 : mid_line0
        fprintf(fid,'%s\n',cac{jj})       
    end
    
%     fprintf(fid,'SUPORT = %i\r\n',201)
    fprintf(fid,'SPC = %i\r\n',202)
    
    for ii =  mid_line0 : mid_line 
        fprintf(fid,'%s\n',cac{ii})       
    end
    
    
    %write spc part
    line='$.1.....2.......3.......4.......5.......6.......7.......8.......9.......10......\r\n';
    fprintf(fid,line);
    
    dof = 35;
    spc = 246;
    SPC_id = ID0;
    TrimData.RefGrid = RefGrid;
    TrimData.SPC_id  = SPC_id;
    TrimData.SPCdof  = spc;
    TrimData.SUPdof  = dof;
                        
%     spc_format='%-8s%-8i%-8i%-8i\r\n';
%     fprintf(fid,spc_format,'SPC1',1,123456,1005);
    
    fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPCADD', TrimData.SPC_id+2, TrimData.SPC_id, TrimData.SPC_id+1);
    %   - SPC
    rows=floor((numel(TrimData.RefGrid)-3)/5);
    remd=rem((numel(TrimData.RefGrid)-3),5);
    fprintf(fid, '%-8s%-8i%-8i%-8i%-8i%-8i\r\n', 'SPC1', TrimData.SPC_id, TrimData.SPCdof, TrimData.RefGrid(1:3).GID);
    fprintf(fid, '        %-8i%-8i%-8i%-8i%-8i\r\n', TrimData.RefGrid(4:4+rows*5-1).GID);
    format_last=strcat(repmat('%-8i',1,remd+1),'\r\n');
    fprintf(fid, format_last, '', TrimData.RefGrid(4+rows*5:end).GID);
    
    %charles added SPC 2
    % -SPC 2 at COG
    fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', TrimData.SPC_id+1, 1246, TrimData.RefGrid(34).GID);
    fprintf(fid, '%-8s%-8i%-8i\r\n', 'SUPORT', TrimData.RefGrid(34).ID, TrimData.SUPdof);
    
    % write the end 
    fprintf(fid,'%s\n',cac{end});
    
%     NastranMethods1.runNastran(trimFile);
    
%      NastranMethods1 = awi.methods.Nastran;
% 
%     
% %     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\A320_half_model_SOL144*.*')
%    
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.xdb');
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.h5');
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.log');
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.f06');
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.f04');
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.h5.*');
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.log.*');
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.f06.*');
%     delete('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test0\A320_half_model_SOL144*.f04.*');
%       
%     NastranMethods1.runNastran('');


    %% Run the analysis - SOL 145 flutter 
    
%     FlightPoint=awi.model.FlightPoint;
%     
%     FlightPoint.Mach=0.78;
% %     FlightPoint.AcvELOCITY=50;
%     FlightPoint.Altitude = 10000;
%     
%     getFlightPointData(FlightPoint)
%     
%     flutterFile = NastranMethods1.writeFlutterFile(Aircraft, FlightPoint, MassCases, '123456', run_folder, 'RequestModeshapes',true,'FlutterMethod','pk');
%     
%     NastranMethods1.runNastran(flutterFile);
% 
% 
% 
% %% flutter results Vg Vf;
% 
% flutter_data = h5read('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\normal_mode\flutter_analysis.h5','/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');
% 
% modes=[1:35];
% 
% for i=1:35
%     Modes_pt=modes(i);
%     [index1,~]=find(flutter_data.POINT==Modes_pt);
%     
%     velocity=flutter_data.VELOCITY(index1);
%     frequency=flutter_data.FREQUENCY(index1);
%     
%     
%     plot(velocity,frequency,'s-')
%        
%     set(gcf,'Color','w')
%     xlabel('Velocity','Interpreter','latex','FontSize',12)
%     ylabel('Frequency','Interpreter','latex','FontSize',12)
%     hold on
%     
% end
%     
% for i=1:35
%     Modes_pt=modes(i);
%     [index1,~]=find(flutter_data.POINT==Modes_pt);
%     
%     velocity=flutter_data.VELOCITY(index1);
%     
%     damping=flutter_data.DAMPING(index1);
%     
%     plot(velocity,damping,'s-')
%     set(gcf,'Color','w')
%     xlabel('Velocity','Interpreter','latex','FontSize',12)
%     ylabel('Damping','Interpreter','latex','FontSize',12)
%     hold on
%     
% end    
%     
    
    %% Example
    
%        fid = fopen( old_filespec );
%     cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true );
%     fclose( fid )
%     fid = fopen( new_filespec, 'w' );
%     for jj = 1 : insert_here
%         fprintf( fid, '%s\n', cac{jj} );
%     end
%     fprintf( fid, '%s\n', new_line );
%     for jj = insert_here+1 : length(cac)
%         fprintf( fid, '%s\n', cac{jj} );
%     end
%     fclose( fid );
%     
%     
    
    
    
    