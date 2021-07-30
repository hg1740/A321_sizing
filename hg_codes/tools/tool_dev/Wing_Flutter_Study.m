% Wing model


Wing = awi.model.LiftingSurface;
Wing.Name = 'A320Wing_right';
Wing.Origin=[0,0,0]; %15

%Use the Leading/Trailing edge sweep to define the planform
Wing.ActiveSet = 'sSet';

%Tail wing dimensions
Wing.SpanVector  = 'Y';
Wing.Span        = 9;
Wing.LESweep     = [0, 0];
Wing.LESweep_eta = [0, 1];
Wing.TESweep     = [0,  0];
Wing.TESweep_eta = [0,  1];
Wing.RootChord   = 1;

%Make sure the beam is at the midchord
all_eta           = Wing.Eta_;
Wing.BeamLoc     = repmat(0.4, size(all_eta));
Wing.BeamLoc_eta = all_eta;


% Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
FrontSpar_root_right = awi.model.Spar;
FrontSpar_root_right.XLoc = [0.15, 0.15];
FrontSpar_root_right.Eta  = [0   , 1];
RearSpar_root_right = awi.model.Spar;
RearSpar_root_right.XLoc = [0.65, 0.65];
RearSpar_root_right.Eta  = [0   , 1];
Wing.add([FrontSpar_root_right, RearSpar_root_right]);

%Define internal layout
Wing.RibPitch      = 0.65;
Wing.StringerPitch = 0.15;

%Make the connector material
E0  = 70e9; %[N/m^2], typical YM of IM7 composite
nu0 = 0.333;
rho0=2800;
Mat_conn = awi.model.Material;
Mat_conn.E  = E0;
Mat_conn.Nu = nu0;
Mat_conn.G  = E0 / (2 * (1 + nu0));
Mat_conn.Rho=rho0;

% material properties
Wing.Material_eta = [0, 1];
Wing.Material     = [Mat_conn, Mat_conn];

% Define box beam corss section
Wing_box=awi.model.BoxBeam;
Wing_box.BoxType='SymmetricBox';
Wing_box.Height=0.15;
Wing_box.Width=0.5;
Wing_box.CoverThickness=0.002;
Wing_box.SparThickness=0.020;
getGeometricProps(Wing_box)

Wing.BoxBeam = Wing_box;
Wing.A   = Wing_box.Abb;
Wing.I11 = Wing_box.Izz;
Wing.I22 = Wing_box.Ixx;
Wing.J   = Wing_box.Jbb;

% 
% for i=1:1:5
%     Inboard_Mass_Nm=strcat('PM_tail_R','i');
%     Inboard_Mass=awi.model.PointMass;
%     Inboard_Mass.SOffset=i*0.2;
%     Inboard_Mass.Mass=10;
%     Inboard_Mass.MassGroup='Group3';
%     Wing.add(Inboard_Mass);
% end

% Aeropanel definition
Wing.AeroPanelLength=[];
Wing.NumAeroPanel=10;
Wing.AeroPanelAR=1;

% Jig shape
Wing.Twist=[0,0];

Wing.Twist_eta=[0,1];

build(Wing);


FWT = insertWingFold(Wing, 'FlareAngle', 12.5, 'FoldAngle',0,'EtaFold',0.75);
% FWT.HingeStiffness = [1e14 1e14 1e14 1e14 100 1e14];
FWT.HingeStiffness = [0 0 0 0 1 0];

FWT.AeroPanelLength=[];
FWT.NumAeroPanel=10;
FWT.AeroPanelAR=1;


%Make the material for FWT

E_fwt  = 70e9; %[N/m^2], typical YM of aluminium
nu_fwt = 0.333;
rho_fwt=2800;
Mat_fwt = awi.model.Material;
Mat_fwt.E  = E_fwt;
Mat_fwt.Nu = nu_fwt;
Mat_fwt.G  = E_fwt / (2 * (1 + nu_fwt));
Mat_fwt.Rho=rho_fwt;

FWT.Material_eta = [0, 1];
FWT.Material     = [Mat_fwt, Mat_fwt];


FWT_box=awi.model.BoxBeam;
FWT_box.BoxType='SymmetricBox';
FWT_box.Height=0.15;
FWT_box.Width=0.5;
FWT_box.CoverThickness=0.002;
FWT_box.SparThickness=0.002;
getGeometricProps(FWT_box)

FWT_eta_=[0,1];

FWT.A   =  FWT_box.Abb;
FWT.A_eta=FWT_eta_;

FWT.I11 = FWT_box.Izz;
FWT.I11_eta=FWT_eta_;

FWT.I22 = FWT_box.Ixx;
FWT.I22_eta = FWT_eta_;

FWT.J   = FWT_box.Jbb;
FWT.J_eta= FWT_eta_;


% % tip mass comment: it makes flutter worse!!!


% FWT_Mass=awi.model.PointMass;
% FWT_Mass.SOffset=0.5;
% FWT_Mass.Mass=100;
% FWT_Mass.MassGroup='Group3';
% FWT.add(FWT_Mass);


Aircraft = awi.model.Aircraft;

Aircraft.add(Wing);

Aircraft.RefArea = Wing.SurfaceArea + FWT.SurfaceArea;
Aircraft.RefSpan  = Wing.Span + FWT.Span;
Aircraft.RefChord = Wing.RootChord; %mean aerodynamic chord = 0.697 for A321 wing;
       
FEM_full=convertToFE(Aircraft);

% run folder
run_folder='C:\Git\A321_sizing\hg_codes\results\flutter';

% %Export it to a file
export(FEM_full, run_folder);

draw(FEM_full)


%% Flutter loadcase

NastranMethods1 = awi.methods.Nastran;
NastranMethods1.AnalysisModel = FEM_full;
MassCases=awi.model.MassCases.empty;

FlightPoint=awi.model.FlightPoint;

FlightPoint.Mach=0.1;
FlightPoint.AcVelocity=50;
FlightPoint.Altitude = 20000;

getFlightPointData(FlightPoint)

flutterFile = NastranMethods1.writeFlutterFile(Aircraft, FlightPoint, MassCases, '123456', run_folder, 'RequestModeshapes',true,'FlutterMethod','pk');

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

NastranMethods1.runNastran(flutterFile);




%% flutter results Vg Vf;

% run_folder='D:\pdf_nastran\pdf_nastran\FFAST_Test_FWT\FFAST Test\Source\Wing_only';

flutter_data = h5read(strcat(run_folder,'\flutter_analysis.h5'),'/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');

modes=[1:35];

figure

for i=1:7
    Modes_pt=modes(i);
    [index1,~]=find(flutter_data.POINT==Modes_pt);
    
    velocity=flutter_data.VELOCITY(index1);
    frequency=flutter_data.FREQUENCY(index1);
    
    
    plot(velocity,frequency,'s-','LineWidth',1)
       
    set(gcf,'Color','w')
    xlabel('Velocity','Interpreter','latex','FontSize',12)
    ylabel('Frequency','Interpreter','latex','FontSize',12)
    hold on
%     axis([260 400 0 25])
%     
end

legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','Interpreter','latex','FontSize',10)
 
 
figure
for i=1:7
    Modes_pt=modes(i);
    [index1,~]=find(flutter_data.POINT==Modes_pt);
    
    velocity=flutter_data.VELOCITY(index1);
    
    damping=flutter_data.DAMPING(index1);
    
    plot(velocity,damping,'s-','LineWidth',1)
    set(gcf,'Color','w')
    xlabel('Velocity','Interpreter','latex','FontSize',12)
    ylabel('Damping','Interpreter','latex','FontSize',12)
    hold on
    
end    

hold on 
plot([0,300],[0,0],'k--')
legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','Interpreter','latex','FontSize',10)



    
    %% Run the analysis- SOL 103 
% %       insert SPC
%     fid = fopen(strcat(run_folder,'\NastranHeaderFile.dat'));
%     
%     cac = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true );
%     cac=cac{1};
%     fclose( fid )
%     
%     mid_line0=22;
%     mid_line1=33;
%     
%     end_line=length(cac)-1; 
%     fid = fopen( strcat(run_folder,'\Modal_analysis.dat'), 'w' );
%     
%     % write existing head file
%     for jj = 1 : mid_line0
%         fprintf(fid,'%s\n',cac{jj})       
%     end
%     
%     fprintf(fid,'SPC = 1 \r\n')
%     
%     for ii =  mid_line0 : mid_line1
%         fprintf(fid,'%s\n',cac{ii})       
%     end
%     
%     % correct eigen value extraction methods
%     fprintf(fid, '%-8s%-8i%-8s%-#8.3g%-16s%-8i\r\n', 'EIGR', 117, 'MGIV', 0, blanks(16),30);
%     
%     
%     for ii =  mid_line1+2 : end_line
%         fprintf(fid,'%s\n',cac{ii})
%     end
%     
% 
%     %write spc part
%     line='$.1.....2.......3.......4.......5.......6.......7.......8.......9.......10......\r\n';
%     fprintf(fid,line);
%     
%     spc_format='%-8s%-8i%-8i%-8i\r\n';
%     fprintf(fid,spc_format,'SPC1',1,123456,1005);
%     
%     % write the end 
%     fprintf(fid,'%s\n',cac{end})
%     
%     Sol103_file='C:\Git\A321_sizing\hg_codes\results\flutter\Modal_analysis.dat';
% 
%     delete(strcat(run_folder, '\Modal_analysis.xdb'));
%     delete(strcat(run_folder, '\Modal_analysis.h5'));
%     delete(strcat(run_folder, '\Modal_analysis.log'));
%     delete(strcat(run_folder, '\Modal_analysis.f06'));
%     delete(strcat(run_folder, '\Modal_analysis.f04'));
%     delete(strcat(run_folder, '\Modal_analysis.op2'));
%     
%     delete(strcat(run_folder, '\Modal_analysis.xdb.*'));
%     delete(strcat(run_folder, '\Modal_analysis.h5.*'));
%     delete(strcat(run_folder, '\Modal_analysis.log.*'));
%     delete(strcat(run_folder, '\Modal_analysis.f06.*'));
%     delete(strcat(run_folder, '\Modal_analysis.f04.*'));
%     delete(strcat(run_folder, '\Modal_analysis.op2.*'));
% 
%     
%     NastranMethods1.runNastran(Sol103_file);
%     
%     
%  %% Result plot   
% % load the model
% model = mni.import_matran(fullfile(run_folder,'modal_analysis.dat'));
% model.draw
% 
% % get modal data
% res_modeshape = mni.result.f06(fullfile(run_folder,'modal_analysis.f06')).read_modeshapes;
% res_freq = mni.result.f06(fullfile(run_folder,'modal_analysis.f06')).read_modes;
% %% apply deformation result
% %pick which mode to plot (1->6)
% modeshape_num = 3;
% 
% % apply the modal deformation
% [~,i] = ismember(model.GRID.GID,res_modeshape.GID(modeshape_num,:));
% model.GRID.Deformation = [res_modeshape.T1(modeshape_num,i);...
%     res_modeshape.T2(modeshape_num,i);res_modeshape.T3(modeshape_num,i)];
% 
% % animate the plot to show the mode shape. 'scale' is the scale of the
% % displacement. 'Frequency' is the frequency.
% %
% % PRESS SPACE TO END THE ANIMATION
% %
% model.update()
% model.animate('Frequency',2,'Scale',1) 












% 
% %% result 101
% 
% % load the model
% model = mni.import_matran(fullfile(run_folder,'SOL101.dat'));
% model.draw
% % extract the displacment data
% res_disp = mni.result.f06(fullfile(run_folder,'SOL101.f06')).read_disp;
% 
% % apply deformation results
% [~,i] = ismember(model.GRID.GID,res_disp.GP);
% model.GRID.Deformation = [res_disp.dX(:,i);res_disp.dY(:,i);res_disp.dZ(:,i)];
% 
% % update the figure (Scale is the normalised scale of displacements 
% % to plot)
% model.update('Scale',0.000001)




