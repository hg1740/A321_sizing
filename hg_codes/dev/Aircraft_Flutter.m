%% Import sized aircraft model

Aircraft_Model=load('C:\Git\A321_sizing\hg_codes\Data\Sizing_with_upated_hinge_lock\A126\AR10\Res_AR10_Eta_100_Model.mat');

Param=Aircraft_Model.Param;


if isfield(Param,'FWT')
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT]=Aircraft_Models_v1(Param);
    
    
else
    
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v1(Param);
    
end


draw(FEM_full);

%% run folder

run_folder='C:\Git\A321_sizing\hg_codes\results\Aircraft_Flutter';

export(FEM_full,run_folder)

%% Flutter loadcase

NastranMethods1 = awi.methods.Nastran;
NastranMethods1.AnalysisModel = FEM_full;
MassCases=awi.model.MassCases.empty;

FlightPoint=awi.model.FlightPoint;

FlightPoint.Mach=0.78;
FlightPoint.AcVelocity=50;
FlightPoint.Altitude = 36000;

getFlightPointData(FlightPoint)

flutterFile = NastranMethods1.writeFlutterFile(Aircraft, FlightPoint, MassCases, '1246', run_folder, 'RequestModeshapes',true,'FlutterMethod','pk');

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

flutter_data = h5read(strcat(run_folder,'\flutter_analysis.h5'),'/NASTRAN/RESULT/AERODYNAMIC/FLUTTER/SUMMARY');

modes=[1:35];

figure

for i=1:20
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

legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','Interpreter','latex','FontSize',10)

figure
for i=1:20
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
legend('1st Mode','2nd Mode','3rd Mode','4th Mode','5th Mode','6th Mode','7th Mode','8th Mode','9th Mode','10th Mode','11th Mode','12th Mode','Interpreter','latex','FontSize',10)%     %% Run the analysis- SOL 103
%     % insert SPC
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
% %     fprintf(fid,'SPC = 1 \r\n')
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
% %     spc_format='%-8s%-8i%-8i%-8i\r\n';
% %     fprintf(fid,spc_format,'SPC1',1,123456,1008);
%     
%     % write the end
%     fprintf(fid,'%s\n',cac{end})
%     
%     Sol103_file=strcat(run_folder,'\Modal_analysis.dat');
%     
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
%     NastranMethods1.runNastran(Sol103_file);
%     
%     
%     %% Result plot
%     % load the model
%     model = mni.import_matran(fullfile(run_folder,'modal_analysis.dat'));
%     model.draw
%     
%     % get modal data
%     res_modeshape = mni.result.f06(fullfile(run_folder,'modal_analysis.f06')).read_modeshapes;
%     res_freq = mni.result.f06(fullfile(run_folder,'modal_analysis.f06')).read_modes;
%     %% apply deformation result
%     %pick which mode to plot (1->6)
%     modeshape_num = 5;
%     
%     % apply the modal deformation
%     [~,i] = ismember(model.GRID.GID,res_modeshape.GID(modeshape_num,:));
%     model.GRID.Deformation = [res_modeshape.T1(modeshape_num,i);...
%         res_modeshape.T2(modeshape_num,i);res_modeshape.T3(modeshape_num,i)];
%     
%     % animate the plot to show the mode shape. 'scale' is the scale of the
%     % displacement. 'Frequency' is the frequency.
%     %
%     % PRESS SPACE TO END THE ANIMATION
%     %
%     model.update()
%     model.animate('Frequency',2,'Scale',30)
%     
    
    
    
    
