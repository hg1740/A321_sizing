

% initial value 
t1=0.005; t2=0.005;
sigma_s=10E9; sigma_w=10E9;
sigmac=5E8;

thickness=zeros(100,2);
counter=1;
while max(sigma_w, sigma_s) > sigmac
    
    thickness(counter,:)=[t1,t2];
    
    [sigma1,sigma2]=CalcStress(t1,t2);
    
    sigma_w=sigma1;
    sigma_s=sigma2;
    
    if sigma_w >= sigmac
        
        t1=t1*1.1;
        
    elseif sigma_s >= sigmac
        
        t2=t2*1.1;
    end
    
    counter=counter+1;
    
    disp([sigma_w,sigma_s])
    
    
end

%% plot results 
% NastranMethods1 = awi.methods.Nastran;
F061=NastranMethods1.extractNastranResults('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144.f06','ReadF06',true,'ReadHDF5',false);
% 
data = h5read('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test1\A320_half_model_SOL144.h5','/NASTRAN/INPUT/NODE/GRID');
Y=data.X(2,346:369);
% %including REB2 nodes
% % Wing_nodes=FEM_half.Children.Children(1, 1).Nodes;
% % 
% % Coords=zeros(3,numel(Wing_nodes));
% % 
% % for i=1:numel(Wing_nodes)
% %     Coords(:,i)=Wing_nodes(i).X;  
% % end
% % %  Coords(:,1)=Wing_nodes(1).X
% % Y=Coords(2,1:24);
%  
figure % bending moment
plot(Y,F061.f06data.Bendingmoment.UMPLN2(69:92),'b-s')

xlabel('Span distance (m)','Interpreter','latex')
ylabel('Bending moment (Nm)','Interpreter','latex')
set(gcf,'color','w')

figure %shear force
plot(Y, F061.f06data.Bendingmoment.USPLN2(69:92),'b-s')
xlabel('Span distance (m)','Interpreter','latex')
ylabel('Vertical shear force (N)','Interpreter','latex')
set(gcf,'color','w')

figure % torque
plot(Y, F061.f06data.Bendingmoment.UTORQUE1(69:92),'b-s')
xlabel('Span distance (m)','Interpreter','latex')
ylabel('Torque (Nm)','Interpreter','latex')
set(gcf,'color','w')

disp(F061.f06data.TrimAngle.AoA)
