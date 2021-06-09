
% 1-25 thickness 1 spar
% 26-50 thickness 2 skin
% 51-75 Astrg
% initial values
x=[0.0316,0.0318,0.0330,0.0318,0.033,0.034,0.037,0.0388,0.039,0.038,0.04,...
    0.0390,0.04055,0.039,0.04,0.0385,0.036,0.037,0.035,0.033,0.0287,0.02,0.01...
    0.01966,0.0278,0.0089,0.0089,0.008892,0.008892,0.008828,0.0087,0.00869...
    0.0085,0.00826,0.0080,0.00779,0.0075,0.0072437,0.00688,0.00658...
    0.006118,0.00568,0.005,0.00456,0.00388,0.0029,0.00123,0.002192,0.0024293,0.0010457,2.48148e-06...
    1.09216e-06,1.1e-06,2.87e-06,1.27e-06,1.977e-06,3.76e-06,1.26e-06,3.55e-06...
    1.54e-06,1.545e-06,1.037e-06,3.3e-06,1.79e-06,3.476e-06,3.5e-06...
    3.1377e-06,1.8923e-06,1.1496e-06,1.3319e-06,9.57e-07,9.57e-07,2.813914e-07,6.44241e-07,1.58e-06];

%initial condition 
cond_set1=[1.1,1.1,1.1,1.1];
cond_set2=[1.1,1.1,1.1,1.1];

record=zeros(100,75);
counter=1;

while max(cond_set1)>1 || min(cond_set2)<0.96
    
    
    record(counter,:)=x(1,:);
    
    [RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=BeamStress_calc_v0(x);
    
    % skin check 
    Check_von_skin1=RVon_skn/5e8;
    Check_von_skin2=Rsigma_pp./Rsigmab_skn;
    
    [~,index1]=find(Check_von_skin1>1);
    [~,index2]=find(Check_von_skin1<0.95);
    [~,index3]=find(Check_von_skin2>1);
    [~,index4]=find(Check_von_skin2<0.95);
    
    skin_inc=unique([index1,index3]); % one of the condition suggest increase then need to increase t
    
    skin_dec=intersect(index2,index4); % only both condition suggest to decrease then decrease
    
    % find the max. coeff. to increase 
%     increase_coeff1_=Check_von_skin1(skin_inc);
%     increase_coeff2_=Check_von_skin2(skin_inc);
%     increase_coeff=max([increase_coeff1_; increase_coeff2_]);
     
%     % find max. to decrease, decrease slowly 
%     
%     decrease_coeff1_=Check_von_skin1(skin_dec);
%     decrease_coeff2_=Check_von_skin2(skin_dec);
%     decrease_coeff=max([decrease_coeff1_; decrease_coeff2_]);
    
    x(1,25+skin_inc)=x(1,25+skin_inc)*1.1;
    x(1,25+skin_dec)=x(1,25+skin_dec)*0.95;
    
    % spar check 
    Check_von_spar=RVon_spr/5e8;
   
    [~,sp_index1]=find(Check_von_spar>1);
    [~,sp_index2]=find(Check_von_spar<0.95);
       
    spar_inc=sp_index1;
    spar_dec=sp_index2;
    
    sp_inc_co=Check_von_spar(spar_inc);
    sp_dec_co=Check_von_spar(spar_dec);
    
    x(1,spar_inc)=x(1,spar_inc)*1.1;
    x(1,spar_dec)=x(1,spar_dec)*0.95;
    
    % stringers check 
    Check_strg=Rsigma_strg./Rsigma_col;
    
    [~,sg_index1]=find(Check_strg>1);
    [~,sg_index2]=find(Check_strg<0.95);
    
    strg_inc=sg_index1;
    strg_dec=sg_index2;
    
%     strg_inc_co=Check_strg(strg_inc);
%     strg_dec_co=Check_strg(strg_dec);
    
    x(1,50+strg_inc)=x(1,50+strg_inc)*1.1;
    x(1,50+strg_dec)=x(1,50+strg_dec)*0.95;
 
    % end condition check 
    cond1=max(Check_von_skin1);
    cond2=max(Check_von_skin2);
    cond3=max(Check_von_spar);
    cond4=max(Check_strg);
    
    cond5=min(Check_von_skin1);
    cond6=min(Check_von_skin2);
    cond7=min(Check_von_spar);
    cond8=min(Check_strg);
    
    cond_set1=[cond1,cond2,cond3,cond4];
    cond_set2=[cond6,cond7,cond8];
    
    disp(cond_set1);
    disp(cond_set2);
    
    counter=counter+1;
    
    
end


%% plotting 

NastranMethods1 = awi.methods.Nastran;
F061=NastranMethods1.extractNastranResults('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test12\A320_half_model_SOL144.f06','ReadF06',true,'ReadHDF5',false);

% extract forces
M_P1=[F061.f06data.Bendingmoment.UMPLN1(69:92),F061.f06data.Bendingmoment.LMPLN1(92)]; % in-plane moment
M_P2=[F061.f06data.Bendingmoment.UMPLN2(69:92),F061.f06data.Bendingmoment.LMPLN2(92)]; % out of plane moment
T=[F061.f06data.Bendingmoment.UTORQUE1(69:92),F061.f06data.Bendingmoment.LTORQUE1(92)];% torque

S_P1=[F061.f06data.Bendingmoment.USPLN1(69:92),F061.f06data.Bendingmoment.LSPLN1(92)]; % in plane shear
S_P2=[F061.f06data.Bendingmoment.USPLN2(69:92),F061.f06data.Bendingmoment.LSPLN2(92)]; % out of plane shear


data = h5read('D:\MATLAB_workspace\ALENA-master\ALENA-master\hg_codes\Sizing_analysis\Result\test10\A320_half_model_SOL144.h5','/NASTRAN/INPUT/NODE/GRID');
Y=data.X(2,346:370);

figure % bending moment
plot(Y,M_P2,'b-s','MarkerFaceColor','b')

xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Bending moment (Nm)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

figure %shear force
plot(Y, S_P2,'b-s','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Vertical shear force (N)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

figure % torque
plot(Y, T,'b-s','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Torque (Nm)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

figure % thickness1: spar
plot(Y, x(1:25),'bs','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Spar thickness (m)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

figure % thickness1: skin
plot(Y, x(26:50),'bs','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Skin thickness (m)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

figure % area: stringers
plot(Y, x(51:75),'bs','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Stringers Area (m$^2$)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')


% calculate bending stiffness distribution 
thickness1=x(1:25);
thickness2=x(26:50);
Astrg=x(51:75);
Bheight=[0.616,0.596,0.577,0.56,0.54,0.52,0.50,0.48,0.46,0.44,0.42,0.398,0.38,0.36,0.34,0.318,0.298,0.278,0.258,0.238,0.218,0.198,0.179,0.1586,0.139];
Bwidth=[3,2.83, 2.66, 2.5, 2.32, 2.15, 1.98,1.8,1.76,1.7,1.66,1.6,1.56,1.51,1.46,1.42,1.37,1.32,1.274,1.22,1.171,1.122,1.07,1.024,0.975];
d_strg=sqrt(Astrg/0.36);
t_strg=0.12*d_strg;

%intialise data array
A_val=zeros(1,numel(thickness1));
Ixx_val=zeros(1,numel(thickness1));
Izz_val=zeros(1,numel(thickness1));
J_val=zeros(1,numel(thickness1));
NumSec=25;
strg_n=0.24;

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
    
    Istrg_xx=((ts*ds^3/12)+(ds*ts^3/12 + ds*ts*(ds/2)^2)*2 + 3*ds*ts*(hs/2-ds/2))*NumStrg*2;
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
    
end

figure % area: stringers
plot(Y, Ixx_val*70e9,'b-s','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('EI (Nm$^2$)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

% figure % area: stringers
% plot(Y, Ixx_val,'b-s','MarkerFaceColor','b')
% xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
% ylabel('EI (Nm$^2$)','FontSize',12,'Interpreter','latex')
% set(gcf,'color','w')
