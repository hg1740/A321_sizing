
% 1-22 thickness 1 spar
% 23-44 thickness 2 skin
% 45-66 Astrg
% initial values

%x=[0.0281069740000000,0.0296910682200000,0.0315902680000000,0.0329900750000000,0.0349788989834126,0.0407432015100000,0.0522362440000000,0.0533170687296945,0.0534612050000000,0.0520904050000000,0.0513023855039625,0.0505684787078318,0.0510259603846507,0.0491945770800000,0.0489910261500000,0.0478704126867382,0.0454181202927953,0.0391673700392247,0.0318407601757024,0.0206794977833985,0.00723994082570000,0.00209542828364976,0.00833284700000000,0.00845500000000000,0.00857430500000000,0.00882753300000000,0.00895443626617861,0.00926972442000000,0.00956242700000000,0.00930490427905254,0.00902012700000000,0.00873620000000000,0.00831744168735230,0.00803746068790621,0.00756715136396738,0.00712025523000000,0.00666034974000000,0.00607674121913327,0.00545812128662870,0.00472215184406769,0.00381994440403732,0.00268435630828117,0.00146771047060000,0.000743420672046869,8.44000000000000e-05,9.25630315236000e-05,9.63000000000000e-05,0.000101983000000000,0.000111430814664192,0.000129693960000000,0.000140640000000000,0.000122837192903585,0.000109740000000000,9.59310000000000e-05,8.41617959003195e-05,7.48748576082306e-05,6.79591281144757e-05,6.02910000000000e-05,5.13572400000000e-05,4.14703900706295e-05,3.28192862601615e-05,2.34623919146144e-05,1.44824055083883e-05,6.13746939102839e-06,2.67391268159400e-06,6.84812073889158e-07];
x=[0.01*ones(1,28),0.001*ones(1,28),1e-5*ones(1,28)];
%initial condition 
cond_set1=[1.1,1.1,1.1,1.1];
cond_set2=[1.1,1.1,1.1,1.1];

Numsec=28;

record=zeros(100,3*Numsec);
counter=1;

while max(cond_set1)>1 || min(cond_set2)<0.95
    
    increase_co=1.08;
    decrease_co=0.96;
    record(counter,:)=x(1,:);
    
    [RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=Static_stress_calc_v1(x);
    
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
    
    x(1,Numsec+skin_inc)=x(1,Numsec+skin_inc)*increase_co;
    x(1,Numsec+skin_dec)=x(1,Numsec+skin_dec)*decrease_co;
    
    % spar check 
    Check_von_spar=RVon_spr/5e8;
   
    [~,sp_index1]=find(Check_von_spar>1);
    [~,sp_index2]=find(Check_von_spar<0.95);
       
    spar_inc=sp_index1;
    spar_dec=sp_index2;
    
%     sp_inc_co=Check_von_spar(spar_inc);
%     sp_dec_co=Check_von_spar(spar_dec);
    
    x(1,spar_inc)=x(1,spar_inc)*increase_co;
    x(1,spar_dec)=x(1,spar_dec)*decrease_co;
    
    % stringers check 
    Check_strg=Rsigma_strg./Rsigma_col;
    
    [~,sg_index1]=find(Check_strg>1);
    [~,sg_index2]=find(Check_strg<0.95);
    
    strg_inc=sg_index1;
    strg_dec=sg_index2;
    
%     strg_inc_co=Check_strg(strg_inc);
%     strg_dec_co=Check_strg(strg_dec);
    
    x(1,Numsec*2+strg_inc)=x(1,Numsec*2+strg_inc)*increase_co;
    x(1,Numsec*2+strg_dec)=x(1,Numsec*2+strg_dec)*decrease_co;
 
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

run_folder = [
        'D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\sizing_2504_extend']; %[-], folder for exporting the NASTRAN model

NastranMethods1 = awi.methods.Nastran;
F061=NastranMethods1.extractNastranResults(strcat(run_folder,'\A320_half_model_SOL144.f06'),'ReadF06',true,'ReadHDF5',false);

% extract forces
M_P1=[F061.f06data.Bendingmoment.UMPLN1(69:95),F061.f06data.Bendingmoment.LMPLN1(95)]; % in-plane moment
M_P2=[F061.f06data.Bendingmoment.UMPLN2(69:95),F061.f06data.Bendingmoment.LMPLN2(95)]; % out of plane moment
T=[F061.f06data.Bendingmoment.UTORQUE1(69:95),F061.f06data.Bendingmoment.LTORQUE1(95)];% torque

S_P1=[F061.f06data.Bendingmoment.USPLN1(69:95),F061.f06data.Bendingmoment.LSPLN1(95)*0.1]; % in plane shear
S_P2=[F061.f06data.Bendingmoment.USPLN2(69:95),F061.f06data.Bendingmoment.LSPLN2(95)*0.1]; % out of plane shear


Grid_coord = h5read(strcat(run_folder,'\A320_half_model_SOL144.h5'),'/NASTRAN/INPUT/NODE/GRID');
Displacement = h5read(strcat(run_folder,'\A320_half_model_SOL144.h5'),'/NASTRAN/RESULT/NODAL/DISPLACEMENT');

Y=Grid_coord.X(2,346:373);
Displacement_Z=Displacement.Z(346:373);

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
plot(Y, x(1:28),'bs','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Spar thickness (m)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

figure % thickness1: skin
plot(Y, x(29:56),'bs','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Skin thickness (m)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')

figure % area: stringers
plot(Y, x(57:84),'bs','MarkerFaceColor','b')
xlabel('Span distance (m)','FontSize',12,'Interpreter','latex')
ylabel('Stringers Area (m$^2$)','FontSize',12,'Interpreter','latex')
set(gcf,'color','w')


% calculate bending stiffness distribution 
thickness1=x(1:28);
thickness2=x(29:56);
Astrg=x(57:84);
Bheight=[0.711000000000000,0.663706105896848,0.616412211793696,0.569118317690545,0.521824423587393,0.474530529484241,0.427236635381089,0.379942741277937,0.368953368744175,0.357963996210414,0.346974623676652,0.335985251142890,0.324995878609128,0.314006506075366,0.303017133541604,0.292027761007842,0.281038388474080,0.270049015940318,0.259059643406556,0.248070270872794,0.237080898339032,0.226091525805270,0.215102153271508,0.204112780737746,0.193123408203984,0.182134035670222,0.171144663136460,0.160155290602698];
Bwidth=[3,2.85770248740050,2.71540497480099,2.57310746220149,2.43080994960198,2.28851243700248,2.14621492440297,2.00391741180347,1.94979613574802,1.89567485969257,1.84155358363712,1.78743230758167,1.73331103152622,1.67918975547077,1.62506847941532,1.57094720335987,1.51682592730442,1.46270465124897,1.40858337519352,1.35446209913807,1.30034082308262,1.24621954702717,1.19209827097172,1.13797699491627,1.08385571886082,1.02973444280537,0.975613166749917,0.921491890694467];


d_strg=sqrt(Astrg/0.36);
t_strg=0.12*d_strg;

%intialise data array
A_val=zeros(1,numel(thickness1));
Ixx_val=zeros(1,numel(thickness1));
Izz_val=zeros(1,numel(thickness1));
J_val=zeros(1,numel(thickness1));
strg_n=0.24;

for ii=1:Numsec
    
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
