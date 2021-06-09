
% 1-22 thickness 1 spar
% 23-44 thickness 2 skin
% 45-66 Astrg
% initial values

x=[0.02*ones(1,NumSec),0.005*ones(1,NumSec),1e-7*ones(1,NumSec)];

%initial condition 
cond_set1=[1.1,1.1,1.1,1.1];
cond_set2=[1.1,1.1,1.1,1.1];

Numsec=25;

record=zeros(100,75);
counter=1;

while max(cond_set1)>1 || min(cond_set2)<0.95
    
    increase_co=1.1;
    decrease_co=0.95;
    record(counter,:)=x(1,:);
    
    [RVon_skn, RVon_spr, Rsigmab_skn, Rtaub_skn, Rsigma_pp, Rsigma_strg, Rsigma_crip, Rsigma_col]=BeamStress_calc_v2(x);
    
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

NastranMethods1 = awi.methods.Nastran;
F061=NastranMethods1.extractNastranResults(strcat(run_folder,'\A320_half_model_SOL144.f06'),'ReadF06',true,'ReadHDF5',false);

% extract forces
M_P1=[F061.f06data.Bendingmoment.UMPLN1(69:89),F061.f06data.Bendingmoment.LMPLN1(89)]; % in-plane moment
M_P2=[F061.f06data.Bendingmoment.UMPLN2(69:89),F061.f06data.Bendingmoment.LMPLN2(89)]; % out of plane moment
T=[F061.f06data.Bendingmoment.UTORQUE1(69:89),F061.f06data.Bendingmoment.LTORQUE1(89)];% torque

S_P1=[F061.f06data.Bendingmoment.USPLN1(69:89),F061.f06data.Bendingmoment.LSPLN1(89)*0.1]; % in plane shear
S_P2=[F061.f06data.Bendingmoment.USPLN2(69:89),F061.f06data.Bendingmoment.LSPLN2(89)*0.1]; % out of plane shear


data = h5read('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\a321_cruise_2p5g.h5','/NASTRAN/INPUT/NODE/GRID');
Y=data.X(2,346:367);

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
thickness1=x(1:22);
thickness2=x(23:44);
Astrg=x(45:66);
Bheight=[0.711000000000000,0.654695466486852,0.598390932973703,0.542086399460555,0.485781865947406,0.429477332434258,0.373172798921109,0.357676186633189,0.342179574345268,0.326682962057348,0.311186349769428,0.295689737481508,0.280193125193587,0.264696512905667,0.249199900617747,0.233703288329826,0.218206676041906,0.202710063753985,0.187213451466065,0.171716839178145,0.156220226890225,0.140723614602304];
Bwidth=[3,2.82803516079563,2.65607032159126,2.48410548238689,2.31214064318252,2.14017580397815,1.96821096477378,1.89097603915990,1.81374111354603,1.73650618793215,1.65927126231827,1.58203633670439,1.50480141109052,1.42756648547664,1.35033155986276,1.27309663424888,1.19586170863501,1.11862678302113,1.04139185740725,0.964156931793372,0.886922006179495,0.809687080565617];
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

%% Drag calc

% function res = get_drag_moment(fold_angle,twist_angle,flare_angle,origin,root_aoa,V)
% 
%     Ajj = op4ToMatrix('ajj.op4');
% 
%     FFAj = op4ToMatrix('ffaj.op4');
% 
%     res_aeroF = mni.result.f06('sol144.f06').read_aeroF;
% 
%     
% 
%     q = (0.5*1.225*V^2);
% 
%     WJ = Ajj*(FFAj./q);
% 
%     
% 
%     panels = [10,15];
% 
%     % get moment arm for each panel
% 
%     xx = linspace(0,0.333849,panels(2)*2+1);
% 
%     xx = xx(2:2:panels(2)*2)*cosd(flare_angle);
% 
%     mom_arm = reshape(repmat(xx,panels(1),1),length(xx)*panels(1),[]);
% 
%     
% 
%     % get FWT surface normal
% 
%     model = gen.WT_model(fold_angle,twist_angle,flare_angle,origin,root_aoa);
% 
%     surface_normal = model.fwt_normal_vector();
% 
%  
% 
%     %calc total drag moment
% 
%     idx = 401:550;
% 
%     fwt_drag = sin(WJ(idx)).*res_aeroF.aeroFz(idx)';
% 
%     drag_mag = dot([1 0 0],surface_normal);
% 
%     res = sum(fwt_drag.*drag_mag.*mom_arm);
% 
% end
% 
% 
% function [matrix, varargout] = op4ToMatrix(filename, varargin)
% 
% %op4ToMatrix Returns a Matlab matrix based on the ASCII data in the .op4
% 
% %defined in 'filename'.
% 
%  
% 
% if nargin < 1
% 
%     [name, path] = uigetfile({'*.op4', 'OUTPUT4 File'}, 'Select an OUTPUT 4 file (.op4)');
% 
%     if isnumeric(name) || isnumeric(path)
% 
%         return
% 
%     end
% 
%     filename = fullfile(path, name);
% 
% end
% 
%  
% 
% %Parse inputs
% 
% p = inputParser;
% 
% addRequired(p, 'filename', @ischar);
% 
% addOptional(p, 'matrixNames', [], @iscell);
% 
% parse(p, filename, varargin{:});
% 
%  
% 
% %Extract data from the .op4 file as raw text
% 
% data = extractDataFromFile(filename);
% 
% if isempty(data)
% 
%     error('File ''%s'' is empty. Unable to extract any data.', filename);
% 
% end
% 
%  
% 
% %Extract matrix header data
% 
% [nCol, nRow, name, maxEntryPerLine, numFormat] = extractHeaderData(data{1});
% 
% data  = data(2:end); %No longer need the header line!
% 
% nData = numel(data);
% 
%  
% 
% %Preallocate matrix
% 
% matrix = zeros(nRow, nCol);
% 
%  
% 
% %Parse data and populate the matrix
% 
% isMatrixEnd = false;
% 
% index = 1;
% 
% while ~isMatrixEnd
% 
%     %Read entry
% 
%     headerData = sscanf(data{index}, '%f %f %f');
% 
%     colID  = headerData(1);
% 
%     rowID  = headerData(2);
% 
%     nEntry = headerData(3);    
% 
%     %How many lines of data need to be read?
% 
%     nLine  = ceil(nEntry/maxEntryPerLine);
% 
%     %Check index has not exceeded matrix size
% 
%     if colID > nCol || rowID > nRow
% 
%         %Update counter
% 
%         index = index + nLine + 1;
% 
%         %Check if we have exceeded the end of the data
% 
%         if index > nData
% 
%             isMatrixEnd = true;
% 
%         end
% 
%         continue
% 
%     end
% 
%     %Read data
% 
%     rowData = data(index+1 : index + nLine);
% 
%     rowData = sscanf(cat(2, rowData{:}), '%f');
% 
%     %Assign data to the matrix
% 
%     matrix(rowID : rowID + nEntry - 1,colID) = rowData;
% 
%     %Update counter
% 
%     index = index + nLine + 1;
% 
%     %Check if we have exceeded the end of the data
% 
%     if index > nData
% 
%         isMatrixEnd = true;
% 
%     end
% 
% end
% 
%  
% 
% if nargout > 1
% 
%    varargout{1} = [nRow, nCol];
% 
% end
% 
%  
% 
% end
% 
%  
% 
% function data = extractDataFromFile(filename)
% 
% %extractDataFromFile Extracts the data from the file 'filename' and returns
% 
% %the data as a cell array of strings.
% 
%  
% 
% %Get file identifier
% 
% fid = fopen(filename, 'r');
% 
% if fid == -1
% 
%     error('Unable to open file ''%s''.', filename);
% 
% end
% 
%  
% 
% %Import data as a cell array of strings
% 
% data = textscan(fid, '%s', ...
% 
%     'delimiter'    , '\n', ...
% 
%     'CommentStyle' , '$' , ...
% 
%     'whitespace'   , '');
% 
% data = data{:};
% 
%  
% 
% %Close file
% 
% fclose(fid);
% 
%  
% 
% end
% 
%  
% 
% function [nCol, nRow, name, maxEntryPerLine, numFormat] = extractHeaderData(headerLine)
% 
% %extractHeaderData Extracts the meta data for the matrix from the header
% 
% %line (the first line) of the .op4 file.
% 
% %
% 
% % Meta data includes:
% 
% %   - The number of columns (nCol)
% 
% %   - The number of rows    (nRow)
% 
% %   - The name of the matrix (name)
% 
% %   - The max number of entries per line (maxEntryPerLine)
% 
% %   - The format spec for the matrix data as printed in the op4 (numFormat)
% 
%  
% 
% headerData = strsplit(headerLine); % -> split by blankspace
% 
%  
% 
% nCol = str2double(headerData{2});
% 
% nRow = str2double(headerData{3});
% 
%  
% 
% name = headerData{5};
% 
% name = name(isletter(name));
% 
%  
% 
% numFormat = strsplit(headerData{end}, ',');
% 
% numFormat = strsplit(numFormat{2}, 'E');
% 
%  
% 
% maxEntryPerLine = str2double(numFormat{1});
% 
%  
% 
% numFormat = ['%', numFormat{2}, 'e'];
% 
%  
% 
%  
% 
% end


