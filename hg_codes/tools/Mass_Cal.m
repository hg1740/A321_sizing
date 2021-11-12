function Mass = Mass_Cal(Param)

% density
rho=2810;

if isfield(Param,'FWT')
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta, FWT_R]=Aircraft_Models_v3(Param);
    
else
    [Aircraft, FEM_full, Wingbox_right, Box_dimensions, Box_CrossSec, Y_eta]=Aircraft_Models_v3(Param);
end


% cross sectional areas
Area=Box_CrossSec.A;
A1=Area(1:end-1);
A2=Area(2:end);

% Am=(A1+A2)/2;


% segment lengths

if isfield(Param,'FWT')
    
    Y_fwt=Y_eta(end) + linspace(FWT_R.RData(1),FWT_R.RData(2),11);
    
    Y=[Y_eta,Y_fwt(2:end)];
    
else
    
    Y=Y_eta;
    
end
    

% Y=Param.Y';

Y1=Y(1:end-1);
Y2=Y(2:end);
Seg_len=(Y2-Y1);



% Volume
% Volume=sum(Seg_len.*Am);
Volume_=(A1+A2+sqrt(A1.*A2)).*Seg_len/3;

Volume=sum(Volume_);

% Structural mass
Mass=Volume*rho;



end

