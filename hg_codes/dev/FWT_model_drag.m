%% Lift dis

    fold_angle  = -10;   %[deg],
    flare_angle = 25;   %[deg],
    fold_eta=0.65;
    hinge_stiffness=1e-4;
    
%% Wing configurations for starboard wing
  
    Aspect_ratio=19; % Aspect ratio = 10.172 for A321 model
    
    Total_area=126;         % include two wing surface areas + floor size on the fuselage
    Fuselage_width=4;       % dimeter of the fuselage
  
    Wing_span = sqrt(Aspect_ratio*Total_area);
    BeamLoc = 0.4;          % choose a spar location: 0 --> 1
    Semi_span=(Wing_span-Fuselage_width)/2; % length of one wing: 16m for A321 model
    
    Root_chord =  Total_area/(1.064*Semi_span + 4);
    LE_sweep=27;            % deg
    
    Wing_area = (Total_area - Fuselage_width*Root_chord)/2;
    
    Mid_chord=0.63685*Root_chord;
    Tip_chord=0.2248*Root_chord;
    
    X0=Root_chord; 
    X1=0.27*Semi_span*tan(27*pi/180) + Mid_chord;
    X2=Semi_span*tan(27*pi/180) + Tip_chord;
    
    tan_TE_sweep1=(X1-X0)/(0.27*Semi_span);
    tan_TE_sweep2=(X2-X1)/(0.73*Semi_span);
    
    TE_sweep1=atan(tan_TE_sweep1)*180/pi; % deg
    TE_sweep2=atan(tan_TE_sweep2)*180/pi; % deg
      
 
    Taper_ratio=Tip_chord/Root_chord;
    
    Mean_cord_coefficient=(2/3)*(1+Taper_ratio+Taper_ratio^2)/(1+Taper_ratio);
    
    
    %% obtain wingbox geometric properties 

    Wingbox = awi.model.LiftingSurface;
    Wingbox.Origin=[20,2,0];

    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox.ActiveSet = 'sSet';

    % Num of element
    Wingbox.NumBeamElem = 23;

    %Wing dimensions
    Wingbox.SpanVector  = 'Y';
    Wingbox.Span        = Semi_span;   %34.1/2;
    Wingbox.LESweep     = [LE_sweep, LE_sweep];
    Wingbox.LESweep_eta = [0, 1];
    Wingbox.TESweep     = [TE_sweep1, TE_sweep2, TE_sweep2];
    Wingbox.TESweep_eta = [0, 0.27, 1];
    Wingbox.RootChord   = Root_chord;
    
    build(Wingbox)
    
    
    FWT = insertWingFold(Wingbox, 'FlareAngle', flare_angle, 'FoldAngle', fold_angle,'EtaFold',fold_eta);
    FWT.HingeStiffness = [1e14 1e14 1e14 1e14 hinge_stiffness 1e14];
    
    


res_aeroF = mni.result.f06(strcat('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\AR19_FWT_eta65','/A321_36000ft_1g.f06')).read_aeroF;

% % index for each lifting surfaces: 0.85
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:660; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=911:1010; %tail wing _right

% % index for each lifting surfaces: 0.80
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:640; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=871:990; %tail wing _right

% % index for each lifting surfaces: 0.75
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:620; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=831:970; %tail wing _right

% % index for each lifting surfaces: 0.70
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:610; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=811:970; %tail wing _right

% index for each lifting surfaces: 0.65

idx_w85_bef =411:470; % wing_right_bef_kink
idx_w85_aft =471:590; % wing_right_aft_kink
idx_c85=1:50; %conn_right
idx_fwt85=771:940; %tail wing _right

% % index for each lifting surfaces: 0.60
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:580; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=751:930; %tail wing _right




% right side of AC
lift_wing85_bef=res_aeroF.aeroFz(idx_w85_bef)';
lift_wing85_aft=res_aeroF.aeroFz(idx_w85_aft)';
lift_conn85=res_aeroF.aeroFz(idx_c85)';
lift_fwt85=res_aeroF.aeroFz(idx_fwt85)';

% calculate panel width for each segment: conn + wing_bef_kink +
% wing_aft_kink + FWT

R_chord=Wingbox.Chord(1);
semi_span=Wingbox.Span+FWT.Span;
conn=2;
kink_pos=Wingbox.YData(2);

width_conn=2/5;
width_wing1=kink_pos/6;

wing_sec2_num=ceil((Wingbox.YData(3)-Wingbox.YData(2))/(2.5*Wingbox.Chord(2)/10));
width_wing2=(Wingbox.YData(3)-Wingbox.YData(2))/wing_sec2_num;


fwt_num=ceil(FWT.Span/(2.5*Wingbox.Chord(3)/10));
width_fwt=FWT.Span/fwt_num;


% normalise lift by panel width to obtain lift per unit span
lift_wing=[lift_conn85;lift_wing85_bef;lift_wing85_aft;lift_fwt85];

lift_wing_matrix=reshape(lift_wing, 10, numel(lift_wing)/10);

% abs. value of the lift 
lift_wing_abs=sum(lift_wing_matrix);

lift_wing_matrix(:,1:5)=lift_wing_matrix(:,1:5)/width_conn;
lift_wing_matrix(:,6:11)=lift_wing_matrix(:,6:11)/width_wing1;
lift_wing_matrix(:,12:12+wing_sec2_num-1)=lift_wing_matrix(:,12:12+wing_sec2_num-1)/width_wing2;
lift_wing_matrix(:,12+wing_sec2_num:12+wing_sec2_num+fwt_num-1)=lift_wing_matrix(:,12+wing_sec2_num:12+wing_sec2_num+fwt_num-1)/width_fwt;

% lift per unit span
lift_wing_var=sum(lift_wing_matrix);

% find corresponding Y positions
Y_conn=0.5*width_conn:width_conn*0.999:2;

Y_wing1=2+0.5*width_wing1:width_wing1*0.999:2+Wingbox.YData(2);

Y_wing2=2+Wingbox.YData(2)+ 0.5*width_wing2:width_wing2*0.999:2+Wingbox.YData(3);

Y_fwt=2+Wingbox.YData(3)+0.5*width_fwt:width_fwt*0.999:2+Wingbox.Span+FWT.Span;

Y=[Y_conn,Y_wing1,Y_wing2,Y_fwt];

% whole wing span lift distribution 
Y_left=sort(-Y);
lift_wing_var_left=flip(lift_wing_var);
lift_wing_abs_left=flip(lift_wing_abs);

Y_all=[Y_left,Y]';
Lift_all=[lift_wing_var_left,lift_wing_var];
Lift_all_abs=[lift_wing_abs_left,lift_wing_abs];

% calculate vorticity
% Air conditions
FlightPoint=awi.model.FlightPoint;
FlightPoint.Mach=0.78;
FlightPoint.AcVelocity=0.78*340;
FlightPoint.Altitude = 36000;
getFlightPointData(FlightPoint,'ISA');
DynPressure = FlightPoint.DynPressure;

Gamma_all=Lift_all/(FlightPoint.AirDensity*FlightPoint.AcVelocity);

% calculate the derivative 
dGdy= gradient(Gamma_all(:)) ./ gradient(Y_all(:));

% downwash 
wj=zeros(1,length(Y_all));
alphai=zeros(1,length(Y_all));


for i=1:length(Y_all)
    
    w=-(1/(4*pi))*dGdy./(Y_all(i)-Y_all);
    
    w=w( ~any( isnan( w ) | isinf( w ), 2 ),: );
    
    wj(i)=sum(w);
    
    alphai(i)=wj(i)/FlightPoint.AcVelocity;  
    
end

Dragi_var=Lift_all_abs.*sin(alphai);
Drag_force=sum([Dragi_var]);
Cdi=Drag_force/(DynPressure*126);


% % dydx = diff(lift_wing_var(:))./diff(Y(:));
% DYDX= gradient(lift_wing_var(:)) ./ gradient(Y(:));

figure 
Y=Y';
lift_wing_var=lift_wing_var';
plot(Y,lift_wing_var,'s-')


figure 
plot(Y_all,Lift_all,'b.')

% figure 
% plot(Y_all,dGdy,'b.')

figure 
plot(Y_all,wj,'b.')


figure 
plot(Y_all,Dragi_var,'b.')

figure 
plot(Y_all,Lift_all_abs,'b.')


%% lift and drag calculateion   
 
%  % Induced drag and lift coefficient
% 
% Ajj=mni.result.op4(strcat(run_folder,'/ajj.op4')).read_matrix();
% FFaj=mni.result.op4(strcat(run_folder,'/ffaj.op4')).read_matrix();
% 
% res_aeroF = mni.result.f06(strcat(run_folder,'/A321_range_calc.f06')).read_aeroF;
% 
% % Air conditions
% FlightPoint=awi.model.FlightPoint;
% FlightPoint.Mach=TrimLoadcase.Mach;
% FlightPoint.AcVelocity=FlightPoint.Mach*340;
% FlightPoint.Altitude = 36000;
% getFlightPointData(FlightPoint,'ISA');
% 
% DynPressure = FlightPoint.DynPressure;
% 
% % downwash 
% WJ = Ajj*(FFaj./DynPressure);
% 
% % surface_normal = model.fwt_normal_vector();
% % drag_mag = dot([1 0 0],surface_normal);
% 
% 
% % index for each lifting surfaces
% idx =407:958; % wing_right
% idx_c=1:48; %conn_right
% idx_t=97:184; %tail wing _right
% 
% idx1 =959:1510; % wing_left
% idx_c1=49:96; %conn_right
% idx_t1=188:278; %tail wing _left
% 
% % right side of AC
% lift_wing=res_aeroF.aeroFz(idx)';
% lift_conn=res_aeroF.aeroFz(idx_c)';
% lift_tail=res_aeroF.aeroFz(idx_t)';
% 
% % left side of AC
% lift_wing_L=res_aeroF.aeroFz(idx1)';
% lift_conn_L=res_aeroF.aeroFz(idx_c1)';
% lift_tail_L=res_aeroF.aeroFz(idx_t1)';
% 
% total_lift=sum(lift_wing)+sum(lift_conn)+sum(lift_tail)-sum(lift_wing_L)-sum(lift_conn_L)-sum(lift_tail_L);
% 
% total_area=(Wingbox_right.SurfaceArea+Connector_right.SurfaceArea+ Tailwing_right.SurfaceArea)*2;
% 
% % calculation of the lift coefficient: CL 
% Cl=total_lift/(DynPressure*total_area);
% 
% 
% % Induced drag - wing
% 
% % wing total on the wing
% lift_wing_left=[lift_conn_L;lift_wing_L];
% 
% lift_wing_left_matrix=reshape(lift_wing_left, 12, numel(lift_wing_left)/12);
% 
% % abs lift distribution on the span
% abs_lift_left=sum(lift_wing_left_matrix);
% 
% % normalise lift by panel width to obtain lift per unit span
% lift_wing_left_matrix(:,1:4)=lift_wing_left_matrix(:,1:4)/0.5;
% lift_wing_left_matrix(:,5:13)=lift_wing_left_matrix(:,5:13)/0.477;
% lift_wing_left_matrix(:,14:end)=lift_wing_left_matrix(:,14:end)/0.31371;
% 
% % lift per unit span
% lift_wing_left=sum(lift_wing_left_matrix);
% 
% %find Y-coords
% y0=0.25:0.5:1.75;
% y1=2.2385:0.477:6.0545;
% y2=6.449915:0.3137:17.74335;
% Y=[y0,y1,y2];
% 
% % Y VS lift
% Y_all=[-Y,Y];
% 
% %absolute lift
% abs_Lift_wing_all=[-abs_lift_left,-abs_lift_left];
% 
% %lift per unit span 
% Lift_wing_all=[-lift_wing_left,-lift_wing_left];
% 
% % circulation
% Gamma=Lift_wing_all/(FlightPoint.AirDensity*FlightPoint.AcVelocity);
% 
% % curve fitting
% Y_fit=linspace(-17.9,17.9,300);
% 
% P1=polyfit(Y_all,Gamma,10);
% % gamma_fit=polyval(P1,Y1_fit);
% P2 = polyder(P1);
% 
% wj=zeros(numel(Y_all),1);
% Fd=zeros(numel(Y_all),1);
% alphai=zeros(numel(Y_all),1);
% 
% for i=1:1:numel(Y_all)
%     
%     x0=Y_all(i);
%     fun_wj = @(x) (-1/(4*pi))*(P2(1)*x.^9+P2(2)*x.^8+P2(3)*x.^7+P2(4)*x.^6+ P2(5)*x.^5+ P2(6)*x.^4+ P2(7)*x.^3+ P2(8)*x.^2+ P2(9)*x +  P2(10))./(x0-x);
%     
%     wj_l = integral(fun_wj,-17.9, x0-0.00001);
%     wj_r = integral(fun_wj,x0+0.00001, 17.9);
%     
%     wj(i)= wj_l+wj_r;
%     alphai(i)=wj(i)/FlightPoint.AcVelocity;
%     Fd(i)=abs_Lift_wing_all(i)*sin(alphai(i));
%      
% end
% 
% % rm trivial data at the tip 
% %TODO: improve the way of integration/curve fiting
% 
% Fd(Fd>10)=0;
% Drag_force=sum(Fd);
% 
% Cdi=Drag_force/(DynPressure*126);
% 
% figure 
% plot(Y_all,alphai,'b.')
% 
% figure 
% plot(Y_all,Lift_wing_all,'b.')
% set(gcf,'Color','w')
% xlabel('Span (m)','Interpreter','latex','FontSize',12)
% ylabel('Lift per unit span (N/m)','Interpreter','latex','FontSize',12)
% % 
% % figure
% % plot(Y_all,abs_Lift_wing_all,'b.')
% figure 
% plot(WJ([idx_c,idx]),'b.')