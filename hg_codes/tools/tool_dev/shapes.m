t1 = 2;
t2 = 1;
Nsec = 25;
yend = 18;
x = linspace(0,yend,Nsec);
xcr=5;

figure;hold on;
plot([0,yend],[t1,t2],'ro');

% % Linear
% A_l = (t2-t1)/yend;
% B_l = t1*ones(size(x));
% y_1 = A_l*x + B_l;
% plot(x,y_1,'b-');
% 
% % Quadratic
% A_q = 0.001;
% B_q = A_l - A_q*yend;
% C_q = B_l;
% y_2 = A_q*x.^2 + B_q*x + C_q; 
% plot(x,y_2,'g-');
% 
% % Cubic
% A_c = 0.00009;
% B_c = -0.004;
% C_c = A_l - A_c*yend^2 - B_c*yend;
% D_c = B_l;
% y_3 = A_c*x.^3 + B_c*x.^2 + C_c*x + D_c; 
% plot(x,y_3,'k-');
% 
% 
% % Quadratic 2
% A_q2=4*A_l/yend;
% B_q2=-4*A_l;
% C_q2=A_l*yend+B_l;
% y_4=A_q2*x.^2 + B_q2*x + C_q2; 
% plot(x,y_4,'r-');


% % test
% yt=0.1*y_1+0.4*y_2+0.2*y_3+0.3*y_4;
% 
% plot(x,yt,'b--')

%% long shape


% conshape=ones(1,20);
% linearshape=1-(0.05:0.05:1);
% 
% mlinear=-1*linearshape+1;%another linear shape
% 
% % quadraticshape=((0.1:0.1:1)-1).*(2*(0.1:0.1:1)-1);
% quadraticshape=-1*(0.05:0.05:1).*(0.05:0.05:1)+1;%quadratic shape
% mquadratic=-1*quadraticshape+1;%another quadratic shape
% 
% for iiii=1:10000
% shape1=xt12(iiii,1)*conshape+xt12(iiii,2)*linearshape+xt12(iiii,3)*quadraticshape+xt12(iiii,4)*mlinear+xt12(iiii,5)*mquadratic;
% 
%  
% 
% shape2=xt12(iiii,6)*conshape+xt12(iiii,7)*linearshape+xt12(iiii,8)*quadraticshape+xt12(iiii,9)*mlinear+xt12(iiii,10)*mquadratic;

%% my shapes

% inboard - flat
[~,ind_i]=find(x<xcr);
[~,ind_o]=find(x>xcr);

y_in_f=ones(1,length(x));
y_in_f(ind_i)=y_in_f(ind_i)*t1;
y_in_f(ind_o)=y_in_f(ind_o)*t2;

plot(x,y_in_f,'r-')

% inboard - linear increase

y_in_l=ones(1,length(x));
A_in_l=(t1-t2)/x(ind_i(end));

y_in_l(ind_i)=x(ind_i)*A_in_l+t2;
y_in_l(ind_o)=t2;

plot(x,y_in_l,'b-')

% inboard - cubic increase 

y_in_c=ones(1,length(x));
B_i_c=3*(t1-t2)/x(ind_i(end))^2;
A_i_c=-2*B_i_c/(3*x(ind_i(end)));

C_i_c=0;
D_i_c=t2;

y_in_c(ind_i)=A_i_c*x(ind_i).^3 + B_i_c*x(ind_i).^2 + C_i_c*x(ind_i) + D_i_c;
y_in_c(ind_o)=t2;

plot(x,y_in_c,'k-')

% outboard - linear drop 

x1=x(ind_o(1));
x2=x(ind_o(end));

A_o_l=(t2-t1)/(x2-x1);
B_o_l=t2-A_o_l*x2;

y_out_l(ind_i)=t2;
y_out_l(ind_o)=A_o_l*x(ind_o)+B_o_l;

plot(x,y_out_l,'k-')

% outboard - quadratic drop 
A_o_q=(t2-t1)/(x2-x1)^2;
B_o_q=t1;

y_out_q(ind_i)=t2;
y_out_q(ind_o)=A_o_q*(x(ind_o)-x1).^2+B_o_q;

plot(x,y_out_q,'k-')

