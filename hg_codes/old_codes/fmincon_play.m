

%% -----------------------------------------------------------

x0=[2*ones(1,25),3*ones(1,25),3*ones(1,25)];

lb=[1*ones(1,25),2*ones(1,25),3*ones(1,25)];
ub=10.0*ones(1,75);

% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% [xopt, fva,ef,output,lamda]=fmincon(@objective,x0,[],[],[],[],lb,[],@constraint,[]);

x=fmincon(@objective,x0,[],[],[],[],lb,ub,@constraint,[]);


function obj=objective(x)

a=weight_calc(x);
obj=a;


% obj=x(1)*x(4)*(x(1)+x(2)+x(3))+x(3);
% obj=t(1)*t(4)*(t(1)+t(2)+t(3))+t(3);

end


function [c,ceq]=constraint(x)


c=stress_calc(x)-100;

ceq=[];


% t=[x(1),x(2),x(3),x(4)];
% c=25-t(1)*t(2)*t(3);
% 
% ceq=sum(t.^2)-40;

end



function weight=weight_calc(x)

    t1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
        x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
        x(21),x(22),x(23),x(24),x(25)];

    t2=[x(26),x(27),x(28),x(29),x(30),x(31),x(32),x(33),x(34),x(35)...
        x(36),x(37),x(38),x(39),x(40),x(41),x(42),x(43),x(44),x(45)...
        x(46),x(47),x(48),x(49),x(50)];

    A=[x(51),x(52),x(53),x(54),x(55),x(56),x(57),x(58),x(59),x(60)...
        x(61),x(62),x(63),x(64),x(65),x(66),x(67),x(68),x(69),x(70)...
        x(71),x(72),x(73),x(74),x(75)];
    
    f1=0;f2=0;f3=0;
    
    for k=1:25
        
        f1=f1+exp(t1(k))*k;
        f2=f2+exp(t2(k))*k;
        f3=f3+exp(A(k))*k;
        
    end
    
    weight=abs(f1+f2+f3);
    
    
end
    
    
function sigma=stress_calc(x)

    t1=[x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10)...
        x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20)...
        x(21),x(22),x(23),x(24),x(25)];

    t2=[x(26),x(27),x(28),x(29),x(30),x(31),x(32),x(33),x(34),x(35)...
        x(36),x(37),x(38),x(39),x(40),x(41),x(42),x(43),x(44),x(45)...
        x(46),x(47),x(48),x(49),x(50)];

    A=[x(51),x(52),x(53),x(54),x(55),x(56),x(57),x(58),x(59),x(60)...
        x(61),x(62),x(63),x(64),x(65),x(66),x(67),x(68),x(69),x(70)...
        x(71),x(72),x(73),x(74),x(75)];
    
    
    sigma=sum(t1)+sum(t2)+sum(A);




end
    











