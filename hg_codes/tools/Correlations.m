be_t=[13,20,30,40,50];

Fe=[67,62,42,24,15.5];


be_t=-7.743*Fe^4 + 0.0006387*Fe^3 + 0.007084*Fe^2 - 1.966*Fe + 76.83;



% test code skin thickness iteration

sigma_psi=200000;
L_inch=20;

% intial guess
t_s=0.15; %(inch)

delta=1; % (inch)

while delta>0.001
    
    N=sigma_psi*t_s;
    Fe=2000*(N/L_inch)^0.5;
    
    t_sk=N/(Fe*1.5);
    
    delta=abs(t_sk-t_s);
    
    t_s=t_sk;
    
end



N_L=[50,100,200,400,800];


b_b=[4,3.5,2,1.25,1];