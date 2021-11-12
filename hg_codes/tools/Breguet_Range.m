function Range=Breguet_Range(Aspect_ratio, Wing_mass, Wing_area, Payload, Fuel_mass)

% Mass without wing (kg)
Mass0=40000;

% Take off weight
TOW= Mass0 + Wing_mass + Payload + Fuel_mass;

% dynamic pressure 
rho=0.4;
U=220;
q=0.5*rho*U^2;

% Lift coefficient
TOW0 = 48000 + Payload + Fuel_mass;
Cl = TOW0*9.81/(q*Wing_area);

% Drag coefficient
CD0=0.025;
CDi=Cl^2/(pi*Aspect_ratio);

CD=CD0+CDi;

% L/D
L_D=Cl/CD;
% L_D=22;

% Range 
SFC=0.6;
Range=(U*3.6/SFC)*L_D*log(TOW/(TOW - Fuel_mass)); 


end