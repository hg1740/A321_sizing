function [Fuel_mass,Fuselage_total_mass]= Lumped_masses(Fuel_fraction,Fuel_capacity,Payload,Fuselage_structure_mass)

    % fuel density is assumed to be 840 kg/m^3
    Fuel_mass=0.5*Fuel_fraction*Fuel_capacity*840/1e3; %kg
    
    % masses on Fuselage: structural weight + Payload 
    Fuselage_total_mass = Fuselage_structure_mass + Payload;  

end




