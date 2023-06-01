
%simulation configuration
nu0 = [0; 0; 0];
eta0 = [-4; -3.2; pi/2]; 
z0 = [eta0; nu0];
tau0 = zeros(3,1); 

sim_time = 50; 
sim_timestep = 0.1; 

ship = ShipControl(z0, tau0);

ship.dyn.nu_current = [-0.05; -0.01; 0]; 

ship.eta_measurement_variance = 0.1*ones(3,1);
ship.nu_measurement_variance = 0.1*ones(3,1); 
ship.dyn.eta_model_distrubance = 1*ones(3,1);  
ship.dyn.nu_model_distrubance = 1*ones(3,1);  

ship.simRK4( ...
    sim_time, ...
    sim_timestep); 

