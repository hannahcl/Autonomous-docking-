
ship = ShipControll(); 

%simulation configuration
nu0 = [0; 0; 0];
eta0 = [-4; -4; pi/8]; 
z0 = [eta0; nu0];

sim_time = 500; 
sim_timestep = 0.2; 

ship.eta_measurement_variance = [0.01; 0.01; 0.01];
ship.nu_measurement_variance =[0.01; 0.01; 0.01];
ship.dyn.eta_model_distrubance =[0.001; 0.001; 0.001]; 
ship.dyn.nu_model_distrubance = [0.1; 0.1; 0.1]; 

ship.simRK4( ...
    z0, ...
    sim_time, ...
    sim_timestep); 

