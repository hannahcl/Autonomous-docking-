
%simulation configuration
nu0 = [0; 0; 0];
eta0 = [4; 1; -pi/8]; 
z0 = [eta0; nu0];
tau0 = zeros(3,1); 

sim_time = 10; 
sim_timestep = 0.1; 

ship = ShipControll(z0, tau0);

ship.eta_measurement_variance = [0.01; 0.01; 0.01];
ship.nu_measurement_variance =[0.01; 0.01; 0.01];
ship.dyn.eta_model_distrubance =[0.001; 0.001; 0.001]; 
ship.dyn.nu_model_distrubance = [0.01; 0.01; 0.01]; 

ship.simRK4( ...
    sim_time, ...
    sim_timestep); 

