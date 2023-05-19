
ship = ShipControll(); 

%initial values 
nu0 = [0; 0; 0];
eta0 = [-4; -4; pi/8]; 
z0 = [eta0; nu0];

sim_time = 10; 
sim_timestep = 0.2; 

ship.simRK4(z0, sim_time, sim_timestep); 