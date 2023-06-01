classdef ShipDynamics

   properties
        m = 75; %[kg]
        Iz = 22; %[kg*m^2] % second moment of interitia = 0.021; %[m^4]
        U = 2; 
        x_g = 0;

        X_u_dot = -5; %[kg/m^2]
        Y_u_dot = -20; %[kg/m^2]
        X_v_dot = -10; %[kg/m^2]
        Y_v_dot = -30; %[kg/m^2]
        Y_r_dot = -10; %[kg·m^2/rad]

        X_u = -1;
        Y_v = -30; 
        Y_r = -7; 
        N_v = -7; 
        N_r = -2;

        X_uu = 1; 
        Y_vv = 5; 
        Y_rv = 0.1; 
        Y_vr = 0.1; 
        Y_rr = 0.1; 
        N_vv = 5; 
        N_rv = 0.1; 
        N_vr = 0.1; 
        N_rr = 1.1; 

        A_11 = 10; 
        A_22 = 10; 
        A_26 = 10;
        A_62 = 0.1; 
        A_66 = 0.1; 

        B_11v = 1; 
        B_22v = 1; 
        B_66v = 1; 
        
        M_RB
        MA
        M 
        D


        C_RB_lin
        C_A_lin
        D_lin
        N_lin
        A_nu
        B_nu
        C_nu

        f_symbolic
        g_symbolic

        tau_wind
        tau_wave
        nu_current

        eta_model_distrubance
        nu_model_distrubance


   end
   methods(Access = public)
       function obj = ShipDynamics()

       %nonlinar model

         obj.M_RB = [
             obj.m 0 0; 
             0 obj.m obj.m*obj.x_g; 
             0 obj.m*obj.x_g obj.Iz];
         
         obj.MA = [
             obj.A_11 0 0; 
             0 obj.A_22 obj.A_26; 
             0 obj.A_62 obj.A_66]; 

         obj.M = obj.M_RB + obj.MA; 

         obj.D = [
             obj.B_11v 0 0; 
             0 obj.B_22v 0; 
             0 0 obj.B_66v]; 

         % disturbance 
      
         obj.tau_wave = zeros(3,1); 
         obj.tau_wind = zeros(3,1); 
         obj.nu_current = zeros(3,1);

         obj.eta_model_distrubance = zeros(3,1); 
         obj.nu_model_distrubance = zeros(3,1); 

         % linear model

         obj.C_RB_lin = [
             0 0 0;
             0 0 obj.m*obj.U;
             0 0 obj.m*obj.x_g*obj.U];

         obj.C_A_lin = [
             0 0 0;
             0 0 -obj.X_u_dot*obj.U;
             0 (obj.X_u_dot-obj.Y_u_dot)*obj.U -obj.Y_u_dot*obj.U];

         obj.D_lin = -[
             obj.X_u 0 0;
             0 obj.Y_v obj.Y_r;
             0 obj.N_v obj.N_r];

         obj.N_lin = obj.C_RB_lin + obj.C_A_lin + obj.D_lin; 

         obj.A_nu = -obj.M\obj.N_lin; 
         obj.B_nu = inv(obj.M); 
         obj.C_nu = eye(3);


         fhs = obj.compute_f_g_symbolic();
         obj.f_symbolic = fhs{1}; 
         obj.g_symbolic = fhs{2}; 
         
       end

     function fhs = compute_f_g_symbolic(obj)
        % z' = f(z) + g(z)*tau, where z = [eta; nu]
        z = sym('z', [6 1]);
        nu = z(4:6); 

        f = 0; 
        
        g = obj.compute_R(z(1)); 

        fhs = {f; g};

     end

     function R = compute_R(obj, psi)
        R = [cos( psi) -sin( psi) 0; 
             sin( psi) cos( psi) 0; 
             0 0 1];
      end


      function eta_dot = model_eta(obj, eta, nu)
         R = obj.compute_R(eta(3)); 
         eta_dot = R*nu + randn(1)*obj.eta_model_distrubance; 
      end

      function nu_dot = model_nu(obj, nu, tau)
        nu_relative = nu - obj.nu_current; 

        C_RB = obj.compute_C_RB(nu);
        N = obj.compute_N(nu_relative); 

        nu_dot = -obj.M\(C_RB*nu + N*nu_relative) + obj.M\(tau + obj.tau_wave + obj.tau_wind) + randn(1)*obj.nu_model_distrubance;  
      end

      function nu_dot = linear_model_nu(obj, nu, tau)  
        nu_dot = -obj.M\obj.N_lin*nu + obj.M\(tau + obj.tau_wave + obj.tau_wind) + randn(1)*obj.nu_model_distrubance;
      end

   end

   methods(Access = private)


      function C_RB = compute_C_RB(obj, nu)
        r = nu(3); 
        C_RB = [0 -obj.m*r -obj.m*obj.x_g*r; 
                obj.m*r 0 0; 
                obj.m*obj.x_g*r 0 0]; 
      end

      function C_A = compute_C_A(obj, nu)
        u = nu(1); 
        v = nu(2); 
        r = nu(3); 
        C_A = [0 0 (obj.Y_v_dot*v + obj.Y_r_dot*r); 
               0 0 -obj.X_u_dot*u; 
               (-obj.Y_v_dot*v-obj.Y_r_dot*r) obj.X_u_dot*u 0]; 
      end

      function N = compute_N(obj, nu)
        CA = obj.compute_C_A(nu); 
        N = CA + obj.D;  
      end

   end

   methods
      function obj = set.tau_wind(obj, val)
         obj.tau_wind = val; 
      end

      function obj = set.tau_wave(obj, val)
         obj.tau_wave = val; 
      end
    
   end

end
