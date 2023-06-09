classdef ShipDynamics

   properties
       %nonlinear model, constant matricies

        
        M_RB
        MA
        M 
        D

        %linear model


        C_RB_lin
        C_A_lin
        D_lin
        N_lin
        A_nu
        B_nu
        C_nu

        %symbolic model

        f_symbolic
        g_symbolic

        %disturbance

        tau_wind
        tau_wave
        nu_current

        eta_model_distrubance
        nu_model_distrubance

        %model coefficients
        m = 75;
        Iz = 22; 
        U = 2; 
        x_g = 0;

        X_u_dot = -5; 
        Y_u_dot = -20; 
        X_v_dot = -10; 
        Y_v_dot = -30;
        Y_r_dot = -10; 

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

   end
   methods(Access = public)
       function obj = ShipDynamics()

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
      
         obj.tau_wave = zeros(3,1); 
         obj.tau_wind = zeros(3,1); 
         obj.nu_current = zeros(3,1);

         obj.eta_model_distrubance = zeros(3,1); 
         obj.nu_model_distrubance = zeros(3,1); 

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

        nu_relative = nu - obj.nu_current; 
        C_RB = obj.compute_C_RB(nu);
        N = obj.compute_N(nu_relative);

        f = [
            obj.compute_R(z(3))*z(4:6); 
            -obj.M\(C_RB*nu + N*nu_relative)]; 
        
        g = zeros(6,3); 
        g(4:6, 1:3) = obj.B_nu; 

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
