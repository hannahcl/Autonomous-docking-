classdef ShipDynamics

   properties
        m = 5; 
        Iz = 1; 
        U = 10; 
        x_g = 1;

        X_u_dot = 1; 
        Y_u_dot = 1; 
        X_v_dot = 1; 
        Y_v_dot = 1; 
        Y_r_dot = 1; 

        X_u = 1;
        Y_v = 1; 
        Y_r = 1; 
        N_v = 1; 
        N_r = 1;

        X_uu = 1; 
        Y_vv = 1; 
        Y_rv = 1; 
        Y_vr = 1; 
        Y_rr = 1; 
        N_vv = 1; 
        N_rv = 1; 
        N_vr = 1; 
        N_rr = 1; 

        A_11 = 1; 
        A_22 = 1; 
        A_26 = 1;
        A_62 = 1; 
        A_66 = 1; 

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

        tau_wind
        tau_wave


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
         
       end

     function R = compute_R(obj, eta)
        R = [cos(eta(3)) -sin(eta(3)) 0; 
             sin(eta(3)) cos(eta(3)) 0; 
             0 0 1];
      end


      function eta_dot = model_eta(obj, eta, nu)
         R = obj.compute_R(eta); 
         eta_dot = R*nu; 
      end

      function nu_dot = model_nu(obj, nu, tau)
        N = obj.compute_N(nu); 
        nu_dot = obj.M\(tau + obj.tau_wave + obj.tau_wind - N*nu); 
      end

      function nu_dot = linear_model_nu(obj, nu, tau) 
        nu_dot = -obj.M\obj.N_lin*nu + tau+ obj.tau_wave + obj.tau_wind;
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

      function Dn = compute_Dn(obj, nu)
        u = abs(nu(1)); 
        v = abs(nu(2)); 
        r = abs(nu(3)); 

        Dn = [
            obj.X_uu*u 0 0; 
            0 obj.Y_vv*v+obj.Y_rv*r obj.Y_vr*v+obj.Y_rr*r; 
            0 obj.N_vv*v+obj.N_rv*r obj.N_vr*v+obj.N_rr*r]; 
      end

      function N = compute_N(obj, nu)
        N = (obj.compute_C_RB(nu) + obj.compute_C_A(nu)) ...
            + obj.D ...
            + obj.compute_Dn(nu); 
      end

   end

end
