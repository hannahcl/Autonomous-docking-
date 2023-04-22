classdef Ship

   properties
        m = 5; 
        Iz = 1; 
        U = 10; 
        x_g = 1;
        X_u_dot = 1; 
        Y_u_dot = 1; 
        X_u = 1;
        Y_v = 1; 
        Y_r = 1; 
        N_v = 1; 
        N_r = 1; 
        
        M_RB
        MA
        M 
        C_RB
        C_A
        D
        N

        B
        C
        A

        Q
        R
        K

        alpha
 
        x0 
   end
   methods(Access = public)
      function obj = Ship()

         %Matricies defining the linear model
         obj.M_RB = [obj.m 0 0; 0 obj.m obj.m*obj.x_g; 0 obj.m*obj.x_g obj.Iz];
         obj.MA = [obj.X_u_dot 0 0; 0 obj.Y_v obj.Y_r; 0 obj.N_v obj.N_r];
         obj.M = obj.M_RB + obj.MA; 
         obj.C_RB = [0 0 0; 0 0 obj.m*obj.U; 0 0 obj.m*obj.x_g*obj.U];
         obj.C_A = [0 0 0; 0 0 -obj.X_u_dot*obj.U; 0 (obj.X_u_dot-obj.Y_u_dot)*obj.U -obj.Y_u_dot*obj.U];
         obj.D = -[obj.X_u 0 0; 0 obj.Y_v obj.Y_r; 0 obj.N_v obj.N_r];
         obj.N = obj.C_RB + obj.C_A + obj.D;

         %LInear model with normal notation
         obj.A = -inv(obj.M)*obj.N;
         obj.B = eye(3); 
         obj.C = [0 -1 0];

         %Nominal controler with feedforward
         obj.Q = eye(3); 
         obj.R = eye(3); 
         obj.K = lqr(obj.A, obj.B, obj.Q, obj.R); 

         %CBF controller
         obj.alpha = 1; 

         %paramters
         obj.x0 = [10; 10; -10];
         
      end

      function u = controller(obj, x)
        %Linear model. No nominal control objective. Allowalbe set: y < 0. 
%         CB = obj.C*obj.B; 
%         nev = (obj.C*obj.A*x + obj.alpha*obj.C*x); 
%         u = nev/(-CB); 
          u = zeros(3,1); 

      end

      function x_dot = model(obj, t, x, u)
        x_dot = obj.A*x + obj.B*u;
      end

      function x_dot = closed_loop_model(obj, t, x, u)
          x_dot = (obj.A-obj.B*obj.K)*x; 
      end

   end

  methods

      function sim(obj)
          
        %RK simulator
        T = 5; 
        ts = 0.4; 
     
        a = [0, 0.5, 0.5, 1]; 
        b = [1/6, 1/3, 1/3, 1/6]; 

        stages = size(b); stages = stages(2); 
        sim_steps = round(T/ts)-1; 
    
        x = zeros(3, sim_steps); 
        x(:, 1) = obj.x0; 
    
        for k = 1:sim_steps
    
            K = zeros(3,stages + 1); 
            sum_b = zeros(3,1); 
            for s = 1:stages
                u = obj.controller(x(:, k)); 
                K(:, s+1) = obj.closed_loop_model(0, x(:, k) + ts*a(s)*K(:, s), u); 
                sum_b  = sum_b + b(s)*K(:,s+1); 
            end
    
            x(:, k+1) = x(:, k) + ts*sum_b; 
    
        end

        t_arr = 0:ts:T; 
        subplot(1,3,1); 
        plot(t_arr,x(1, :), '-o')  
        hold on
        plot(t_arr,x(2, :), '-o')
        hold on
        plot(t_arr,x(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('State Variables')
        legend('x_1', 'x_2', 'x_3')
       
      end
   end
end