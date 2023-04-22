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
        Kr

        alpha
        H
 
        x0 
        r
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
         obj.C = eye(3);

         %Nominal controler with feedforward
         obj.Q = eye(3); 
         obj.R = eye(3); 
         obj.K = lqr(obj.A, obj.B, obj.Q, obj.R); 
         obj.Kr = obj.B*(inv(obj.C*inv(obj.B*obj.K-obj.A)*obj.B)); 

         %CBF controller
         obj.alpha = 1; 
%          obj.H = [0 0 0; 
%                 0 -1 0; 
%                 0 0 0];  
        obj.H = eye(3); 

         %paramters
         obj.x0 = [3; 3; -3];
         obj.r = [1; 0; 0]; 
         
      end

      function u_cbf = cbf_controller(obj, x)
        u_cbf = inv(-obj.H*obj.B)*(obj.H*obj.A + obj.alpha*obj.H)*x; 

      end

      function x_dot = model(obj, t, x, u)
        x_dot = obj.A*x + obj.B*u;
      end

      function x_dot = closed_loop_model(obj, t, x, u_cbf)
          x_dot = (obj.A-obj.B*obj.K)*x + obj.Kr*obj.r + obj.B*u_cbf; 
      end

      function x_dot = tot_closed_loop_model(obj, t, x)
          u_cbf = obj.cbf_controller(x); 
          x_dot = obj.closed_loop_model(0, x, u_cbf) 
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
        sim_steps = floor(T/ts); 
    
        x = zeros(3, sim_steps);  
        x(:, 1) = obj.x0; 
        u = zeros(3, sim_steps);
    
        for k = 1:sim_steps 
            store = zeros(3,stages + 1); 
            sum_b = zeros(3,1); 
            for s = 1:stages
                store(:, s+1) = tot_closed_loop_model(obj, 0, x(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            x(:, k+1) = x(:, k) + ts*sum_b; 
            u(:, k+1) = obj.K*x(:, k); 
    
        end

        x_cbf = zeros(3, sim_steps);
        u_cbf = zeros(3, sim_steps);
        x_cbf(:, 1) = obj.x0;

        for k = 1:sim_steps 
            store = zeros(3,stages + 1); 
            sum_b = zeros(3,1); 
            for s = 1:stages
                store(:, s+1) = tot_closed_loop_model(obj, 0, x_cbf(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            x_cbf(:, k+1) = x_cbf(:, k) + ts*sum_b; 
            u_cbf(:, k+1) = obj.cbf_controller(x_cbf(:, k)); 
    
        end

        t_arr = 0:ts:T; 
        
        % Plot 1
        subplot(2,2,1);
        plot(t_arr,x(1, :), '-o')  
        hold on
        plot(t_arr,x(2, :), '-o')
        plot(t_arr,x(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('State Variables')
        legend('x_1', 'x_2', 'x_3')

        subplot(2,2,2);
        plot(t_arr,u(1, :), '-o')  
        hold on
        plot(t_arr,u(2, :), '-o')
        plot(t_arr,u(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('Input')
        legend('u_1', 'u_2', 'u_3')
        
        % Plot 2
        subplot(2,2,3);
        plot(t_arr,x_cbf(1, :), '-o')  
        hold on
        plot(t_arr,x_cbf(2, :), '-o')
        plot(t_arr,x_cbf(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('State Variables')
        legend('x_1', 'x_2', 'x_3')

        subplot(2,2,4);
        plot(t_arr,u_cbf(1, :), '-o')  
        hold on
        plot(t_arr,u_cbf(2, :), '-o')
        plot(t_arr,u_cbf(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('Input')
        legend('u_cbf_1', 'u_cbf_2', 'u_cbf_3')
        
        % Set common title for both plots
        sgtitle('State Variables over Time')
       
      end
   end
end