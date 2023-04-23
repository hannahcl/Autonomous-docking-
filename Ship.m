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
        k_h
 
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
         %obj.A = -inv(obj.M)*obj.N;
         obj.A = -ones(3,3); 
         obj.B = eye(3); 
         obj.C = eye(3);

         %Nominal controler with feedforward
         obj.Q = 10*eye(3); 
         obj.R = eye(3); 
         obj.K = lqr(obj.A, obj.B, obj.Q, obj.R); 
         obj.Kr = obj.B*(inv(obj.C*inv(obj.B*obj.K-obj.A)*obj.B)); 

         %CBF controller
         obj.alpha = 1;   
%         obj.H = [1 0 0]; 
%         obj.k_h = 0; 
         obj.H = eye(3); 
         obj.k_h = -1*[1; 1; 1]; 

         %paramters
         obj.x0 = [-3; 2; -3];
         obj.r = [-1; -1; 0]; 
         
      end

      function u = safe_ctrl_analytical(obj, x)

          u_lower_bound = (obj.H*obj.B)\(-obj.H*obj.A*x - obj.alpha*(obj.H*x+obj.k_h)); 
          u = obj.nominell_controller(x);

%           if (u(1) < u_lower_bound(1)) || (u(2) < u_lower_bound(2)) || (u(3) < u_lower_bound(3))
%             u = u_lower_bound; 
%             disp('On boundary of safe set.')
%           end

          if (u(1) < u_lower_bound(1)) 
            u(1) = u_lower_bound(1); 
            disp('On boundary of safe set, for x1.')
          end

          if (u(2) < u_lower_bound(2)) 
            u(2) = u_lower_bound(2); 
            disp('On boundary of safe set, for x2.')
          end

          if (u(3) < u_lower_bound(3)) 
            u(3) = u_lower_bound(3); 
            disp('On boundary of safe set, for x3.')
          end

      end

      function [c, ceq] = check_barrier_func(obj, x, u)
        
         K_cbf_left = obj.H*obj.A*x + obj.H*obj.B*u; 
         K_cbf_right = -obj.alpha*(obj.H*x + obj.k_h); 

         c = K_cbf_right - K_cbf_left; 
         ceq = 0;   
      end

     function u_safe = safe_ctrl_fmincon(obj, x)
        f = @(u) norm(u - obj.nominell_controller(x))^2;
        nonlcon = @(u) obj.check_barrier_func(x, u);

        u0 = obj.nominell_controller(x);  
        
        options = optimoptions(@fmincon,'Display','none');
        u_safe = fmincon(f, u0, [], [], [], [], [], [], nonlcon, options);
     end

      function u = nominell_controller(obj, x)
        u = obj.Kr*obj.r -obj.K*x; 
      end

      function x_dot = model(obj, x, u)
        x_dot = obj.A*x + obj.B*u;
      end

      function x_dot = closed_loop_model_nomiell(obj, x)
          u = obj.nominell_controller(x); 
          x_dot = obj.model(x, u); 
      end

      function x_dot = closed_loop_model_cbf(obj, x)
          %u = obj.safe_ctrl_analytical(x); 
          u = obj.safe_ctrl_fmincon(x); 
          x_dot = obj.model(x, u); 
      end

   end

  methods

      function sim(obj)
          
        T = 10; 
        ts = 0.2; 
     
        a = [0, 0.5, 0.5, 1]; 
        b = [1/6, 1/3, 1/3, 1/6]; 

        stages = size(b); stages = stages(2); 
        sim_steps = floor(T/ts); 
    
        %% Run simulation without cbf
        x = zeros(3, sim_steps);  
        x(:, 1) = obj.x0; 
        u = zeros(3, sim_steps);
    
        for k = 1:sim_steps 
            store = zeros(3,stages + 1); 
            sum_b = zeros(3,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_nomiell(x(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            x(:, k+1) = x(:, k) + ts*sum_b; 
            u(:, k+1) = obj.nominell_controller(x(:, k)); 
    
        end

        %% Run simulation with cbf
        x_cbf = zeros(3, sim_steps);
        u_cbf = zeros(3, sim_steps);
        x_cbf(:, 1) = obj.x0;

        for k = 1:sim_steps 
            store = zeros(3,stages + 1); 
            sum_b = zeros(3,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_cbf(x_cbf(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            x_cbf(:, k+1) = x_cbf(:, k) + ts*sum_b; 
            u_cbf(:, k+1) = obj.safe_ctrl_fmincon(x_cbf(:, k)); 
    
        end

        %% Plot

        t_arr = 0:ts:T; 

        subplot(2,2,1);
        plot(t_arr,x(1, :), '-o')  
        hold on
        plot(t_arr,x(2, :), '-o')
        plot(t_arr,x(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('State Variables with nominal ctrl')
        legend('x_1', 'x_2', 'x_3')

        subplot(2,2,2);
        plot(t_arr,u(1, :), '-o')  
        hold on
        plot(t_arr,u(2, :), '-o')
        plot(t_arr,u(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('Nominal ctrl input')
        legend('u_1', 'u_2', 'u_3')
       
        subplot(2,2,3);
        plot(t_arr,x_cbf(1, :), '-o')  
        hold on
        plot(t_arr,x_cbf(2, :), '-o')
        plot(t_arr,x_cbf(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('State Variables with CBF')
        legend('x_1', 'x_2', 'x_3')

        subplot(2,2,4);
        plot(t_arr,u_cbf(1, :), '-o')  
        hold on
        plot(t_arr,u_cbf(2, :), '-o')
        plot(t_arr,u_cbf(3, :), '-o')
        hold off
        xlabel('Time')
        ylabel('Ctrl input when using CBF')
        legend('u_cbf_1', 'u_cbf_2', 'u_cbf_3')

        sgtitle('Tilte')
       
      end
   end
end