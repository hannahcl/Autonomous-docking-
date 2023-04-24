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

        Q_lqr
        R_lqr
        K
        Kr

        alpha
        H
        k_h
 
        v0 
        x0
        z0
        ref_v
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

         %LInear model with normal notation for v dynamics
         %obj.A = -inv(obj.M)*obj.N;
         obj.A = -ones(3,3); 
         obj.B = eye(3); 
         obj.C = eye(3);

         %Nominal controler with feedforward
         obj.Q_lqr = 10*eye(3); 
         obj.R_lqr= eye(3); 
         obj.K = lqr(obj.A, obj.B, obj.Q_lqr, obj.R_lqr); 
         obj.Kr = obj.B*(inv(obj.C*inv(obj.B*obj.K-obj.A)*obj.B)); 

         %CBF controller
         obj.alpha = 1;   
%         obj.H = [1 0 0]; 
%         obj.k_h = 0; 
         obj.H = eye(3); 
         obj.k_h = -1*[1; 1; 1]; 

         %paramters
         obj.v0 = [-2; 2; -3];
         obj.x0 = [0; 1; 1]; 
         obj.z0 = [obj.x0; obj.v0];
         obj.ref_v= [0; 0; 0]; 
         
      end

      %% Controllers

      function u = safe_ctrl_analytical(obj, v)

          u_lower_bound = (obj.H*obj.B)\(-obj.H*obj.A*v - obj.alpha*(obj.H*v+obj.k_h)); 
          u = obj.nominell_controller(v);

%           if (u(1) < u_lower_bound(1)) || (u(2) < u_lower_bound(2)) || (u(3) < u_lower_bound(3))
%             u = u_lower_bound; 
%             disp('On boundary of safe set.')
%           end

          if (u(1) < u_lower_bound(1)) 
            u(1) = u_lower_bound(1); 
            disp('On boundary of safe set, for v1.')
          end

          if (u(2) < u_lower_bound(2)) 
            u(2) = u_lower_bound(2); 
            disp('On boundary of safe set, for v2.')
          end

          if (u(3) < u_lower_bound(3)) 
            u(3) = u_lower_bound(3); 
            disp('On boundary of safe set, for v3.')
          end

      end

      function [c, ceq] = check_barrier_func(obj, v, u)
        
         K_cbf_left = obj.H*obj.A*v + obj.H*obj.B*u; 
         K_cbf_right = -obj.alpha*(obj.H*v + obj.k_h); 

         c = K_cbf_right - K_cbf_left; 
         ceq = 0;   
      end

     function u_safe = safe_ctrl_fmincon(obj, v)
        f = @(u) norm(u - obj.nominell_controller(v))^2;
        nonlcon = @(u) obj.check_barrier_func(v, u);

        u0 = obj.nominell_controller(v);  
        
        options = optimoptions(@fmincon,'Display','none');
        u_safe = fmincon(f, u0, [], [], [], [], [], [], nonlcon, options);
     end

      function u = nominell_controller(obj, v)
        u = obj.Kr*obj.ref_v-obj.K*v; 
      end

      %% Models

      function v_dot = model_v_dyn(obj, v, u)
        v_dot = obj.A*v + obj.B*u;
      end

      function x_dot = model_x_dyn(obj, x, v)
         R = [cos(x(3)) -sin(x(3)) 0; 
             sin(x(3)) cos(x(3)) 0; 
             0 0 1]; 

         x_dot = R*v; 
      end

      function v_dot = closed_loop_model_v_dyn_nomiell(obj, v)
          u = obj.nominell_controller(v); 
          v_dot = obj.model_v_dyn(v, u); 
      end

      function v_dot = closed_loop_model_v_dyn_cbf(obj, v)
          %u = obj.safe_ctrl_analytical(v); 
          u = obj.safe_ctrl_fmincon(v); 
          v_dot = obj.model_v_dyn(v, u); 
      end

      function z_dot = closed_loop_model_z_dyn_nomiell(obj, z)
          x = z(1:3); 
          v = z(4:6); 

          u = obj.nominell_controller(v); 
          v_dot = obj.model_v_dyn(v, u); 
          x_dot = obj.model_x_dyn(x, v); 

          z_dot = [x_dot; v_dot]; 
      end

      function z_dot = closed_loop_model_z_dyn_cbf(obj, z)
          x = z(1:3); 
          v = z(4:6); 

          %u = obj.safe_ctrl_analytical(v); 
          u = obj.safe_ctrl_fmincon(v); 
          v_dot = obj.model_v_dyn(v, u); 
          x_dot = obj.model_x_dyn(x, v); 

          z_dot = [x_dot; v_dot]; 
      end



   end

  methods

      function sim(obj)
          
        T = 3; 
        ts = 0.2; 
     
        a = [0, 0.5, 0.5, 1]; 
        b = [1/6, 1/3, 1/3, 1/6]; 

        stages = size(b); stages = stages(2); 
        sim_steps = floor(T/ts); 
    
        %% Run simulation without cbf
        z = zeros(6, sim_steps);  
        z(:, 1) = obj.z0; 
        u = zeros(3, sim_steps);
    
        for k = 1:sim_steps 
            store = zeros(6,stages + 1); 
            sum_b = zeros(6,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_z_dyn_nomiell(z(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            z(:, k+1) = z(:, k) + ts*sum_b; 
            u(:, k+1) = obj.nominell_controller(z(4:6, k)); 
    
        end

        %% Run simulation with cbf
        z_cbf = zeros(6, sim_steps);
        z_cbf(:, 1) = obj.z0;
        u_cbf = zeros(3, sim_steps);

        for k = 1:sim_steps 
            store = zeros(6,stages + 1); 
            sum_b = zeros(6,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_z_dyn_cbf(z_cbf(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            z_cbf(:, k+1) = z_cbf(:, k) + ts*sum_b; 
            u_cbf(:, k+1) = obj.safe_ctrl_fmincon(z_cbf(4:6, k)); 
    
        end

        %% Plot

        t_arr = 0:ts:T;
        
        x = z(1:3, :);
        x_cbf = z_cbf(1:3, :);
        v = z(4:6, :);
        v_cbf = z_cbf(4:6, :);
        
        color_v = [0 0.4470 0.7410];
        color_u = [0.8500 0.3250 0.0980];
        color_r = [0.9290 0.6940 0.1250];
        color_cbf_traj = [0 0.7290 0.5000]; 

        subplot(1,3,1);
        plot(t_arr,v(1, :),'Color',color_v)  
        hold on 
        plot(t_arr,v(2, :),'Color',color_u)
        plot(t_arr,v(3, :),'Color',color_r)
        plot(t_arr,v_cbf(1, :),'--','Color',color_v)  % stippled
        plot(t_arr,v_cbf(2, :),'--','Color',color_u)  % stippled
        plot(t_arr,v_cbf(3, :),'--','Color',color_r)  % stippled
        hold off 
        xlabel('Time')
        ylabel('State Variables nu = [v, u, r]')
        legend('v', 'u', 'r', 'v_cbf', 'u_cbf', 'r_cbf')
        
        subplot(1,3,2);
        plot(t_arr,u(1, :),'Color',color_v)  
        hold on
        plot(t_arr,u(2, :),'Color',color_u)
        plot(t_arr,u(3, :),'Color',color_r)
        plot(t_arr,u_cbf(1, :),'--','Color',color_v)  % stippled
        plot(t_arr,u_cbf(2, :),'--','Color',color_u)  % stippled
        plot(t_arr,u_cbf(3, :),'--','Color',color_r)  % stippled
        hold off 
        xlabel('Time')
        ylabel('Ctrl input')
        legend('u_1', 'u_2', 'u_3', 'u_cbf_1', 'u_cbf_2', 'u_cbf_3')
        
        % PLotting trajectories

        subplot(1,3,3);
        plot(x(1, :), x(2, :),'Color',color_v)
        hold on
        plot(x_cbf(1, :), x_cbf(2, :),'--','Color',color_cbf_traj)  % stippled
        
        % compute arrow lengths and directions
        arrow_frac = 0.2;  % fraction of distance between adjacent points to use for arrow length
        d = diff(x(1:2,:)).^2;
        dist = sqrt(sum(d,1));
        dist_frac = arrow_frac * dist;
        angles = x(3,:);
        angles_cbf = x_cbf(3, :);
        
        % compute arrow displacements
        delta_x = dist_frac .* cos(angles);
        delta_y = dist_frac .* sin(angles);
        delta_x_cbf = dist_frac .* cos(angles_cbf);
        delta_y_cbf = dist_frac .* sin(angles_cbf);
        
        % draw arrows
        quiver(x(1,1:end-1), x(2,1:end-1), delta_x(1:end-1), delta_y(1:end-1), 'Color', color_v, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        quiver(x_cbf(1,1:end-1), x_cbf(2,1:end-1), delta_x_cbf(1:end-1), delta_y_cbf(1:end-1), 'Color', color_cbf_traj, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        
        hold off
        xlabel('x')
        ylabel('y')
        legend('x', 'x_cbf', 'Orientation', 'horizontal')
        
        hold on
        scatter(x(1,1), x(2,1), 'x', 'LineWidth', 2, 'Color', color_v, 'DisplayName', 'x_i')
        scatter(x(1,end), x(2,end), 'o', 'LineWidth', 2, 'Color', color_v, 'DisplayName', 'x_f')
        scatter(x_cbf(1,end), x_cbf(2,end), 'o', 'LineWidth', 2, 'Color', color_v, 'DisplayName', 'x_f_cbf')
        hold off



        sgtitle('Tilte')

       
      end
   end
end