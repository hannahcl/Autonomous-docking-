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

        cbf
 
        nu0 
        eta0
        z0

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

         %Linear model with normal notation for nu dynamics
         obj.A = -inv(obj.M)*obj.N;
         obj.B = eye(3); 
         obj.C = eye(3);

         %Nominal controler with feedforward
         obj.Q_lqr = 10*eye(3); 
         obj.R_lqr= eye(3); 
         obj.K = lqr(obj.A, obj.B, obj.Q_lqr, obj.R_lqr); 
         obj.Kr = obj.B*(inv(obj.C/(obj.B*obj.K-obj.A)*obj.B)); 

         %CBF controller
%          obj.alpha = 1;   
%          obj.H = [0 1 0]; 
%          obj.k_h = 0; 

         %paramters
         obj.nu0 = [3; 1; pi/4];
         obj.eta0 = [5; -3.5; 1]; 
         obj.z0 = [obj.eta0; obj.nu0];

         obj.cbf = cbf(); 
         
      end

      %% Controllers

      function nu_ref = nominell_ctrl_nu_ref(obj, eta)
        nu_ref = -eta; 
      end

      function [c, ceq] = check_barrier_func(obj, nu, eta)

        c = zeros(4, 1); 
        
        for i=1:4
            c(i) = -(obj.cbf.grad_h1_fh(eta)*obj.g(eta)*nu + obj.cbf.alpha*obj.cbf.h1_fh(eta)); 

        end

         ceq = 0;   
      end


      function nu_ref_safe = safe_ctrl_nu_ref(obj, eta)
        nu_ref = obj.nominell_ctrl_nu_ref(eta); 

        f = @(nu_ref_safe) norm(nu_ref_safe - nu_ref)^2;
        nonlcon = @(nu_ref_safe) obj.check_barrier_func(nu_ref_safe, eta);

        options = optimoptions(@fmincon,'Display','none');
        nu_ref_safe = fmincon(f, nu_ref, [], [], [], [], [], [], nonlcon, options);

      end


      function tau = ctrl_nu(obj, nu, nu_ref)
        tau = obj.Kr*nu_ref-obj.K*nu; 
      end

      %% Models

      function g_val = g(obj, eta)
        g_val = [cos(eta(3)) -sin(eta(3)) 0; 
             sin(eta(3)) cos(eta(3)) 0; 
             0 0 1] ;     
      end

      function nu_dot = model_nu_dyn(obj, nu, u)
        nu_dot = obj.A*nu + obj.B*u;
      end

      function nu_dot = closed_loop_model_nu_dyn(obj, nu, nu_ref)
          u = obj.ctrl_nu(nu, nu_ref); 
          nu_dot = obj.model_nu_dyn(nu, u); 
      end

      function eta_dot = model_eta_dyn(obj, eta, nu)
         R = [cos(eta(3)) -sin(eta(3)) 0; 
             sin(eta(3)) cos(eta(3)) 0; 
             0 0 1]; 
         eta_dot = R*nu; 
      end

      function z_dot = closed_loop_model_z_dyn_nomiell(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 

          nu_ref = obj.nominell_ctrl_nu_ref(eta); 
          nu_dot = obj.closed_loop_model_nu_dyn(nu, nu_ref); 
          eta_dot = obj.model_eta_dyn(eta, nu);

          z_dot = [eta_dot; nu_dot]; 
      end

      function z_dot = closed_loop_model_z_dyn_cbf(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 

          disp('h(eta)')
          obj.cbf.h1_fh(eta)
          
          disp('lie bound')
          obj.cbf.grad_h1_fh(eta)*obj.g(eta)*nu + obj.cbf.alpha*obj.cbf.h1_fh(eta)

          nu_ref_safe = obj.safe_ctrl_nu_ref(eta); 
          nu_dot = obj.closed_loop_model_nu_dyn(nu, nu_ref_safe); 
          eta_dot = obj.model_eta_dyn(eta, nu); 

          z_dot = [eta_dot; nu_dot]; 
      end



   end

  methods

      function sim(obj)
          
        T = 5; 
        ts = 0.2; 
     
        a = [0, 0.5, 0.5, 1]; 
        b = [1/6, 1/3, 1/3, 1/6]; 

        stages = size(b); stages = stages(2); 
        sim_steps = floor(T/ts); 
    
        %% Run simulation without cbf
        z = zeros(6, sim_steps);  
        z(:, 1) = obj.z0; 
        nu_ref = zeros(3, sim_steps); 
    
        for k = 1:sim_steps 
            store = zeros(6,stages + 1); 
            sum_b = zeros(6,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_z_dyn_nomiell(z(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            z(:, k+1) = z(:, k) + ts*sum_b; 
            nu_ref(:, k+1) = obj.nominell_ctrl_nu_ref(z(1:3, k)); 
    
        end

        %% Run simulation with cbf
        z_cbf = zeros(6, sim_steps);
        z_cbf(:, 1) = obj.z0;
        nu_ref_safe = zeros(3, sim_steps);

        for k = 1:sim_steps 
            store = zeros(6,stages + 1); 
            sum_b = zeros(6,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_z_dyn_cbf(z_cbf(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            z_cbf(:, k+1) = z_cbf(:, k) + ts*sum_b; 
            nu_ref_safe(:, k+1) = obj.safe_ctrl_nu_ref(z(1:3, k)); 
    
        end

        %% Plot

        t_arr = 0:ts:T;
        
        eta = z(1:3, :);
        eta_cbf = z_cbf(1:3, :);
        nu = z(4:6, :);
        nu_cbf = z_cbf(4:6, :);
        
        color_v = [0 0.4470 0.7410];
        color_u = [0.8500 0.3250 0.0980];
        color_r = [0.9290 0.6940 0.1250];
        color_cbf_traj = [0 0.7290 0.5000]; 

        subplot(1,3,1);
        plot(t_arr,nu_ref(1, :),'Color',color_v)  
        hold on
        plot(t_arr,nu_ref(2, :),'Color',color_u)
        plot(t_arr,nu_ref(3, :),'Color',color_r)
        plot(t_arr,nu_ref_safe(1, :),'--','Color',color_v)  % stippled
        plot(t_arr,nu_ref_safe(2, :),'--','Color',color_u)  % stippled
        plot(t_arr,nu_ref_safe(3, :),'--','Color',color_r)  % stippled
        hold off 
        xlabel('Time')
        ylabel('Ctrl input')
        legend('u ref', 'v ref', 'r ref', 'u ref safe', 'v ref safe', 'r ref safe')
        
        
        subplot(1,3,2);
        plot(t_arr,nu(1, :),'Color',color_v)  
        hold on 
        plot(t_arr,nu(2, :),'Color',color_u)
        plot(t_arr,nu(3, :),'Color',color_r)
        plot(t_arr,nu_cbf(1, :),'--','Color',color_v)  % stippled
        plot(t_arr,nu_cbf(2, :),'--','Color',color_u)  % stippled
        plot(t_arr,nu_cbf(3, :),'--','Color',color_r)  % stippled
        hold off 
        xlabel('Time')
        ylabel('State Variables nu = [v, u, r]')
        legend('v', 'u', 'r', 'v_cbf', 'u_cbf', 'r_cbf')

        % --- Plotting trajectory of eta ---

        % Plot x vs y
        subplot(1,3,3);
        plot(eta(1, :), eta(2, :),'Color',color_v)
        hold on
        plot(eta_cbf(1, :), eta_cbf(2, :),'--','Color',color_cbf_traj)  % stippled
        hold off

                
        % Draw vectors showwing heading
        
        arrow_length = 0.2;  
        psi = eta(3,1:end-1);
        psi_safe = eta_cbf(3, 1:end-1);

        hold on
        quiver(eta(1,1:end-1), eta(2,1:end-1), arrow_length*cos(psi), arrow_length*sin(psi), 'Color', color_v, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        quiver(eta_cbf(1,1:end-1), eta_cbf(2,1:end-1), arrow_length*cos(psi_safe), arrow_length*sin(psi_safe), 'Color', color_cbf_traj, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        hold off
        xlabel('x')
        ylabel('y')
        legend('eta', 'eta_safe')
        
        
        % Draw starting and end points
        hold on
        scatter(eta(1,1), eta(2,1), 'x', 'LineWidth', 2, 'Color', color_v, 'DisplayName', 'eta_i')
        scatter(eta(1,end), eta(2,end), 'o', 'LineWidth', 2, 'Color', color_v, 'DisplayName', 'eta_f')
        scatter(eta_cbf(1,end), eta_cbf(2,end), 'o', 'LineWidth', 2, 'Color', color_v, 'DisplayName', 'eta_f_cbf')
        hold off

        % Draw contour of doc
        
        x_plot= -5:0.1:5;
        f = zeros(1, length(x_plot)); 
        
        for i =1:length(x_plot)
           
            f(i, :) = obj.cbf.f(x_plot(i)); 
        end
        hold on
        plot(x_plot, f);
        hold off

        legend('eta', 'eta_safe', 'eta_i', 'eta_f', 'eta_f_cbf', 'doc')

        sgtitle('Tilte')

       
      end
   end
end

