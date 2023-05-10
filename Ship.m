classdef Ship

   properties

        B
        C
        A
        Q_lqr
        R_lqr
        K
        Kr

        cbf
        dyn
 
        nu0 
        eta0
        z0

   end
   methods(Access = public)
      function obj = Ship()

         obj.cbf = cbf(); 
         obj.dyn = ShipDynamics(); 

         %Nominal controler with feedforward
         obj.A = -obj.dyn.M\obj.dyn.N_lin;
         obj.B = eye(3); 
         obj.C = eye(3); 
         obj.Q_lqr = 10*eye(3); 
         obj.R_lqr= eye(3); 
         obj.K = lqr(obj.A, obj.B, obj.Q_lqr, obj.R_lqr); 
         obj.Kr = obj.B*(inv(obj.C/(obj.B*obj.K-obj.A)*obj.B)); 

         %initial values 
         obj.nu0 = [-1; 1; 0];
         obj.eta0 = [-1.2; -4.2; 0]; 
         obj.z0 = [obj.eta0; obj.nu0];

         
         
      end

      %% Controllers

      function nu_ref = ctrl_eta(obj, eta)
        nu_ref = -eta; 
      end

      function tau = ctrl_nu_nominell(obj, nu, nu_ref)
        tau = obj.Kr*nu_ref-obj.K*nu; 
      end

      function [c, ceq] = compute_constraints(obj,tau, nu, eta)
        z = [eta; nu]; 
        c = zeros(4, 1); 
        for i=1:4
           % c(i) = -(obj.cbf.grad_h1_fh(eta)*obj.dyn.compute_R(eta)*nu + obj.cbf.alpha*obj.cbf.h1_fh(eta)); 
            c(i) = -(obj.cbf.Lf2_h_fh(z) + obj.cbf.LfLg_h_fh(z)*tau + obj.cbf.K_alpha*[obj.cbf.h_fh(z); obj.cbf.Lf_h_fh(z)]); 
        end
         ceq = 0;   
      end

      function tau_safe = ctrl_nu_safe(obj, tau_nominell, nu, eta)

        f = @(tau_safe) norm(tau_safe - tau_nominell)^2;
        nonlcon = @(tau_safe) obj.compute_constraints(tau_safe, nu, eta);

        options = optimoptions(@fmincon,'Display','none');
        tau_safe = fmincon(f, tau_nominell, [], [], [], [], [], [], nonlcon, options);

      end

      %% Closed loop models

      function z_dot = closed_loop_model_z_nomiell(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 

          nu_ref = obj.ctrl_eta(eta); 
          tau = obj.ctrl_nu_nominell(nu, nu_ref); 

          nu_dot = obj.dyn.linear_model_nu(nu, tau); %OBS change to nonlienar model
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
      end

      function z_dot = closed_loop_model_z_cbf(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 
          
          nu_ref = obj.ctrl_eta(eta); 
          tau_nominell = obj.ctrl_nu_nominell(nu, nu_ref); 
          tau_safe = obj.ctrl_nu_safe(tau_nominell, nu, eta); 

          nu_dot = obj.dyn.linear_model_nu(nu, tau_safe); %OBS change to nonlienar model
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
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
        z = zeros(6, sim_steps);  
        z(:, 1) = obj.z0; 
        %tau = zeros(3, sim_steps); 
    
        for k = 1:sim_steps 
            store = zeros(6,stages + 1); 
            sum_b = zeros(6,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_z_nomiell(z(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            z(:, k+1) = z(:, k) + ts*sum_b; 
            %tau(:, k+1) = obj.ctrl_nu_nominell(z(1:3, k)); 
    
        end

        %% Run simulation with cbf
        z_cbf = zeros(6, sim_steps);
        z_cbf(:, 1) = obj.z0;
        %tau_safe = zeros(3, sim_steps);

        for k = 1:sim_steps 
            store = zeros(6,stages + 1); 
            sum_b = zeros(6,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_z_cbf(z_cbf(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            z_cbf(:, k+1) = z_cbf(:, k) + ts*sum_b; 
            %tau_safe(:, k+1) = obj.ctrl_nu_safe(z(1:3, k)); 
    
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

        %subplot(1,2,1);
        % --- plot tau

%         plot(t_arr,tau(1, :),'Color',color_v)  
%         hold on
%         plot(t_arr,tau(2, :),'Color',color_u)
%         plot(t_arr,tau(3, :),'Color',color_r)
%         plot(t_arr,tau_safe(1, :),'--','Color',color_v)  % stippled
%         plot(t_arr,tau_safe(2, :),'--','Color',color_u)  % stippled
%         plot(t_arr,tau_safe(3, :),'--','Color',color_r)  % stippled
%         hold off 
%         xlabel('Time')
%         ylabel('Tau')
%         legend('tau 1', 'tau 2', 'tau 3', 'tau 1 safe', 'tau 2 safe', 'tau 3 safe')
        
        
        subplot(1,2,1);
        % --- plot nu 

%         plot(t_arr,nu(1, :),'Color',color_v)  
%         hold on 
%         plot(t_arr,nu(2, :),'Color',color_u)
%         plot(t_arr,nu(3, :),'Color',color_r)
%         plot(t_arr,nu_cbf(1, :),'--','Color',color_v)  % stippled
%         plot(t_arr,nu_cbf(2, :),'--','Color',color_u)  % stippled
%         plot(t_arr,nu_cbf(3, :),'--','Color',color_r)  % stippled
%         hold off 
%         xlabel('Time')
%         ylabel('Nu')
%         legend('v', 'u', 'r', 'v_cbf', 'u_cbf', 'r_cbf')


        % --- plot h and gradient of h

        h_arr = zeros(size(t_arr));
        deta_h_arr = zeros(size(t_arr));

        for i = 1:length(t_arr)
            %-(obj.cbf.Lf2_h_fh(z) + obj.cbf.LfLg_h_fh(z)*tau + obj.cbf.K_alpha*[obj.cbf.h_fh(z); obj.cbf.Lf_h_fh(z)]);
            h_arr(i) = obj.cbf.h_fh([eta_cbf(:, i); nu_cbf(:, i)]);
            deta_h_arr(i) = 0;
        end

        plot(t_arr, h_arr, '-o')
        hold on
        plot(t_arr, deta_h_arr, '-o')
        hold off
        xlabel('Time')
        ylabel('')
        legend('h', 'd_{eta} h')

        % --- Plotting trajectory of eta ---

        % Plot x vs y
        subplot(1,2,2);
        plot(eta(1, :), eta(2, :),'Color',color_v)
        hold on
        plot(eta_cbf(1, :), eta_cbf(2, :),'--','Color',color_cbf_traj)  % stippled
        hold off

                
        % Draw vectors showing heading
        
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

