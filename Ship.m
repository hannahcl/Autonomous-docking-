classdef Ship

   properties

        Q_lqr
        R_lqr
        K
        Kr

        cbf
        dyn
        nu_current
 
        nu0 
        eta0
        z0

   end
   methods(Access = public)
      function obj = Ship()

         obj.cbf = cbf(); 
         obj.dyn = ShipDynamics(); 
         obj.nu_current = zeros(3,1); 

         %Nominal controler with feedforward
         obj.Q_lqr = 10*eye(3); 
         obj.R_lqr= eye(3); 
         obj.K = lqr(obj.dyn.A_nu, obj.dyn.B_nu, obj.Q_lqr, obj.R_lqr); 
         obj.Kr = obj.dyn.B_nu*(inv(obj.dyn.C_nu/(obj.dyn.B_nu*obj.K-obj.dyn.A_nu)*obj.dyn.B_nu)); 

         %initial values 
         obj.nu0 = [0; 0; 0];
         obj.eta0 = [-0.0; -4; pi/8]; 
         obj.z0 = [obj.eta0; obj.nu0];
         obj.dyn.tau_wind = [0; 0; 0]; 

         
         
      end

      %% Controllers

      function nu_ref = ctrl_eta(obj, eta, eta_ref)
        nu_ref = -(eta - eta_ref); 
      end

      function tau = ctrl_nu_nominell(obj, nu, nu_ref)
        tau = obj.Kr*nu_ref-obj.K*nu; 
      end

      function [c, ceq] = compute_constraints(obj,tau, nu, eta)
        z = [eta; nu]; 
        c = zeros(5, 1); 
        
        if (abs(z(1)) > 1)
            %Barrier for stage 1, before driving into berth
            c(1) = -(obj.cbf.Lf2_h1_fh(z) + obj.cbf.LgLf_h1_fh(z)*tau + obj.cbf.K1_alpha*[obj.cbf.h1_fh(z); obj.cbf.Lf_h1_fh(z)]); 
            c(2) = 0; 
            c(3) = 0; 
        else
            %Barrier for stage 2, when driving into berth
            c(1) = -(obj.cbf.Lf2_h4_fh(z) + obj.cbf.LgLf_h4_fh(z)*tau + obj.cbf.K4_alpha*[obj.cbf.h4_fh(z); obj.cbf.Lf_h4_fh(z)]); 
            c(2) = -(obj.cbf.Lf2_h2_fh(z) + obj.cbf.LgLf_h2_fh(z)*tau + obj.cbf.K2_alpha*[obj.cbf.h2_fh(z); obj.cbf.Lf_h2_fh(z)]);
            c(3) = -(obj.cbf.Lf2_h3_fh(z) + obj.cbf.LgLf_h3_fh(z)*tau + obj.cbf.K3_alpha*[obj.cbf.h3_fh(z); obj.cbf.Lf_h3_fh(z)]);
        end

        %Barrier for keeping the system observable
%         if (obj.cbf.ho1_fh(z) < obj.cbf.ho2_fh(z))
%             c(4) = -(obj.cbf.Lf2_ho1_fh(z) + obj.cbf.LgLf_ho1_fh(z)*tau + obj.cbf.Ko1_alpha*[obj.cbf.ho1_fh(z); obj.cbf.Lf_ho1_fh(z)]);
%         else
%             c(4) = -(obj.cbf.Lf2_ho2_fh(z) + obj.cbf.LgLf_ho2_fh(z)*tau + obj.cbf.Ko2_alpha*[obj.cbf.ho2_fh(z); obj.cbf.Lf_ho2_fh(z)]);
%         end

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

          if (abs(eta(1)) > 1)
            eta_ref = [0; -3.2; 0]; 
          else
            eta_ref = [0; 0; 0]; 
          end

          nu_ref = obj.ctrl_eta(eta, eta_ref); 
          tau = obj.ctrl_nu_nominell(nu, nu_ref); 

          nu_dot = obj.dyn.model_nu(nu, obj.nu_current, tau); 
          %nu_dot = obj.dyn.linear_model_nu(nu, tau); 
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
      end

      function z_dot = closed_loop_model_z_cbf(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 


          if (abs(eta(1)) > 1)
            eta_ref = [0; -3.2; 0]; 
          else
            eta_ref = [0; 0; 0]; 
          end
          
          nu_ref = obj.ctrl_eta(eta, eta_ref); 
          tau_nominell = obj.ctrl_nu_nominell(nu, nu_ref); 
          tau_safe = obj.ctrl_nu_safe(tau_nominell, nu, eta); 

%           % ----
%           if (abs(obj.cbf.ho1_fh(z))> 15 )
%               z
%               tau_safe
%           end
%           % ----

          nu_dot = obj.dyn.model_nu(nu, nu_ref, tau_safe); %OBS change to nonlienar model
          %nu_dot = obj.dyn.linear_model_nu(nu, tau_safe);
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
      end

   end

  methods

      function sim(obj)
          
        T = 50; 
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

        end

        %% Run simulation with cbf
        z_cbf = zeros(6, sim_steps);
        z_cbf(:, 1) = obj.z0;
        h1 = zeros(1,sim_steps);
        h2 = zeros(1,sim_steps);
        h3 = zeros(1,sim_steps);
        h4 = zeros(1,sim_steps);
        ho1 = zeros(1,sim_steps);
        ho2 = zeros(1,sim_steps);


        for k = 1:sim_steps 
            store = zeros(6,stages + 1); 
            sum_b = zeros(6,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_z_cbf(z_cbf(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            z_cbf(:, k+1) = z_cbf(:, k) + ts*sum_b; 
            h1(k) = obj.cbf.h1_fh(z_cbf(:, k)); 
            h2(k) = obj.cbf.h2_fh(z_cbf(:, k));  
            h3(k) = obj.cbf.h3_fh(z_cbf(:, k)); 
            h4(k) = obj.cbf.h4_fh(z_cbf(:, k));
            ho1(k) = obj.cbf.ho1_fh(z_cbf(:, k));
            ho2(k) = obj.cbf.ho2_fh(z_cbf(:, k));

        end

        %% Plot

        t_arr = 0:ts:T;
        t_arr = t_arr(1:sim_steps); 
        
        eta = z(1:3, :);
        eta_cbf = z_cbf(1:3, :);

        color_v = [0 0.4470 0.7410];
        color_u = [0.8500 0.3250 0.0980];
        color_r = [0.9290 0.6940 0.1250];
        color_cbf_traj = [0 0.7290 0.5000]; 

        subplot(1,2,1);
        % --- plot h 

        plot(t_arr, h1, '-o')
        hold on
        plot(t_arr, h2, '-o')
        hold on
        plot(t_arr, h3, '-o')
        hold on
        plot(t_arr, h4, '-o')
        hold on
        plot(t_arr, ho1, '-o')
        hold on
        plot(t_arr, ho2, '-o')
        hold off
        xlabel('Time')
        ylabel('')
        legend('h1', 'h2', 'h3', 'h4', 'ho1', 'ho2')

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
        quiver(eta(1,1:end-1), eta(2,1:end-1), arrow_length*cos(psi+pi/2), arrow_length*sin(psi+pi/2), 'Color', color_v, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        quiver(eta_cbf(1,1:end-1), eta_cbf(2,1:end-1), arrow_length*cos(psi_safe+pi/2), arrow_length*sin(psi_safe+pi/2), 'Color', color_cbf_traj, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
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

        % Draw contour of doc: TODO: replace with lines
        
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

