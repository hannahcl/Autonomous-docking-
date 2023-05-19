classdef ShipControll

   properties

        Q_lqr
        R_lqr
        K
        Kr

        cbf
        dyn
 
   end
   methods(Access = public)
      function obj = ShipControll()

         %Dynimacs of the ship
         obj.dyn = ShipDynamics(); 

         %Controll barrier functions
         obj.cbf = cbf(); 
        
 
         %Nominal controler with feedforward
         obj.Q_lqr = 10*eye(3); 
         obj.R_lqr= eye(3); 
         obj.K = lqr(obj.dyn.A_nu, obj.dyn.B_nu, obj.Q_lqr, obj.R_lqr); 
         obj.Kr = obj.dyn.B_nu*(inv(obj.dyn.C_nu/(obj.dyn.B_nu*obj.K-obj.dyn.A_nu)*obj.dyn.B_nu)); 

      end

      %% Controllers

      function nu_ref = ctrl_eta(obj, eta, eta_ref)
        nu_ref = -(eta - eta_ref); 
      end

      function tau = ctrl_nu_nominell(obj, nu, nu_ref)
        tau = obj.Kr*nu_ref-obj.K*nu; 
      end

      function [c, ceq] = compute_constraints(obj,tau, nu, eta, stage)
        z = [eta; nu]; 
        c = zeros(4, 1); 
        
        %Barrier for keeping the system from crashing into the dock
        if (stage == 0)
            c(1) = -(obj.cbf.Lf2_h1_fh(z) + obj.cbf.LgLf_h1_fh(z)*tau + obj.cbf.K1_alpha*[obj.cbf.h1_fh(z); obj.cbf.Lf_h1_fh(z)]); 
            c(2) = 0; 
            c(3) = 0; 
        else
            c(1) = -(obj.cbf.Lf2_h4_fh(z) + obj.cbf.LgLf_h4_fh(z)*tau + obj.cbf.K4_alpha*[obj.cbf.h4_fh(z); obj.cbf.Lf_h4_fh(z)]); 
            c(2) = -(obj.cbf.Lf2_h2_fh(z) + obj.cbf.LgLf_h2_fh(z)*tau + obj.cbf.K2_alpha*[obj.cbf.h2_fh(z); obj.cbf.Lf_h2_fh(z)]);
            c(3) = -(obj.cbf.Lf2_h3_fh(z) + obj.cbf.LgLf_h3_fh(z)*tau + obj.cbf.K3_alpha*[obj.cbf.h3_fh(z); obj.cbf.Lf_h3_fh(z)]);
        end

        %Barrier for keeping the system observable
        if (obj.cbf.ho1_fh(z) < obj.cbf.ho2_fh(z))
            c(4) = -(obj.cbf.Lf2_ho1_fh(z) + obj.cbf.LgLf_ho1_fh(z)*tau + obj.cbf.Ko1_alpha*[obj.cbf.ho1_fh(z); obj.cbf.Lf_ho1_fh(z)]);
        else
            c(4) = -(obj.cbf.Lf2_ho2_fh(z) + obj.cbf.LgLf_ho2_fh(z)*tau + obj.cbf.Ko2_alpha*[obj.cbf.ho2_fh(z); obj.cbf.Lf_ho2_fh(z)]);
        end
         ceq = 0;   
      end

      function tau_safe = ctrl_nu_safe(obj, tau_nominell, nu, eta, stage)

        f = @(tau_safe) norm(tau_safe - tau_nominell)^2;
        nonlcon = @(tau_safe) obj.compute_constraints(tau_safe, nu, eta, stage);

        options = optimoptions(@fmincon,'Display','none');
        tau_safe = fmincon(f, tau_nominell, [], [], [], [], [], [], nonlcon, options); 

      end

      function eta_ref = give_eta_ref(obj, stage)
          if (stage == 0)
            eta_ref = [0; -3.2; 0]; 
          else
            eta_ref = [0; 0; 0]; 
          end
      end

      %% Closed loop models

      function z_dot = closed_loop_model_z_nomiell(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 

          if (abs(eta(1)) > 1)
            stage = 0;  %Converge towards dock
          else
            stage = 1;  %Drive into berth
          end

          eta_ref = obj.give_eta_ref(stage); 
          nu_ref = obj.ctrl_eta(eta, eta_ref); 
          tau = obj.ctrl_nu_nominell(nu, nu_ref); 

          nu_dot = obj.dyn.model_nu(nu, tau); 
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
      end

      function z_dot = closed_loop_model_z_cbf(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 

          if (abs(eta(1)) > 1)
            stage = 0;  %Converge towards dock
          else
            stage = 1;  %Drive into berth
          end

          eta_ref = obj.give_eta_ref(stage); 
          nu_ref = obj.ctrl_eta(eta, eta_ref); 
          tau_nominell = obj.ctrl_nu_nominell(nu, nu_ref); 
          tau_safe = obj.ctrl_nu_safe(tau_nominell, nu, eta, stage); 

          nu_dot = obj.dyn.model_nu(nu, tau_safe); 
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
      end

   end

  methods

      function simRK4(obj, z0, T, ts)

        %Butcher tableau
        a = [0, 0.5, 0.5, 1]; 
        b = [1/6, 1/3, 1/3, 1/6]; 

        stages = size(b); stages = stages(2); 
        sim_steps = floor(T/ts); 
    
        %% Run simulation without cbf
        z = zeros(6, sim_steps);  
        z(:, 1) = z0; 

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
        z_cbf(:, 1) = z0;
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
        t_arr = 0:ts:T; t_arr = t_arr(1:sim_steps); 
        
        eta = z(1:3, :);
        eta_cbf = z_cbf(1:3, :);

        color_nom = [0 0.4470 0.7410];
        color_cbf = [0 0.7290 0.5000]; 

        % --- plot control barrier functions 
        subplot(1,2,1); 

        plot(t_arr, h1)
        hold on
        plot(t_arr, h2)
        hold on
        plot(t_arr, h3)
        hold on
        plot(t_arr, h4)
        hold on
        plot(t_arr, ho1)
        hold on
        plot(t_arr, ho2)
        hold off
        xlabel('Time')
        ylabel('')
        legend('h1', 'h2', 'h3', 'h4', 'ho1', 'ho2')

        % --- Plotting trajectory of eta ---
        subplot(1,2,2);
        plot(eta(1, :), eta(2, :),'Color',color_nom)
        hold on
        plot(eta_cbf(1, :), eta_cbf(2, :),'Color',color_cbf)  
        hold off

        % Draw vectors showing heading
        arrow_length = 0.2;  
        psi = eta(3,1:end-1);
        psi_safe = eta_cbf(3, 1:end-1);

        hold on
        quiver(eta(1,1:end-1), eta(2,1:end-1), arrow_length*cos(psi+pi/2), arrow_length*sin(psi+pi/2), 'Color', color_nom, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        quiver(eta_cbf(1,1:end-1), eta_cbf(2,1:end-1), arrow_length*cos(psi_safe+pi/2), arrow_length*sin(psi_safe+pi/2), 'Color', color_cbf, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        hold off
                
        % Draw starting and end points
        hold on
        scatter(eta(1,1), eta(2,1), 'x', 'LineWidth', 2, 'DisplayName', 'eta start')
        scatter(eta(1,end), eta(2,end), 'o', 'LineWidth', 2, 'DisplayName', 'eta finish')
        scatter(eta_cbf(1,end), eta_cbf(2,end), 'o', 'LineWidth', 2, 'DisplayName', 'eta_cbf finish')
        hold off
        
        % Draw contour of the dock
        l1x = [-5, -1]; l1y = [-3; -3]; 
        l2x = [-1, -1]; l2y = [-3; 0]; 
        l3x = [-1, 1]; l3y = [0; 0]; 
        l4x = [1, 1]; l4y = [0; -3]; 
        l5x = [1, 5]; l5y = [-3; -3]; 
        
        hold on
        plot(l1x,l1y, 'Color', [1 0.5 0]); 
        plot(l2x,l2y, 'Color', [1 0.5 0]);
        plot(l3x,l3y, 'Color', [1 0.5 0]);
        plot(l4x,l4y, 'Color', [1 0.5 0]);
        plot(l5x,l5y, 'Color', [1 0.5 0]);
        hold off

        legend('eta', 'eta safe')
        xlabel('x')
        ylabel('y')

        sgtitle('')

       
      end
   end
end

