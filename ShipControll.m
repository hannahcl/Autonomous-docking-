classdef ShipControll

   properties

        Q_lqr
        R_lqr
        K
        Kr

        cbf_ds0_1 %cbf for avioding dock in stage 0
        cbf_ds1_1 %cbf for avioding dock in stage 1
        cbf_ds1_2 %cbf for avioding dock in stage 1
        cbf_ds1_3 %cbf for avioding dock in stage 1
        cbf_o_1   %cbf for keeping the sys observable
        cbf_o_2   %cbf for keeping the sys observable

        cbfs
        dyn

        eta_measurement_variance
        nu_measurement_variance

 
   end
   methods(Access = public)
      function obj = ShipControll()

         %Dynimacs of the ship
         obj.dyn = ShipDynamics(); 

         %Controll barrier function handels
         % z' = f(z) + g(z)*tau, where z = [eta; nu]
         z = sym('z', [6 1]);

         h(z) = [0 -1 0 0 0 0]*z -3; 
         obj.cbf_ds0_1 = cbf(obj.dyn.f_symbolic, obj.dyn.g_symbolic, h, z); 

         h(z) = [0 -1 0 0 0 0]*z;
         obj.cbf_ds1_1 = cbf(obj.dyn.f_symbolic, obj.dyn.g_symbolic, h, z); 

         h(z) = [1 0 0 0 0 0]*z + 1;
         obj.cbf_ds1_2 = cbf(obj.dyn.f_symbolic, obj.dyn.g_symbolic, h, z);

         h(z) = [-1 0 0 0 0 0]*z +1;
         obj.cbf_ds1_3 = cbf(obj.dyn.f_symbolic, obj.dyn.g_symbolic, h, z);

         h(z) = -z(1) +z(2)*atan(z(3) - pi/4);
         obj.cbf_o_1 = cbf(obj.dyn.f_symbolic, obj.dyn.g_symbolic, h, z); 

         h(z) = z(1) +z(2)*atan(-z(3) - pi/4);
         obj.cbf_o_2 = cbf(obj.dyn.f_symbolic, obj.dyn.g_symbolic, h, z);

         obj.cbfs = {obj.cbf_ds0_1; obj.cbf_ds1_1; obj.cbf_ds1_2; obj.cbf_ds1_3; obj.cbf_o_1; obj.cbf_o_2}; 
 
         %Nominal controler with feedforward
         obj.Q_lqr = 10*eye(3); 
         obj.R_lqr= eye(3); 
         obj.K = lqr(obj.dyn.A_nu, obj.dyn.B_nu, obj.Q_lqr, obj.R_lqr); 
         obj.Kr = obj.dyn.B_nu*(inv(obj.dyn.C_nu/(obj.dyn.B_nu*obj.K-obj.dyn.A_nu)*obj.dyn.B_nu)); 

         obj.eta_measurement_variance = zeros(3,1); 
         obj.nu_measurement_variance = zeros(3,1);

      end

      %% Controllers

      function nu_ref = ctrl_eta(obj, eta, eta_ref)
        nu_ref = -(eta - eta_ref); 
      end

      function tau = ctrl_nu_nominell(obj, nu, nu_ref)
        tau = obj.Kr*nu_ref-obj.K*nu; 
      end

      function tau_safe = ctrl_nu_safe(obj, tau_nominell, nu, eta, stage)
        f = @(tau_safe) norm(tau_safe - tau_nominell)^2;
        nonlcon = @(tau_safe) obj.compute_constraints(tau_safe, nu, eta, stage);

        options = optimoptions(@fmincon,'Display','none');
        tau_safe = fmincon(f, tau_nominell, [], [], [], [], [], [], nonlcon, options); 
      end

      function [c, ceq] = compute_constraints(obj,tau, nu, eta, stage)
        z = [eta; nu]; 
        c = zeros(4, 1); 
        
        %Barrier for keeping the system from crashing into the dock
        if (stage == 0)
            c(1) = -obj.cbf_ds0_1.compute_bound(z, tau); 
            c(2) = 0; 
            c(3) = 0; 
        else
            c(1) = -obj.cbf_ds1_1.compute_bound(z, tau);
            c(2) = -obj.cbf_ds1_2.compute_bound(z, tau);
            c(3) = -obj.cbf_ds1_3.compute_bound(z, tau);
        end

        %Barrier for keeping the system observable
        if (obj.cbf_o_1.h(z) < obj.cbf_o_2.h(z))
            c(4) =  -obj.cbf_o_1.compute_bound(z, tau);
        else
            c(4) = -obj.cbf_o_2.compute_bound(z, tau);
        end
         ceq = 0;   
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
          eta_m = eta + randn(1)*obj.eta_measurement_variance; 
          nu_m = nu + randn(1)*obj.nu_measurement_variance;  

          if (abs(eta_m(1)) > 1)
            stage = 0;  %Converge towards dock
          else
            stage = 1;  %Drive into berth
          end

          eta_ref = obj.give_eta_ref(stage); 
          nu_ref = obj.ctrl_eta(eta_m, eta_ref); 
          tau = obj.ctrl_nu_nominell(nu_m, nu_ref); 

          nu_dot = obj.dyn.model_nu(nu, tau);
          %nu_dot = obj.dyn.linear_model_nu(nu, tau); 
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
      end

      function z_dot = closed_loop_model_z_cbf(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 
          eta_m = eta + randn(1)*obj.eta_measurement_variance; 
          nu_m = nu + randn(1)*obj.nu_measurement_variance; 

          if (abs(eta_m(1)) > 1)
            stage = 0;  %Converge towards dock
          else
            stage = 1;  %Drive into berth
          end

          eta_ref = obj.give_eta_ref(stage); 
          nu_ref = obj.ctrl_eta(eta_m, eta_ref); 
          tau_nominell = obj.ctrl_nu_nominell(nu_m, nu_ref); 
          tau_safe = obj.ctrl_nu_safe(tau_nominell, nu_m, eta_m, stage); 

          nu_dot = obj.dyn.model_nu(nu, tau_safe); 
          %nu_dot = obj.dyn.linear_model_nu(nu, tau_safe);
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
      end

   end

  methods
      function measurement = add_noise(obj, z)
        measurement = z + randn(1)*[obj.eta_measurement_variance; obj.nu_measurement_variance];
      end

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
        for i = 1:numel(obj.cbfs)
            h{i} = zeros(1, sim_steps);
        end


        for k = 1:sim_steps 
            store = zeros(6,stages + 1); 
            sum_b = zeros(6,1); 
            for s = 1:stages
                store(:, s+1) = obj.closed_loop_model_z_cbf(z_cbf(:, k) + ts*a(s)*store(:, s));  
                sum_b  = sum_b + b(s)*store(:,s+1); 
            end
    
            z_cbf(:, k+1) = z_cbf(:, k) + ts*sum_b; 

            for i = 1:numel(obj.cbfs)
                h{i}(k) = obj.cbfs{i}.h(z_cbf(:, k));
            end

        end

        %% Plot
        t_arr = 0:ts:T; t_arr = t_arr(1:sim_steps); 
        
        eta = z(1:3, :);
        eta_cbf = z_cbf(1:3, :);

        color_nom = [0 0.4470 0.7410];
        color_cbf = [0 0.7290 0.5000]; 

        % --- plot control barrier functions 
        subplot(1,2,1); 

        for i = 1:numel(obj.cbfs)
            plot(t_arr, h{i});
            hold on
        end
        hold off
        xlabel('Time')
        ylabel('')
        legend('cbf ds0 1', 'cbf ds1 1','cbf ds1 2', 'cbf ds1 3','cbf o 1', 'cbf o 2')

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
        
        step = 10;
        indices = 1:step:length(eta)-1;
        
        quiver(eta(1,indices), eta(2,indices), arrow_length*cos(psi(indices)+pi/2), arrow_length*sin(psi(indices)+pi/2), 'Color', color_nom, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        quiver(eta_cbf(1,indices), eta_cbf(2,indices), arrow_length*cos(psi_safe(indices)+pi/2), arrow_length*sin(psi_safe(indices)+pi/2), 'Color', color_cbf, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        
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

