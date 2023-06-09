classdef ShipControl < handle
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
        active_cbfs_stage0
        active_cbfs_stage1

        waypoint_stage0
        waypoint_stage1

        dyn

        z0
        tau0

        eta_measurement_variance
        nu_measurement_variance

        n_active_cbfs
        n_ctrl_inputs
        dock_width
        dock_length

        stage
        tentative_switch_stage

 
   end
   methods(Access = public)
      function obj = ShipControl(z0, tau0)

         %Initial conditions
         obj.z0 = z0; 
         obj.tau0 = tau0; 

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

         h(z) = -z(1) +20*z(2)*atan(0.1*(z(3) - pi/3));
         obj.cbf_o_1 = cbf(obj.dyn.f_symbolic, obj.dyn.g_symbolic, h, z); 

         h(z) = z(1) +20*z(2)*atan(0.1*(-z(3) - pi/3));
         obj.cbf_o_2 = cbf(obj.dyn.f_symbolic, obj.dyn.g_symbolic, h, z);

         obj.cbfs = {obj.cbf_ds0_1; obj.cbf_ds1_1; obj.cbf_ds1_2; obj.cbf_ds1_3; obj.cbf_o_1; obj.cbf_o_2}; 
         obj.active_cbfs_stage0 = {obj.cbf_ds0_1;obj.cbf_o_1; obj.cbf_o_2};
         obj.active_cbfs_stage1 = {obj.cbf_ds1_1; obj.cbf_ds1_2; obj.cbf_ds1_3; obj.cbf_o_1; obj.cbf_o_2}; 

         obj.n_active_cbfs = 4;
         obj.n_ctrl_inputs = 3; 
 
         %Nominal controler with feedforward
         obj.Q_lqr = 10*eye(3);   
         obj.R_lqr= eye(3); 
         obj.K = lqr(obj.dyn.A_nu, obj.dyn.B_nu, obj.Q_lqr, obj.R_lqr); 
         obj.Kr = obj.dyn.B_nu*(inv(obj.dyn.C_nu/(obj.dyn.B_nu*obj.K-obj.dyn.A_nu)*obj.dyn.B_nu)); 

         obj.eta_measurement_variance = zeros(3,1); 
         obj.nu_measurement_variance = zeros(3,1);

         obj.dock_width = 1; 
         obj.dock_length = 3; 

         obj.stage = 0; 
         obj.tentative_switch_stage = false; 

         obj.waypoint_stage0 = [0; -(obj.dock_length+0.2); 0];
         obj.waypoint_stage1 = [0; 0; 0]; 

      end

      %% Controllers

      function nu_ref = ctrl_eta(obj, eta, eta_ref)
        nu_ref = -(eta - eta_ref); 
      end

      function tau = ctrl_nu_nominell(obj, nu, nu_ref)
        tau = obj.Kr*nu_ref-obj.K*nu; 
      end

      function tau_safe = ctrl_nu_safe(obj, tau_nominell, nu, eta, stage)
        %Create matricies for tau_safe = quadprog(H,f,lin_con_A, lin_con_b).
        H = 2*eye(3); 
        f = -tau_nominell; 

        lin_con_A = zeros(obj.n_active_cbfs, obj.n_ctrl_inputs); 
        lin_con_b = zeros(obj.n_active_cbfs, 1); 

        z = [eta; nu]; 

        %Barrier for avioding collision
        if (stage == 0)
            lin_con_s01 = obj.cbf_ds0_1.compute_linear_constraints(z); 
            lin_con_A(1,:) = lin_con_s01(1:3); 
            lin_con_b(1) = lin_con_s01(4); 
        else
            lin_con_s11 = obj.cbf_ds1_1.compute_linear_constraints(z); 
            lin_con_A(1,:) = lin_con_s11(1:3); 
            lin_con_b(1) = lin_con_s11(4);

            lin_con_s12 = obj.cbf_ds1_2.compute_linear_constraints(z); 
            lin_con_A(2,:) = lin_con_s12(1:3); 
            lin_con_b(2) = lin_con_s12(4);

            lin_con_s13 = obj.cbf_ds1_3.compute_linear_constraints(z); 
            lin_con_A(3,:) = lin_con_s13(1:3); 
            lin_con_b(3) = lin_con_s13(4);
        end

        %Barrier for keeping the system observable
        if (obj.cbf_o_1.h(z) < obj.cbf_o_2.h(z))
            lin_con_o1 = obj.cbf_o_1.compute_linear_constraints(z); 
            lin_con_A(4,:) = lin_con_o1(1:3); 
            lin_con_b(4) = lin_con_o1(4);

        else
            lin_con_o2 = obj.cbf_o_2.compute_linear_constraints(z); 
            lin_con_A(4,:) = lin_con_o2(1:3); 
            lin_con_b(4) = lin_con_o2(4);
        end

        options = optimoptions('quadprog', 'Display', 'off');
        tau_safe = quadprog(H,f,lin_con_A, lin_con_b, [], [], [], [], tau_nominell, options); 

      end

      function eta_ref = give_eta_ref(obj, stage)
          if (stage == 0)
            eta_ref = obj.waypoint_stage0; 
          else
            eta_ref = obj.waypoint_stage1; 
          end 
      end

      function tau_safe = ctrl_z_safe(obj, z)
          eta_ref = obj.give_eta_ref(obj.stage); 
          nu_ref = obj.ctrl_eta(z(1:3), eta_ref); 
          tau_nominell = obj.ctrl_nu_nominell(z(4:6), nu_ref); 
          tau_safe = obj.ctrl_nu_safe(tau_nominell, z(4:6), z(1:3), obj.stage); 
      end

      function tau_nominell = ctrl_z_nominell(obj, z)
          eta_ref = obj.give_eta_ref(obj.stage); 
          nu_ref = obj.ctrl_eta(z(1:3), eta_ref); 
          tau_nominell = obj.ctrl_nu_nominell(z(4:6), nu_ref); 
      end

      function s = update_stage(obj, eta)
          if (abs(eta(1) + obj.eta_measurement_variance(1,1)) > obj.dock_width)
            s = 0;  %Converge towards dock
          else
            s = 1;  %Drive into berth
          end

          if ~(s == obj.stage)
              obj.tentative_switch_stage = true;
              disp('tentative switch')
          end
          obj.stage = s; 
      end

      %% Closed loop models

      function z_dot = closed_loop_model_z_nomiell(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 
          z_measured = [eta + randn(1)*obj.eta_measurement_variance;  nu + randn(1)*obj.nu_measurement_variance]; 

          obj.update_stage(z_measured(1:3)); 

          tau = obj.ctrl_z_nominell(z_measured); 

          nu_dot = obj.dyn.model_nu(nu, tau);
          %nu_dot = obj.dyn.linear_model_nu(nu, tau); 
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot]; 
      end

      function z_dot = closed_loop_model_z_cbf(obj, z)
          eta = z(1:3); 
          nu = z(4:6); 
          z_measured = [eta + randn(1)*obj.eta_measurement_variance;  nu + randn(1)*obj.nu_measurement_variance]; 

          obj.update_stage(z_measured(1:3));
          tau_safe = obj.ctrl_z_safe(z_measured); 

          if obj.tentative_switch_stage
              obj.tentative_switch_stage = false; 
              if (obj.stage == 1)
                for i = 1:numel(obj.active_cbfs_stage1)
                    cbf_valid = obj.active_cbfs_stage1{i}.check_initial_conditions(z, tau_safe);
                    if ~cbf_valid
                        obj.stage = 0; 
                        disp('Could not switch stage.')
                        tau_safe = obj.ctrl_z_safe(z_measured);
                    end
                    break
                end
              end   
          end

          nu_dot = obj.dyn.model_nu(nu, tau_safe);
          %nu_dot = obj.dyn.linear_model_nu(nu, tau_safe); 
          eta_dot = obj.dyn.model_eta(eta, nu);
          z_dot = [eta_dot; nu_dot];
      end

   end

  methods


      function simRK4(obj, T, ts)

        %Butcher tableau
        a = [0, 0.5, 0.5, 1]; 
        b = [1/6, 1/3, 1/3, 1/6]; 

        stages = size(b); stages = stages(2); 
        sim_steps = floor(T/ts); 
    
        %% Run simulation without cbf
        z = zeros(6, sim_steps);  
        z(:, 1) = obj.z0; 
        obj.stage = 0;
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
        for i = 1:numel(obj.cbfs)
            h{i} = zeros(1, sim_steps);
        end
        stage_arr = zeros(1, sim_steps);

        obj.stage = 0;

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
            stage_arr(k) = obj.update_stage(z_cbf(1:3, k)); 

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
        plot(t_arr, stage_arr); 
        hold off 

        xlabel('Time')
        ylabel('')
        legend('cbf ds0 1', 'cbf ds1 1','cbf ds1 2', 'cbf ds1 3','cbf o 1', 'cbf o 2', 'stage')

        % --- Plotting trajectory of eta ---
        subplot(1,2,2);
        plot(eta(1, :), eta(2, :),'Color',color_nom)
        hold on
        plot(eta_cbf(1, :), eta_cbf(2, :),'Color',color_cbf)  
        hold off

        % Draw vectors showing heading
        arrow_length = 0.2;  
        psi = eta(3,1:end-1); 
        psi_safe = eta_cbf(3, 1:end-1) ;
        
        hold on
        
        step = 10;
        indices = 1:step:length(eta)-1;
        
        quiver(eta(1,indices), eta(2,indices), arrow_length*cos(-psi(indices)+pi/2), arrow_length*sin(-psi(indices)+pi/2), 'Color', color_nom, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        quiver(eta_cbf(1,indices), eta_cbf(2,indices), arrow_length*cos(-psi_safe(indices)+pi/2), arrow_length*sin(-psi_safe(indices)+pi/2), 'Color', color_cbf, 'MaxHeadSize', 0.5, 'AutoScale', 'off')
        
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

