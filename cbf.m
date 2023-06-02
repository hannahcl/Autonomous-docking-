classdef cbf

    properties

        dyn

        h  
        Lf_h 
        Lg_h
        Lf2_h  
        LgLf_h 

        K_alpha 
    end

    methods(Access = public)

        function obj = cbf(f, g, h, z)

            obj.dyn = ShipDynamics(); 

            fhs = obj.create_fhs_for_2order_hi(f, g, h, z);
            obj.h= fhs{1}; 
            obj.Lf_h = fhs{2};
            obj.Lg_h = fhs{3}; 
            obj.Lf2_h= fhs{4}; 
            obj.LgLf_h= fhs{5}; 

            obj.K_alpha = [1 2]; 
   
        end

        function fhs = create_fhs_for_2order_hi(obj, f, g, h, z)
            %compute lie derivatives
            Lf_h = simplify((gradient(h,z).')*f);
            Lg_h = simplify((gradient(h,z).')*g);
            grad_Lf_h = gradient(Lf_h, z).'; 

            Lf2_h = simplify(grad_Lf_h*f); 
            LgLf_h = simplify(grad_Lf_h*g); 

            h = matlabFunction(h, 'Vars', {z}); 
            Lf_h = matlabFunction(Lf_h, 'Vars', {z});
            Lg_h = matlabFunction(Lg_h, 'Vars', {z});
            Lf2_h = matlabFunction(Lf2_h, 'Vars', {z}); 
            LgLf_h = matlabFunction(LgLf_h, 'Vars', {z}); 

            fhs = {h; Lf_h; Lg_h; Lf2_h; LgLf_h};
        end

        function lin_con = compute_linear_constraints(obj, z)
            A = -obj.Lg_h(z);
            b = obj.Lf2_h(z) + obj.K_alpha*[obj.h(z); obj.Lf_h(z)];
            lin_con = [A, b];  
        end

        function cbf_valid = check_initial_conditions(obj, z0, tau0)
            cbf_valid = true; 

            F = [0 1; 0 0]; 
            G = [0; 1];
            A = F - G*obj.K_alpha; 
            lambda = eig(A); 

            if ~(isreal(lambda(1)) && (lambda(1) < 0))
                cbf_valid = false; 
            end

            if ~(isreal(lambda(2)) && (lambda(2) < 0))
                cbf_valid = false; 
            end

            nu_0 = obj.h(z0); 
            nu_0_dot = obj.Lf_h(z0); 
            nu_1 = nu_0_dot - lambda(1)*nu_0; 
            nu_1_dot = obj.Lf2_h(z0) + obj.LgLf_h(z0)*tau0 -lambda(1)*obj.Lf_h(z0);

            if ~(nu_0 >= 0)
                cbf_valid = false; 
            end

            if ~(-lambda(1) >= (nu_0_dot/nu_0))
                cbf_valid = false; 
            end

            if ~(-lambda(2) >= (nu_1_dot/nu_1))
                cbf_valid = false; 
            end

        end

        function lin_con = compute_linear_constraints_1order(obj, z)
            A = -obj.Lg_h(z); 
            b = obj.K_alpha(1)*obj.h(z); 
            lin_con = [A, b];  
        end

        function cbf_valid = check_initial_conditions_1order(obj, z0, tau0)
            cbf_valid = true; 

            F = 1; 
            G = 1;
            A = F - G*obj.K_alpha(1); 
            lambda = eig(A); 

            if ~(isreal(lambda) && (lambda < 0))
                cbf_valid = false; 
            end

            nu_0 = obj.h(z0) ;  
            nu_0_dot = obj.Lf_h(z0); 
 
            if ~(nu_0 >= 0)
                cbf_valid = false; 
            end

            if ~(-lambda(1) >= (nu_0_dot/nu_0))
                cbf_valid = false; 
            end

        end

    end

end