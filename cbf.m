classdef cbf

    properties

        dyn

        h  
        Lf_h 
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
            obj.Lf2_h= fhs{3}; 
            obj.LgLf_h= fhs{4}; 

            obj.K_alpha = [1 2]; 
   
        end

        function fhs = create_fhs_for_2order_hi(obj, f, g, h, z)

            %compute lie derivatives
            Lf_h = simplify((gradient(h,z).')*f); 
            grad_Lf_h = gradient(Lf_h, z).'; 

            Lf2_h = simplify(grad_Lf_h*f); 
            LgLf_h = simplify(grad_Lf_h*g); 

            h = matlabFunction(h, 'Vars', {z}); 
            Lf_h = matlabFunction(Lf_h, 'Vars', {z}); 
            Lf2_h = matlabFunction(Lf2_h, 'Vars', {z}); 
            LgLf_h = matlabFunction(LgLf_h, 'Vars', {z}); 

            fhs = {h; Lf_h; Lf2_h; LgLf_h};

            res = obj.check_initial_conditions(); 
        end

        function bound = compute_bound(obj, z, tau)
            bound = obj.Lf2_h(z) + obj.LgLf_h(z)*tau + obj.K_alpha*[obj.h(z); obj.Lf_h(z)]; 
        end

        function result = check_initial_conditions(obj)
            %% TODO: fill inn
            result = true; 
        end

    end

end