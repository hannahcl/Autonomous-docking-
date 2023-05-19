classdef cbf

    properties

        dyn

        %lie derivatives
        h  
        Lf_h 
        Lf2_h  
        LgLf_h 

        K_alpha 
    end

    methods(Access = public)

        function obj = cbf(h, z)

            obj.dyn = ShipDynamics(); 

            fhs = obj.create_fhs_for_2order_hi(h, z);
            obj.h= fhs{1}; 
            obj.Lf_h = fhs{2}; 
            obj.Lf2_h= fhs{3}; 
            obj.LgLf_h= fhs{4}; 

            obj.K_alpha = [1 2]; 
   
        end

        function fhs = create_fhs_for_2order_hi(obj, h, z)

            %define f
            f(z) = [
                cos(z(3))*z(4) - sin(z(3))*z(5);
                sin(z(3))*z(4) + cos(z(3))*z(5); 
                z(6); 
                obj.dyn.A_nu(1,1)*z(4) + obj.dyn.A_nu(1,2)*z(5) + obj.dyn.A_nu(1,3)*z(6); 
                obj.dyn.A_nu(2,1)*z(4) + obj.dyn.A_nu(2,2)*z(5) + obj.dyn.A_nu(2,3)*z(6); 
                obj.dyn.A_nu(3,1)*z(4) + obj.dyn.A_nu(3,2)*z(5) + obj.dyn.A_nu(3,3)*z(6)
            ]; 

            %define g
            g= zeros(6,3); 
            g(4:6, 1:3) = obj.dyn.B_nu; 

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
        end

    end

end