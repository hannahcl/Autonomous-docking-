classdef cbf
    %consider creating an instances of cbf for each h, instead of defining
    %each h in cbf. 
    properties

        dyn

        h1_fh  
        Lf_h1_fh  
        Lf2_h1_fh  
        LgLf_h1_fh 


        h2_fh  
        Lf_h2_fh  
        Lf2_h2_fh  
        LgLf_h2_fh 


        h3_fh  
        Lf_h3_fh  
        Lf2_h3_fh  
        LgLf_h3_fh 

        h4_fh  
        Lf_h4_fh  
        Lf2_h4_fh  
        LgLf_h4_fh 

        ho1_fh  
        Lf_ho1_fh  
        Lf2_ho1_fh  
        LgLf_ho1_fh 

        ho2_fh  
        Lf_ho2_fh  
        Lf2_ho2_fh  
        LgLf_ho2_fh 

        K1_alpha 
        K2_alpha
        K3_alpha
        K4_alpha
        Ko1_alpha
        Ko2_alpha


    end



    methods(Access = public)

        function obj = cbf()

            obj.dyn = ShipDynamics(); 

            z = sym('z', [6 1]);

            h1(z) = [0 -1 0 0 0 0]*z -3; 
            fhs = obj.create_fhs_for_2order_hi(h1, z); 
            obj.h1_fh = fhs{1}; 
            obj.Lf_h1_fh = fhs{2}; 
            obj.Lf2_h1_fh = fhs{3}; 
            obj.LgLf_h1_fh = fhs{4}; 


            h2(z) = [1 0 0 0 0 0]*z + 1;
            fhs = obj.create_fhs_for_2order_hi(h2, z); 
            obj.h2_fh = fhs{1}; 
            obj.Lf_h2_fh = fhs{2}; 
            obj.Lf2_h2_fh = fhs{3}; 
            obj.LgLf_h2_fh = fhs{4}; 

            h3(z) = [-1 0 0 0 0 0]*z +1; 
            fhs = obj.create_fhs_for_2order_hi(h3, z); 
            obj.h3_fh = fhs{1}; 
            obj.Lf_h3_fh = fhs{2}; 
            obj.Lf2_h3_fh = fhs{3}; 
            obj.LgLf_h3_fh = fhs{4}; 

            h4(z) = [0 -1 0 0 0 0]*z;
            fhs = obj.create_fhs_for_2order_hi(h4, z); 
            obj.h4_fh = fhs{1}; 
            obj.Lf_h4_fh = fhs{2}; 
            obj.Lf2_h4_fh = fhs{3}; 
            obj.LgLf_h4_fh = fhs{4};

            ho1(z) = -z(1) +z(2)*atan(z(3) - pi/4); 
            fhs = obj.create_fhs_for_2order_hi(ho1, z); 
            obj.ho1_fh = fhs{1}; 
            obj.Lf_ho1_fh = fhs{2}; 
            obj.Lf2_ho1_fh = fhs{3}; 
            obj.LgLf_ho1_fh = fhs{4};

            ho2(z) = z(1) +z(2)*atan(-z(3) - pi/4); 
            fhs = obj.create_fhs_for_2order_hi(ho2, z); 
            obj.ho2_fh = fhs{1}; 
            obj.Lf_ho2_fh = fhs{2}; 
            obj.Lf2_ho2_fh = fhs{3}; 
            obj.LgLf_ho2_fh = fhs{4};

            obj.K1_alpha = [1 2]; 
            obj.K2_alpha = [1 2];
            obj.K3_alpha = [1 2];
            obj.K4_alpha = [1 2];
            obj.Ko1_alpha = [1 2];
            obj.Ko2_alpha = [1 2];

   
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