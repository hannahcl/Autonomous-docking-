classdef cbf
    properties

        %k1 = dedpth of the doc
        %k2 = width of the doc
        %k3 = tollerance for difference between sigmoid func and acctual doc
        %k4 = sparpnes of sigmoid
        
        %k5 = length of boat
        %k6 = width of boat

        k1
        k2
        k3
        k4
        k5
        k6
        d
        alpha

        theta

%         h1_fh
%         grad_h1_fh

        dyn
        A 

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

        K1_alpha 
        K2_alpha
        K3_alpha


    end



    methods(Access = public)

        function obj = cbf()

            obj.k1 = 3; 
            obj.k2 = 1; 
            obj.k3 = obj.k2 - 0.1; 
            obj.k4 = 10e-3; 
            obj.k5 = 2.5; 
            obj.k6 = 1.5; 
            obj.d = sqrt(obj.k5^2 + obj.k6^2);
            obj.alpha = 1; 
          
            obj.theta = obj.compute_theta(); 

            obj.dyn = ShipDynamics(); 
            obj.A = -obj.dyn.M\obj.dyn.N_lin; 

            fhs = obj.create_fhs_for_2order_h1(); 
            obj.h1_fh = fhs{1}
            obj.Lf_h1_fh = fhs{2}; 
            obj.Lf2_h1_fh = fhs{3}; 
            obj.LgLf_h1_fh = fhs{4}; 

            fhs = obj.create_fhs_for_2order_h2(); 
            obj.h2_fh = fhs{1}
            obj.Lf_h2_fh = fhs{2}; 
            obj.Lf2_h2_fh = fhs{3}; 
            obj.LgLf_h2_fh = fhs{4}; 

            fhs = obj.create_fhs_for_2order_h3(); 
            obj.h3_fh = fhs{1} 
            obj.Lf_h3_fh = fhs{2}; 
            obj.Lf2_h3_fh = fhs{3}; 
            obj.LgLf_h3_fh = fhs{4}; 

            obj.K1_alpha = [1 2]; 
            obj.K2_alpha = [1 2];
            obj.K3_alpha = [1 2];

%             fhs = obj.create_func_handles(); 
%             obj.h1_fh = fhs{1};  
%             obj.grad_h1_fh = fhs{2}; 
   
        end

        function fhs = create_fhs_for_2order_h1(obj)

            z = sym('z', [6 1]);

            %define h
%             syms x
%             sig = @(x) 1/(1+exp(-x)); 
%             doc= @(x) -obj.k1 + obj.k1*(sig((x+obj.k3)/obj.k4) - sig((x-obj.k3)/obj.k4)); 
%              
%             h(z) = doc(z(1)) - z(2); 
            %h(z) = (z(1)^2 + (-1 -z(2))^2 - 1); 
            h(z) = [0 -1 0 0 0 0]*z; 

            %define f

            f(z) = [
                cos(z(3))*z(4) - sin(z(3))*z(5);
                sin(z(3))*z(4) + cos(z(3))*z(5); 
                z(6); 
                obj.A(1,1)*z(4) + obj.A(1,2)*z(5) + obj.A(1,3)*z(6); 
                obj.A(2,1)*z(4) + obj.A(2,2)*z(5) + obj.A(2,3)*z(6); 
                obj.A(3,1)*z(4) + obj.A(3,2)*z(5) + obj.A(3,3)*z(6)
            ]; 

            %define g
            g= zeros(6,3); 
            g(4:6, 1:3) = eye(3); 

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

        function fhs = create_fhs_for_2order_h2(obj)

            z = sym('z', [6 1]);

            %define h
            h(z) = [1 0 0 0 0 0]*z + 1; 

            %define f

            f(z) = [
                cos(z(3))*z(4) - sin(z(3))*z(5);
                sin(z(3))*z(4) + cos(z(3))*z(5); 
                z(6); 
                obj.A(1,1)*z(4) + obj.A(1,2)*z(5) + obj.A(1,3)*z(6); 
                obj.A(2,1)*z(4) + obj.A(2,2)*z(5) + obj.A(2,3)*z(6); 
                obj.A(3,1)*z(4) + obj.A(3,2)*z(5) + obj.A(3,3)*z(6)
            ]; 

            %define g
            g= zeros(6,3); 
            g(4:6, 1:3) = eye(3); 

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

        function fhs = create_fhs_for_2order_h3(obj)

            z = sym('z', [6 1]);

            %define h
            h(z) = [-1 0 0 0 0 0]*z +1; 

            %define f

            f(z) = [
                cos(z(3))*z(4) - sin(z(3))*z(5);
                sin(z(3))*z(4) + cos(z(3))*z(5); 
                z(6); 
                obj.A(1,1)*z(4) + obj.A(1,2)*z(5) + obj.A(1,3)*z(6); 
                obj.A(2,1)*z(4) + obj.A(2,2)*z(5) + obj.A(2,3)*z(6); 
                obj.A(3,1)*z(4) + obj.A(3,2)*z(5) + obj.A(3,3)*z(6)
            ]; 

            %define g
            g= zeros(6,3); 
            g(4:6, 1:3) = eye(3); 

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

        function fhs = create_func_handles(obj) 

            syms x
            eta = sym('eta', [3 1]); 

            sig = @(x) 1/(1+exp(-x)); 
            f= @(x) -obj.k1 + obj.k1*(sig((x+obj.k3)/obj.k4) - sig((x-obj.k3)/obj.k4)); 
             
           % theta1 = atan((obj.k5/2)/(obj.k6/2)); 
%             theta4 = theta1 + pi/2;
%             theta3 = theta1 + pi; 
%             theta2 = theta1 - pi/2;

            %create function handle in for loop here afterwards
            
%             extremum1_x(eta) = eta(1) + obj.d*cos(eta(3) + theta1); 
%             extremum1_y(eta) = eta(2) + obj.d*sin(eta(3) + theta1); 

            h1(eta) = f(eta(1)) - eta(2); 
            grad_h1(eta) = simplify(gradient(h1, eta)'); 
  
            h1_fh = matlabFunction(h1, 'Vars', {eta}); 
            grad_h1_fh = matlabFunction(grad_h1, 'Vars', {eta}); 

            fhs = {h1_fh; grad_h1_fh}; 

    end

        function val_hi = hi(obj, eta, i)
            ext_pt = obj.compute_extremum(eta, i); 
            val_hi = obj.f(ext_pt(1)) - ext_pt(2); 
        end

        function f_val = f(obj, x)
            temp1 = (x+obj.k3)/obj.k4; 
            temp2 = (x-obj.k3)/obj.k4;

            temp3 = 1/(1+exp(temp1)); 
            temp4 = 1/(1+exp(temp2)); 
            
            f_val = -obj.k1 - obj.k1*(temp3 - temp4);
        end

    end

    methods(Access = private)

        function theta = compute_theta(obj) 
            theta1 = atan((obj.k5/2)/(obj.k6/2)); 
            theta4 = theta1 + pi/2;
            theta3 = theta1 + pi; 
            theta2 = theta1 - pi/2;  

            theta = [theta1; theta2; theta3; theta4]; 

        end



        function sig_val = sig_func(obj, x)
            sig_val= 1./(1+exp(-x)); 
        end


        function extremum = compute_extremum(obj, eta, i)
            x = eta(1) + obj.k1*cos(eta(3) + obj.theta(i));
            y = eta(2) + obj.k1*sin(eta(3) + obj.theta(i));
            extremum = [x; y]; 
        end


    end
end