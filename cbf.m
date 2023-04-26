classdef cbf
    properties
        k1
        k2
        k3
        k4
        k5
        k6
        d

        theta

    end

    methods(Access = public)

        function obj = cbf()

            obj.k1 = 1; 
            obj.k2 = 1; 
            obj.k3 = 1; 
            obj.k4 = 1; 
            obj.k5 = 1; 
            obj.k6 = 1; 
            obj.d = sqrt(obj.k5^2 + obj.k6^2);
          
            obj.theta = obj.compute_theta(); 
            
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

        function sig = sig_func(obj, x)
            sig= 1/(1+exp(-x)); 
        end

        function f = f_func(obj, x)
            f = obj.k1 + obj.k1*(sig((x+obj.k3)/obj.k4) - sig((x-obj.k3)/obj.k4));
        end

        function extremum = compute_extremum(obj, eta, i)
            x = eta(1) + k1*cos(eta(3) + theta1);
            y = eta(2) + k1*sin(eta(3) + theta1);
            extremum = [x, y]; 
        end

        function val_hi = h1(obj, eta, i)
            [x, y] = obj.compute_extremum(eta, i); 
            val_hi = obj.f(x) - y; 
        end

    end
end