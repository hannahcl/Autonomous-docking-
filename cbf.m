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

            obj.k1 = 3; 
            obj.k2 = 2; 
            obj.k3 = 0.9; 
            obj.k4 = 10e-3; 
            obj.k5 = 2.5; 
            obj.k6 = 1.5; 
            obj.d = sqrt(obj.k5^2 + obj.k6^2);
          
            obj.theta = obj.compute_theta(); 
            
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