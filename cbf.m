classdef cbf
    properties

        %k1 = dedpth of the doc
        %k2 = width of teh doc
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

        h1_fh
        grad_h1_fh


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
            obj.alpha = 10; 
          
            obj.theta = obj.compute_theta(); 

            fhs = obj.create_func_handles(); 
            obj.h1_fh = fhs{1};  
            obj.grad_h1_fh = fhs{2}; 
            
        end

        function fhs = create_func_handles(obj) 
            %does not have to be publisc. will find function handle only
            %once. 

            syms x
            eta = sym('eta', [3 1]); 

            sig = @(x) 1/(1+exp(-x)); 
            f= @(x) obj.k1 + obj.k1*(sig((x+obj.k3)/obj.k4) - sig((x-obj.k3)/obj.k4)); 
             
           % theta1 = atan((obj.k5/2)/(obj.k6/2)); 
%             theta4 = theta1 + pi/2;
%             theta3 = theta1 + pi; 
%             theta2 = theta1 - pi/2;

            %create function handle in for loop here afterwards
            
%             extremum1_x(eta) = eta(1) + obj.d*cos(eta(3) + theta1); 
%             extremum1_y(eta) = eta(2) + obj.d*sin(eta(3) + theta1); 
            h1(eta) = f(eta(1)) - eta(2) -6; 
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