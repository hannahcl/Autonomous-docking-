classdef MyClass < handle
    properties
        myProperty
    end
    
    methods
        function obj = MyClass()
            obj.myProperty = 0;
        end
        
        function obj = modifyProperty(obj, newValue)
            obj.myProperty = newValue;
        end

        function foo(obj)
            for i=1:3
                obj.myProperty
                obj = obj.modifyProperty(i); 
            end
        end
    end
end
