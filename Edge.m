classdef Edge
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        index
        targets
        locations
        enabled
    end
    
    methods
        function obj = Edge(index, targets, locations)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.index = index;
            obj.targets = targets;
            obj.locations = locations;
        end
        
        function outputArg = drawEdge(obj)
            hold on
            p1 = plot(obj.locations(:,1),obj.locations(:,2),'k','LineWidth',2);
            if ~obj.enabled 
                p1.Color(4) = 0.1;%0.2
            end
           
        end
        
        function outputArg = highlightEdge(obj,color)
            hold on
            p1 = plot(obj.locations(:,1),obj.locations(:,2),color,'LineWidth',3);
            if obj.enabled
                p1.Color(4) = 0.5;
            else
                p1.Color(4) = 0.2;
            end
        end
        
        function distanceOutput = getLength(obj)
            distanceOutput = norm(obj.locations(1,:)-obj.locations(2,:),2);
            if ~obj.enabled
%                 distanceOutput = 10*distanceOutput;
                distanceOutput = 100*distanceOutput;
            end
                
        end
        
    end
end

