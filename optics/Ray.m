classdef Ray
    properties
        src_pt; %the starting point for the ray
        dir;  %the direction of the ray, enforced normalization
        tof; %time paremeter to get to end point (time of flight)
        color = [255, 255, 184]/norm([255, 255, 184]);
        type;
    end
    methods
        function obj  = Ray(src_pt, dir)
            %enforce column vectors
            if ~(all(size(src_pt) == [3 1]) && all(size(dir) == [3 1]))
                error('Ray data must be 3x1 column vectors');
            end
            obj.src_pt = src_pt;
            obj.dir = dir/norm(dir);
        end

        function plot(obj)
            if(obj.tof == inf)
                end_pt = obj.src_pt + 5*obj.dir; %make it really long
            elseif (obj.tof <= 0)
                error('Bad ray. Time of flight must be positive scalar for plotting. \n %s');
            elseif (obj.tof > 0)
                end_pt = obj.src_pt + obj.tof*obj.dir;
            end
            
            plot3([obj.src_pt(1); end_pt(1)], ...
                  [obj.src_pt(2); end_pt(2)], ...
                  [obj.src_pt(3); end_pt(3)], obj.color);
        end
    end
end