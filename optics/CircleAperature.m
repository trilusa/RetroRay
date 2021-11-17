classdef CircleAperature
    properties
        radius;
        center_pt;
        normal;
    end
    methods
        function obj = CircleAperature(radius, center_pt, normal)
            if ~(all(size(radius) == [1 1]) && all(size(center_pt) == [3 1])  && all(size(normal) == [3 1]))
                error('Must use 3x1 column vectors');
            end
            obj.radius = radius;
            obj.center_pt = center_pt;
            obj.normal = normal;
        end

        function [new_ray, old_ray] = intersect(obj, ray)

            %find intersection point with plane
            t = (dot(obj.normal,obj.center_pt) - dot(obj.normal,ray.src_pt)) / dot(obj.normal, ray.dir);
            int_pt = ray.src_pt + ray.dir * t;

%           plot3(int_pt(1), int_pt(2), int_pt(3), '*'); %for debug
           
            %check if inside circle, return new ray accordingly
            old_ray=ray;
            if norm(int_pt-obj.center_pt) < obj.radius    
                old_ray.tof = inf;
                old_ray.type = 'MISSED'; 
                new_ray = old_ray;
            else
                new_ray=old_ray;
                old_ray.tof = t;  
                old_ray.type = 'ABSORBED'; % absorbed
                new_ray.type = "NULL";
            end
        end
        

        function p = plot(obj)
            w = null(obj.normal'); % Find two orthonormal vectors which are orthogonal to v
            C = obj.center_pt;

            % Make inner and outer boundaries
            t = linspace(0,2*pi);
            rin = obj.radius;
            rout = obj.radius*3;

            xin= C(1,1)+rin*cos(t)*w(1,1)+rin*sin(t)*w(1,2);
            yin= C(2,1)+rin*cos(t)*w(2,1)+rin*sin(t)*w(2,2);
            zin = C(3,1)+rin*cos(t)*w(3,1)+rin*sin(t)*w(3,2);
        
            xout= C(1,1)+rout*cos(t)*w(1,1)+rout*sin(t)*w(1,2);
            yout= C(2,1)+rout*cos(t)*w(2,1)+rout*sin(t)*w(2,2);
            zout = C(3,1)+rout*cos(t)*w(3,1)+rout*sin(t)*w(3,2);

            % Make patch
            p = patch([xout,xin],[yout,yin],[zout,zin],'k','linestyle','none','facealpha',0.25);    
   
        end
    end
end