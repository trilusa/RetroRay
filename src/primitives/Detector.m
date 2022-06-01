classdef Detector
    properties
        center_pt; %pos
        normal; %normal
        r; %radius
        hits;%hits
        invert;
    end

    methods
        %% constructor
        function obj = Detector(cp,n,d,op)
            obj.center_pt = cp;
            obj.normal = n/norm(n);
            obj.r=d/2;
            obj.hits=0;
            if op=="InvertDetection"
                obj.invert = -1;
            else
                obj.invert = 1;
            end
        end

        function hits = hit(obj)
            obj.hits = obj.hits + 1;
            hits = obj.hits;
        end

        %% plot
        % render a circle at point and normal

        % make a circle
        function [new_ray, old_ray] = intersect(obj, ray)
            ray_dir = ray.dir;
            ray_src_pt = ray.src_pt;
            obj_normal = obj.normal;
            obj_center_pt = obj.center_pt;

            %find intersection point with plane
            t = (dot(obj_normal,obj_center_pt) - dot(obj_normal,ray_src_pt)) / dot(obj_normal, ray_dir);
            int_pt = ray_src_pt + ray_dir * t;

%             plot3(int_pt(1), int_pt(2), int_pt(3), 'x'); %for debug

            %check if inside circle, return new ray accordingly
            old_ray=ray;
            if (obj.invert*norm(int_pt-obj.center_pt)) < (obj.r*obj.invert) 
                new_ray = old_ray;
                old_ray.tof=t;   
                old_ray.type = "DETECTED";
                new_ray.type = "NULL";
            else
                old_ray.tof = inf;
                old_ray.type = "MISS";
                new_ray = old_ray;
            end
        end

        function plot(obj)

            w = null(obj.normal'); % Find two orthonormal vectors which are orthogonal to v
            C = obj.center_pt;
            % Make inner and outer boundaries
            t = linspace(0,2*pi);
            rin = obj.r;

            xin= C(1,1)+rin*cos(t)*w(1,1)+rin*sin(t)*w(1,2);
            yin= C(2,1)+rin*cos(t)*w(2,1)+rin*sin(t)*w(2,2);
            zin = C(3,1)+rin*cos(t)*w(3,1)+rin*sin(t)*w(3,2);

            patch(xin,yin,zin,'k','linestyle','none','facealpha',0.50);

        end
    end
end