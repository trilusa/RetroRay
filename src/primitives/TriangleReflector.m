classdef TriangleReflector
    properties
        A; B; C; 
        normal;
        color;
    end
    methods
        function obj = TriangleReflector(A, B, C)
            if ~(all(size(A) == [3 1]) && all(size(B) == [3 1])  && all(size(C) == [3 1]))
                error('Must use 3x1 column vectors');
            end
            obj.A = A;
            obj.B = B;
            obj.C = C;

            obj.normal = -cross(C-A,B-A);
            obj.normal = obj.normal/norm(obj.normal);

            obj.color = [232, 166, 67]/norm([232, 166, 67]);
        end

        function [new_ray, old_ray] = intersect(obj, old_ray)
           a = obj.A;
           b = obj.B;
           c = obj.C;
          
            %calculate intersection point 
            t = (dot(obj.normal, a) - dot(obj.normal, old_ray.src_pt)) / ...
                dot(obj.normal, old_ray.dir);
            int_pt = old_ray.src_pt + old_ray.dir * t;

            %test if intesction point is inside triangle
            X = [dot(obj.normal, cross(b-a, int_pt-a))
                 dot(obj.normal, cross(c-b, int_pt-b)) 
                 dot(obj.normal, cross(a-c, int_pt-c))];
            isInside = all(X >= 0);
           
            if isInside
                if t < .00001 % ray started on surface
                    old_ray.type = "NULL";
                    new_ray = old_ray;
                else
                    old_ray.tof = t;
                    old_ray.type = "REFLECTED";
                    new_dir = (old_ray.dir - 2*obj.normal*(dot(old_ray.dir, obj.normal)));
                    new_ray = Ray(int_pt, new_dir);
                end
            else  % its a miss, return with tof=inf
                old_ray.tof = inf;
                old_ray.type = "MISSED";
                new_ray = old_ray;
            end
        end
        

        function p = plot(obj)
            % Make patch
            all = [obj.A obj.B obj.C]';
            p = patch(all(:,1), all(:,2), all(:,3), obj.color,'facealpha',0.5);    
        end
    end
end