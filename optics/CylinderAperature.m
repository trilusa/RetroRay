classdef CylinderAperature
    properties
        height;
        radius;
        normal;
        center;
        color;
        % utils
        Rx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
        Ry = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
        Rz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
    end

    methods
        function obj = CylinderAperature(radius, height, center, normal)
            if ~(all(size(radius) == [1 1]) && all(size(height) == [1 1]) && all(size(center) == [3 1])  && all(size(normal) == [3 1]))
                error('Must use 3x1 column vectors');
            end
            obj.height = height;
            obj.radius = radius;
            obj.center = center;
            obj.normal = normal / norm(normal);
            obj.color = 'r';%[232, 166, 67]/norm([232, 166, 67]);
        end

        function plot(obj)
           disp( obj.radius )
            [X, Y, Z] = cylinder2P(obj.radius, 100, (obj.center - [0;0;.025])', (obj.center + obj.normal * obj.height)');
            surf(X,Y,Z,"FaceAlpha",.2, "EdgeAlpha",0);
        end

        % calculate azimuth and elevation from cartesion surface normal


        %subtract from 90 deg for proper ref acix

        %transform the ray so the intersection with the unit cylender of
        % h=1,r=1, centered around y axis is the same raletive to the
        % actual cylender orientation

        %rotate

        % solve quadratic to find intersection points with the sides of
        % cylinder. Filter out

        % the complex t's were filtered out

        %it did hit the infinite cylender, so we will calculate the exact
        %intersection point and see if it hit actual cylender
        function [new_ray, old_ray] = intersect(obj, ray)
            
%           convert normal into sperical coordinated, rotate to match my
%           coord systems convention (down = 0 elev
            [a,e,~] = cart2sph(obj.normal(1),obj.normal(2),obj.normal(3));
            e=(pi/2)-e;

            % rotating back so as to hit unit cylender
            sp_hat = obj.Ry(-e)*obj.Rz(-a)*(ray.src_pt-obj.center);
            dir_hat = obj.Ry(-e)*obj.Rz(-a)*ray.dir;

            % scale 
            sp_hat(1:2) = sp_hat(1:2)/obj.radius;
            sp_hat(3) = sp_hat(3)/obj.height;
            dir_hat(1:2) = dir_hat(1:2)/obj.radius;
            dir_hat(3) = dir_hat(3)/obj.height;

            %solve for intersections on inf cylender, filter out negative
            %and imaginary solutions
            t=roots([dir_hat(1)^2+dir_hat(2)^2, ...
                2*(sp_hat(1)*dir_hat(1) + sp_hat(2)*dir_hat(2)), ...
                sp_hat(1)^2+sp_hat(2)^2-1]);
            t=t(imag(t)==0);
            t=t(t>0);


            old_ray=ray;
            if isempty(t)
                old_ray.tof = inf;
                old_ray.type = "MISSED";
                new_ray = old_ray;
            else

                int_pts = ray.src_pt + t'.*ray.dir;
                %                 int_pts = obj.Rz(a)*obj.Ry(e)*int_pts
                %plot points of intersection on inf cylender (for debug)
                %                 plot3(int_pts(:,1), int_pts(:,2), int_pts(:,3), '*',LineStyle='none', Color='c');

                %calculate distance along cylender axis from center point of bottom face of the
                %cylender (obj.center)
                %                 obj.normal*obj.height

                %                 int_pts-obj.center
                
                %distance from center along axiz
                dfca = ((obj.normal)' * (int_pts-obj.center));

                %filter out any t that are not on cyl surf
                t = t( (dfca>0)  & (dfca < obj.height));
                t = min(t); %take smallest t (the one that hit the closest surf)

                if isempty(t) %both inf cylinder hits not on fin cyling
                    old_ray.tof = inf;
                    old_ray.type = "MISSED";
                    new_ray = old_ray;
                else %there is a hit, calc the hit point and length of the preceding ray
                    new_ray=old_ray;
                    old_ray.tof = t;  
                    old_ray.type = 'ABSORBED'; % absorbed
                    new_ray.type = "NULL";
                end

                %                 ray.end_pt = int_pt;
                %                 ray.length = norm(ray.end_pt - ray.src_pt);
                %
            end

            %             v=ray.dir;
            %             p=ray.src_pt;
            %             va=obj.normal;
            %             pa=obj.center;
            %             dp=p-pa;
            %
            %             A= (v-dot(v,va)*va).^2
            %             B= [2 2 2]*dot(v-dot(v,va)*va, dp-dot(dp,va)*va)
            %             C= (dp-dot(dp,va)*va).^2 - obj.radius^2
            %
            %             t = roots([A; B; C])
            %             t = t(imag(t)==0)

        end
    end
end