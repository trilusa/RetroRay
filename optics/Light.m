classdef Light
    properties
        center_pt; %center of light, enforced as [0;0;z]
        normal; %light plane surface normal, enforced as [0;0;-1]

        %Im modeling the light as a ring region
        %set both radii to zero for a point source
        id; %inner diameter (where PD goes inside)
        od; %outer diamter
    end

    methods
        %constructor
        function obj = Light(cp, n, id, od)
            obj.center_pt=cp;
            if all(n/norm(n) ~= [0; 0; -1])
                error('Does not support tilted sources yet');
            end
            if all(cp(1:2)~=0)
                error('Does not suppot sources not on z-axis');
            end

            obj.normal=n/norm(n);
            obj.id=id;
            obj.od=od;
        end

        function [ray, az, el] = genRandomRay(obj)
            % generates src_pt as random point on light ring and a dir vector
            % according to the lambertian ligth source model

            %calc random r withing light region
            % (outer radius > r > inner radius)
            r = obj.id/2 + (obj.od-obj.id)/2 *rand;

            %calc a random angle wrt to center of light, t
            % then convert to cartesian
            th = 2*pi*rand();
            x = r*cos(th);
            y = r*sin(th);

            %put it into a vector with the appropriate z (light height)
            src_pt = [x; y; obj.center_pt(3)];

            % generate random ray direction according to lambertian distro.
            % using inverse transform smapling
            % https://en.wikipedia.org/wiki/Inverse_transform_sampling


%             el=(pi/32)*rand(); %for debug
            el = asin(sqrt(rand()));

            %random azimuth
            az = 2*pi*rand();

            %change to cartesian vect representing direction
            [x,y,z] = sph2cart(az, (-pi/2)-el, 1);
            dir = [x;y;z];

            %construct ray and return it
            ray = Ray(src_pt, dir);
        end

        function plot(obj)
      
                w = null(obj.normal'); % Find two orthonormal vectors which are orthogonal to v
                C = obj.center_pt;
                % Make inner and outer boundaries
                t = linspace(0,2*pi);
                rin = obj.id/2;
                rout = obj.od/2;

                xin= C(1,1)+rin*cos(t)*w(1,1)+rin*sin(t)*w(1,2);
                yin= C(2,1)+rin*cos(t)*w(2,1)+rin*sin(t)*w(2,2);
                zin = C(3,1)+rin*cos(t)*w(3,1)+rin*sin(t)*w(3,2);

                xout= C(1,1)+rout*cos(t)*w(1,1)+rout*sin(t)*w(1,2);
                yout= C(2,1)+rout*cos(t)*w(2,1)+rout*sin(t)*w(2,2);
                zout = C(3,1)+rout*cos(t)*w(3,1)+rout*sin(t)*w(3,2);


                % Make patch
                patch([xout,xin],[yout,yin],[zout,zin],'y','facealpha',0.75);

        end
    end
end