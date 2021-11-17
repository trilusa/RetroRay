clc;clear all;clf;tic;
N=5e4; %change this to match the "sim_length" output of test sim
h=.2;
d=.02;%.01;
th = deg2rad([5,30,60]);
% theta0 = deg2rad(0);
% theta1 = deg2rad(5);
light = Light([0;0;h], [0; 0; -1], 0, 0);
% detector0 = Detector([-h*sin(theta0); 0; -h*(1-cos(theta0))], [-sin(theta0); 0; -cos(theta0)], 1);
% detector1 = Detector([-h*sin(theta1); 0; -h*(1-cos(theta1))], [-sin(theta1); 0; -cos(theta1)], 1);
detectors=Detector.empty(length(th),0);
for i=1:length(th)
    detectors(i) = Detector([h*sin(th(i)); 0; h*(1-cos(th(i)))], [sin(th(i)); 0; -cos(th(i))], d);
%     detectors(i).plot()
end
detector_hits = zeros(1,length(th));
ray = light.genRandomRay();

% hold on

% light.plot();
% detectors(1).plot();
% detector1.plot();

for i=1:N
    ray = light.genRandomRay();
    for j=1:length(detectors)

        [new_ray, old_ray] = detectors(j).intersect(ray);
        if old_ray.type == "DETECTED"
            detector_hits(j) = detector_hits(j) + 1;
%             old_ray.color = 'r';
%             old_ray.plot();
        else
%             old_ray.color = 'k';
%             old_ray.plot();
        end
    end
end

%%Calculating the number of iterations the for the desired sim time
disp(detector_hits)
time=toc
avg_ray_time = time/N
desired_sim_time_s = 7*60*60 %7 hours
desired_N = desired_sim_time_s/avg_ray_time

%%writing results to file
fid = fopen('lambertian_test_data.txt', 'a+');
fprintf(fid, '%d %d %d %d %d %.2e %.1f %.3f %.2f\n', detector_hits,N,h,d,time);
fclose(fid);
% 
% view(3)
% axis('equal')
% zlim([-1, h+1])
% ylim([-h-1,h+1])
% xlim([-h-1,h+1])
% hold off
