% This is the simulation script to calculate the alpha parameter from
% RETRO: Retroreflector based Visible Light Indoor Localization for Real-time Tracking of IoT Devices
% by ray tracing method.

% Code is not optimized for performance, results may vary
clc; clear all; 

%% Set simulation Paramters
% Meta Params
    N = 5e5;         % number of rays to cast per trial
    N_to_plot = N; % number to plot if debug 
    debug = true;   % will plot and print if true

% Scene Params
    v = .250;%0:0.1:.4;     % horizontal displacements to simulate
    R_test_pos = [-v; zeros(2, length(v))]; % (v,0,0) positions to test retroreflector
    h = 1.5;         % height of lamp above retrorefltor (m)

% Retroreflector Params
    R_d = .050;       % front face diamter (r in the paper, 50mm)
    R_ri = 1;        % refractive index of CC material
    R_L = .0357;     % Length of CC portion of retroreflector (35.7mm)
    R_Ls = .0063;    % Length of recession from front face to top of CC (6.3 mm)
    R_az = deg2rad(0); %azimuth angle wrt +x axis
    R_el = deg2rad(0);

% Detector Params
    D_d = 0.01;     % detector diameter (1mm to simulate point source)
    D_pos = [0;0;h]; % position of center of detecor, in the center of Light
    D_norm = [0;0;-1]; % unit normal (straight down)

% Light Params
    L_od = .1;     %  outer diamter of light source [100mm (front face * 2)+ largest PD]
    L_id = D_d;        % inner diameter (0 to make disk source)
    L_pos = [0;0;h]; % position of center of light, straight up by h
    L_norm = [0;0;-1]; % light unit normal (straight down)

%% Setup scene
light = Light(L_pos, L_norm, L_id, L_od);
detector = Detector(D_pos, D_norm, D_d,'NoInversion');
detector_hits = zeros(length(v),1);
ray_log = cell(N_to_plot,6,length(v));
elx = deg2rad(5);
azx1 = deg2rad(357);
azx2 = deg2rad(3);

%% Main raytracing loop
for i=1:length(v) % for every x position of reflector
    logged_ray_cnt=0;
    [refl1, refl2, refl3, cylender, circle] = buildCornerCube(R_d, R_L, R_Ls, R_test_pos(:,i), R_az, R_el);
    tic
    for n = 1:N
        [ray, az, el] = light.genRandomRay();
        if ((el < elx) || (az > azx1) || (az < azx2))  
            done = false;
            k=1;
            while(~done)        
                %% intersect each surface with the incident ray, add all to hits 
                [new_ray1, old_ray1] = refl1.intersect(ray);
                [new_ray2, old_ray2] = refl2.intersect(ray);
                [new_ray3, old_ray3] = refl3.intersect(ray);
                [new_ray4, old_ray4] = circle.intersect(ray);
                [new_ray5, old_ray5] = cylender.intersect(ray);
                [new_ray6, old_ray6] = detector.intersect(ray);
    
                hits = [new_ray1, old_ray1; ...
                        new_ray2, old_ray2; ...
                        new_ray3, old_ray3; ...
                        new_ray4, old_ray4; ...
                        new_ray5, old_ray5; ...
                        new_ray6, old_ray6];
    
                %% find which which object was the one actually hit
    
                %filter out any with NULL type and negative tof
                hits = hits([hits(:,2).type] ~= "NULL", :);
                hits = hits([hits(:,2).tof] > 0, :);
    
                % find the row index of the least tof
                [m, idx] = min([hits(:,2).tof]);
    
                %save the new_ray, old_ray pair to best_hit
                best_new_ray = hits(idx,1);
                best_old_ray = hits(idx,2);
    
                %% decide what to do next based hit type
                if best_old_ray.type == "MISSED"
                    ray = best_old_ray;
                    done = true;
                elseif best_old_ray.type == "REFLECTED"
                    ray = best_new_ray;
                elseif best_old_ray.type == "ABSORBED"
                    done = true;
                elseif best_old_ray.type == "DETECTED"
                    done = true;
                    detector_hits(i) = detector_hits(i) + 1;
                else
                    disp("error: bad return ray");
                end
    
                if debug && logged_ray_cnt <= N_to_plot
%                       if best_old_ray.type == "MISSED" || best_old_ray.type == "DETECTED"
                        ray_log{n,k,i} = best_old_ray;         
                        logged_ray_cnt = logged_ray_cnt + 1;
%                      end
                end
                 k=k+1;
            end        
        else
%             if debug
%                 ray.tof = inf;
%                 ray_log{n,1,i} = ray;         
%                 logged_ray_cnt = logged_ray_cnt + 1;
%             end
        end
    end
    toc
    elapsed=toc;
end

%% Plotting function
if debug
    
    trial_to_plot = 1;
    clf
    hold on
    axis("equal")
    grid minor
    
    xlim([-2.5 .15]);
    ylim([-.2 1]);
    zlim([-.1 1.6]);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    
    light = Light(L_pos, L_norm, L_id, L_od);
    % plot first reflector position
    [refl1, refl2, refl3, cylender, circle] = buildCornerCube( ...
        R_d(1), R_L(1), R_Ls(1), R_test_pos(:, trial_to_plot), R_az, R_el);

    detector.plot();
    light.plot();
    refl1.plot();
    refl2.plot();
    refl3.plot();
    circle.plot();
    cylender.plot();

    colors = ['r','k','k','b','m','c'];
    q=10000;
    for i=q:q+8000%length([ray_log{:,1,trial_to_plot}])
        for j=1:6%length([ray_log{1,}])
            if ~isempty(ray_log{i,j,trial_to_plot}) && ray_log{i,j}.type == "DETECTED" %&& ray_log{i,j}.type ~= "MISSED"
                ray_log{i,j,trial_to_plot}.color = colors(j);
                plot(ray_log{i,j,trial_to_plot});
                
            end
        end
    end
    hold off
end
    view(180,0);

% hold on
% Use fourth input for color scale.
% patch([1 -1 -1 1]*10, [1 1 -1 -1]*10, [-1 -1 -1 -1 ]*.5, 'b')  
% p1=patch([1 -1 -1 1]*10, [-1 -1 -1 -1 ]*.5, [1 1 -1 -1]*10,'b')  
% p1.FaceAlpha=.2;
% hold off

%%
avg_ray_time = elapsed/N
disp(detector_hits');
disp('done')
% 
% fid = fopen('alpha_data.txt', 'a+');
% fprintf(fid, '%d  %d  %d  %d  %d  R_d=%.3f  D_d=%.3f  h=%.2f  N=%.2e  T=%.2fs\n', detector_hits, R_d, D_d,h, N,elapsed);
% fclose(fid);
