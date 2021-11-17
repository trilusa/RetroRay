% This is the simulation script to calculate the alpha parameter from
% RETRO: Retroreflector based Visible Light Indoor Localization for Real-time Tracking of IoT Devices
% by ray tracing method.

% Code is not optimized for performance, results may vary
clc; clear all; 

%% Set simulation Paramters
% Meta Params
    N = 1e4;         % number of rays to cast per trial
    N_to_plot = 1e3; % number to plot if debug 
    debug = false;   % will plot and print if true
    trial_to_plot = 2; % which run should be recorded

% Scene Params
    v = 0:0.1:.4;     % horizontal displacements to simulate
    R_test_pos = [v; zeros(2, length(v))]; % (v,0,0) positions to test retroreflector
    h = 1.5;         % height of lamp above retrorefltor (m)

% Retroreflector Params
    R_d = .050;       % front face diamter (r in the paper, 50mm)
    R_ri = 1;        % refractive index of CC material
    R_L = .0357;     % Length of CC portion of retroreflector (35.7mm)
    R_Ls = .0063;    % Length of recession from front face to top of CC (6.3 mm)
    R_az = deg2rad(0); %azimuth angle wrt +x axis
    R_el = deg2rad(0);

% Detector Params
    D_d = .05;%0.01;     % detector diameter (1mm to simulate point source)
    D_pos = [0;0;h]; % position of center of detecor, in the center of Light
    D_norm = [0;0;-1]; % unit normal (straight down)

% Light Params
    L_od = D_d +.25;     % outer diamter of light source
    L_id = D_d;        % inner diameter (0 to make disk source)
    L_pos = [0;0;h]; % position of center of light, straight up by h
    L_norm = [0;0;-1]; % light unit normal (straight down)

%% Setup scene
light = Light(L_pos, L_norm, L_id, L_od);
detector = Detector(D_pos, D_norm, D_d);
detector_hits = zeros(length(v),1);
ray_log = cell(N_to_plot,6);
tic;
%% Main raytracing loop
for i=1:length(v) % for every x position of reflector
    [refl1, refl2, refl3, cylender, circle] = buildCornerCube(R_d, R_L, R_Ls, R_test_pos(:,i), R_az, R_el);
    for n = 1:N
        ray = light.genRandomRay();
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
            end

            if debug && n <= N_to_plot && i == trial_to_plot
                if best_old_ray.type == "DETECTED"
                    ray_log{n,k} = best_old_ray;
                    k=k+1;
                end
            end
            
        end        
    end
end

%% Plotting function
if debug
    clf
    hold on
    axis("equal")
    grid on
    s=.5;
    xlim(s*[-1 1]);
    ylim(s*[-1 1]);
    zlim(s*[-1 1]);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    view(deg2rad(30),deg2rad(30));
    
    % plot first reflector position
    [refl1, refl2, refl3, cylender, circle] = buildCornerCube( ...
        R_d, R_L, R_Ls, R_test_pos(:, trial_to_plot), R_az, R_el);

    detector.plot();
    light.plot();
    refl1.plot();
    refl2.plot();
    refl3.plot();
    circle.plot();
    cylender.plot();

    colors = ['b','k','k','r','m','c'];
    for i=1:length([ray_log{:,1}])
        for j=1:6%length([ray_log{1,}])
            if ~isempty(ray_log{i,j}) && ray_log{i,j}.type ~= "ABSORBED" && ray_log{i,j}.type ~= "MISSED"
                ray_log{i,j}.color = colors(j);
                plot(ray_log{i,j});
            end
        end
    end
    hold off
end

disp(detector_hits');
elapsed=toc
avg_ray_time = elapsed/N

fid = fopen('alpha_data.txt', 'a+');
fprintf(fid,'%d ', detector_hits);
fprintf(fid,'\n');
fclose(fid);
