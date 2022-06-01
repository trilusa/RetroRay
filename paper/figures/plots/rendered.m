% This is the simulation script to calculate the alpha parameter from
% RETRO: Retroreflector based Visible Light Indoor Localization for Real-time Tracking of IoT Devices
% by ray tracing method.

% Code is not optimized for performance, results may vary
clc; clear all; 
res=load('alpha_PD10mm_H1.5m_Rd50mm_L200mm_N1000000(1).mat');
source_pts  = cell2mat(res.detected(48));

%% Set simulation Paramters
% Meta Params
    N = 1e3;         % number of rays to cast per trial
    N_to_plot = 1e3; % number to plot if debug 
    debug = true;   % will plot and print if true

% Scene Params
    V = res.V;%0:0.1:.4;     % horizontal displacements to simulate
    R_test_pos = [-V; zeros(2, length(V))]; % (v,0,0) positions to test retroreflector
    h = res.Htemp;         % height of lamp above retrorefltor (m)

% Retroreflector Params
r=1;
    R_d = 50e-3%r.Rtemp;       % front face diamter (r in the paper, 50mm)
    R_ri = 1;        % refractive index of CC material
    R_L = .0357;     % Length of CC portion of retroreflector (35.7mm)
    R_Ls = .0063;    % Length of recession from front face to top of CC (6.3 mm)
    R_az = deg2rad(0); %azimuth angle wrt +x axis
    R_el = deg2rad(0);

% Detector Params
    D_d = res.Dtemp;     % detector diameter (1mm to simulate point source)
    D_pos = [0;0;h]; % position of center of detecor, in the center of Light
    D_norm = [0;0;-1]; % unit normal (straight down)

% Light Params
    L_od = res.L_do;     %  outer diamter of light source [100mm (front face * 2)+ largest PD]
    L_do=L_od;
    L_ro =L_do/2;
    L_id = D_d;        % inner diameter (0 to make disk source)
    L_ri = L_id/2;
    L_pos = [0;0;h]; % position of center of light, straight up by h
    L_norm = [0;0;-1]; % light unit normal (straight down)


 
%% Setup scene
light = Light(L_pos, L_norm, L_id, L_od);
detector = Detector(D_pos, D_norm, D_d);
ray_log = cell(N_to_plot,6,length(1));

%% Main raytracing loop
for v=48 % for every x position of reflector
    logged_ray_cnt=0;
    [refl1, refl2, refl3, cylender, circle] = buildCornerCube(R_d, R_L, R_Ls, R_test_pos(:,v), R_az, R_el);
    tic
    for n = 1:N
                
               x = (rand(1,1)-0.5)*L_do;
            y = (rand(1,1)-0.5)*L_do;
            while( (x^2+y^2 > L_ro^2) ||  (x^2+y^2 < L_ri^2))
                x = (rand(1,1)-0.5)*L_do;
                y = (rand(1,1)-0.5)*L_do;
            end
            source_pt  = [x y h];

            x = (rand(1,1)-0.5)*R_d(r);
            y = (rand(1,1)-0.5)*R_d(r);
            while (x^2+y^2 > (R_d(r)/2)^2)
                x = (rand(1,1)-0.5)*R_d(r);
                y = (rand(1,1)-0.5)*R_d(r);
            end
            dest_pt = [x y 0] + R_test_pos(:,v)'; %offset by v in the x direction 
%             all_rays{d,trial} = [all_rays{d,trial}; [source_pt, dest_pt]];
            
            
            %calculate the ray normal so ray object can be generated
            source_ray_norm = (dest_pt - source_pt)/norm(dest_pt - source_pt);
            ray = Ray(source_pt',source_ray_norm');
        if true
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
                else
                    disp("error: bad return ray");
                end
    
                if debug && logged_ray_cnt <= N_to_plot
%                       if best_old_ray.type == "MISSED" || best_old_ray.type == "DETECTED"
                        ray_log{n,k,1} = best_old_ray;         
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
clf
if debug
    trial_to_plot = 1;
    clf
    hold on
    axis("equal")
    grid on
    s=1.5;
    xlim(s*[-1 1]);
    ylim(s*[-1 1]);
    zlim(s*[-1 1]);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    view(deg2rad(30),deg2rad(30));
    
    % plot first reflector position
    [refl1, refl2, refl3, cylender, circle] = buildCornerCube( ...
        R_d, R_L, R_Ls, R_test_pos(:, 48), R_az, R_el);

    detector.plot();
    light.plot();
    refl1.plot();
    refl2.plot();
    refl3.plot();
    circle.plot();
    cylender.plot();

    colors = ['b','r','g','r','m','c'];
    for v=1:length([ray_log{200:800,1}])
        for j=[ 1 2 3]%length([ray_log{1,}])
            if ~isempty(ray_log{v,j})% && ray_log{i,j}.type ~= "ABSORBED" && ray_log{i,j}.type ~= "MISSED"
                ray_log{v,j}.color = colors(j);
                plot(ray_log{v,j});
            end
        end
    end
    hold off
end

