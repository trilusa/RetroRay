% clc; 
% clear all; clf; disp(['Simulation started on: ' char(datetime)]); pwd
% %% Set simulation Paramters
%  %parpool(2);
% % Meta Params
% N=5e4;
% N_to_plot=N;
% debug=false;
% % Scene Params
% % v = 0:0.05:.8;    % horizontal displacements to simulate
% % V = 0:0.025:1;    % horizontal displacements to simulate
% % v = [v; v];
% % v = v(:)';
% V=.25;
% R_test_pos = [-V; zeros(2, length(V))]; % (v,0,0) positions to test retroreflector
% 
% % % Sim 1
% % D_d = [1e-3 10e-3 0.1e-3 ]
% % R_d = .050       % front face diamter (r in the paper, 50mm)
% % R_L = .0357     % Length of CC portion of retroreflector (35.7mm)
% % R_Ls = .0063    % Length of recession from front face to top of CC (6.3 mm)
% % H = [1.5]
% 
% % % %sim 2
% % %  D_d = [10e-3]
% % %  R_d = [ 1e-3 10e-3];% 50e-3 ];    
% % %  R_L = [ .0357/50 .0357/5];% .0357 ];     
% % %  R_Ls = [.0063/50 .0063/5];% .0063 ];
% % %  H = [1.5];
% 
% % %Sim 3
%   D_d = [10e-3]
%   R_d = 50e-3       % front face diamter (r in the paper, 50mm)
%   R_L = .0357     % Length of CC portion of retroreflector (35.7mm)
%   R_Ls = .0063    % Length of recession from front face to top of CC (6.3 mm)
%   H = 1.5          % height of lamp above retrorefltor (m)
% 
% 
% R_az = deg2rad(0); %azimuth angle wrt +x axis
% R_el = deg2rad(0);
% R_ri = 1;        % refractive index of CC material
%%
% elapsed = zeros(length(D_d),length(V));        % tracks run time for each trial 
% detector_hits = zeros(1,length(V));  % tracks number of hits per trial
detected = cell(1,length(V));
ray_log = cell(N_to_plot,6,length(V));
% all_rays = cell(length(D_d),length(v));

for r=1:length(R_d)
for h=1:length(H)
    D_pos = [0;0;H(h)]; % position of center of detecor, in the center of Light
    D_norm = [0;0;-1]; % unit normal (straight down)D_pos = [0;0;h]; % position of center of detecor, in the center of Light
for d=1:length(D_d)
      
    %calculate bounaried of light soruce for monte carlo
%     L_ro = R_d(r) + D_d(d);%R_d; %outer radius of light
    L_do = D_d(d)+2*R_d(r);
    L_ro = L_do/2;
    L_ri = D_d(d)/2; %inner radios of the light srouce (cuts out PD)
        
    disp(['*** PD= ' num2str(D_d(d)*1000) 'mm RR=' num2str(R_d(r)*1e3) 'mm H=' num2str(H(h)) 'm ***'])
    tic;
    for v = 1:length(V)
        disp(['Thread ' num2str(v) ': v = ' num2str(V(v)) 'm']);

        %generate scene components
        detector = Detector(D_pos, D_norm, D_d(d),"NoInvert"); %detector obj for ray tracer
        [refl1, refl2, refl3, cylender, circle] = buildCornerCube(R_d(r), R_L(r), R_Ls(r), [0;0;0] + R_test_pos(:,v), R_az, R_el);
       
        %% raytrace algorithm
        k=1;
        for n = 1:N
%             generate the source pt on light ring       
            x = (rand(1,1)-0.5)*L_do;
            y = (rand(1,1)-0.5)*L_do;
            while( (x^2+y^2 > L_ro^2) ||  (x^2+y^2 < L_ri^2))
                x = (rand(1,1)-0.5)*L_do;
                y = (rand(1,1)-0.5)*L_do;
            end
            source_pt  = [x y H(h)];


        % generate source point on the detector (reverse ray tracing)
%             x = (rand(1,1)-0.5)*D_d(d);
%             y = (rand(1,1)-0.5)*D_d(d);
%             while( (x^2+y^2) > (D_d(d)/2)^2 )
%                 x = (rand(1,1)-0.5)*D_d(d);
%                 y = (rand(1,1)-0.5)*D_d(d);
%             end
%             source_pt  = [x y H(h)-.001];
            
            % generate point on retroreflector
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
            
            done = false;
            i=1;
            while(~done)
                %% intersect each surface with the incident ray, add all to hits
                [new_ray1, old_ray1] = refl1.intersect(ray);
                [new_ray2, old_ray2] = refl2.intersect(ray);
                [new_ray3, old_ray3] = refl3.intersect(ray);
                [new_ray4, old_ray4] = circle.intersect(ray);
%                 [new_ray5, old_ray5] = cylender.intersect(ray);
                [new_ray6, old_ray6] = detector.intersect(ray);
                
                hits = [new_ray1, old_ray1; ...
                    new_ray2, old_ray2; ...
                    new_ray3, old_ray3; ...
                    new_ray4, old_ray4; ...
%                     new_ray5, old_ray5; ...
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
                    last_ray_dest = best_old_ray.src_pt + best_old_ray.dir * best_old_ray.tof;
                    detected{v} = [detected{v}; [source_pt, last_ray_dest']];
                    ray_log{n,i} = best_old_ray;         
                        error(n)

%                     detector_hits(d,v) = detector_hits(d,v) + 1;
                    
                else
                    disp("error: bad return ray");
                end

                if n <= N_to_plot
%                       if best_old_ray.type == "DETECTED"
                        ray_log{n,i} = best_old_ray;  
%                      end
                end
                i=i+1;
            end
            
        end
        %elapsed(d,v) = toc;
    end
    %%
    suffix=1;
    fn =['alpha_PD' num2str(D_d(d)*1000) 'mm_H' num2str(H(h)) 'm_Rd' num2str(R_d(r)*1000) 'mm_L' num2str(L_do*1000) 'mm_N' num2str(N), '(' num2str(suffix) ').mat'];
    while exist(fn, 'file') == 2 
        suffix = suffix+1;
    fn =['alpha_PD' num2str(D_d(d)*1000) 'mm_H' num2str(H(h)) 'm_Rd' num2str(R_d(r)*1000) 'mm_L' num2str(L_do*1000) 'mm_N' num2str(N), '(' num2str(suffix) ').mat'];
    end
    
%     detected_this_trial = cell(1,length(V));
%     for k = 1:length(V)
%         detected_this_trial{k} = detected{d,k};
% %         num_detected_this_trial = detector_hits(d,k);
%     end
    Rtemp=R_d(r);
    Dtemp=D_d(d);
    Htemp=H(h);
    save(fn, 'N', 'V', 'detected','Rtemp','Dtemp','Htemp', 'L_do');
    toc;
end
end
end

%% Plotting function
% Light Params
%     L_od = .1;     %  outer diamter of light source [100mm (front face * 2)+ largest PD]
%     L_id = D_d;        % inner diameter (0 to make disk source)
    L_pos = [0;0;1.5]; % position of center of light, straight up by h
    L_norm = [0;0;-1]; % light unit normal (straight down)
    
    trial_to_plot = 1;
    clf
    hold on
    axis("equal")
    grid minor
    
    xlim([-2.5 1.15]);
    ylim([-1.2 1]);
    zlim([-1.1 1.6]);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    
    light = Light(L_pos, L_norm, L_ri*2, L_do);
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

    colors = ['b','k','k','r','m','c'];
    q=10000;
    for i=44%length([ray_log{:,1,trial_to_plot}])
        if ~isempty(ray_log{i,4})
            if ray_log{i,4}.type == "DETECTED" 
            for j=1:6%length([ray_log{1,}])
                if ~isempty(ray_log{i,j,trial_to_plot})%% && ray_log{i,j}.type == "DETECTED" %&& ray_log{i,j}.type ~= "MISSED"
                    ray_log{i,j,trial_to_plot}.color = colors(j);
                    plot(ray_log{i,j,trial_to_plot});
                    
                end
            end
        end
        end
    end
    hold off
    view(180,0);