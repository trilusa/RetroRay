clc; clear all; clf;
% disp("simulating...")
%% Set simulation Paramters
% Meta Params
    L_Nx = 12; %10e6;         % number of rays to cast per trial
    L_Ny = 12;
    R_N = 12;
    debug = true;
    v = 0:0.05:.15;     % horizontal displacements to simulate
    R_test_pos = [-v; zeros(2, length(v))]; % (v,0,0) positions to test retroreflector
    h = 1.5;         % height of lamp above retrorefltor (m)
    N_to_plot_per_trial = L_Nx*L_Ny*R_N^2;
% Detector Params
    D_d = 0.005;     % detector diameter (1mm to simulate point source)
    D_pos = [0;0;h]; % position of center of detecor, in the center of Light
    D_norm = [0;0;-1]; % unit normal (straight down)

% Retroreflector Params
    R_d = .050;       % front face diamter (r in the paper, 50mm)
    R_ri = 1;        % refractive index of CC material
    R_L = .0357;     % Length of CC portion of retroreflector (35.7mm)
    R_Ls = .0063;    % Length of recession from front face to top of CC (6.3 mm)
    R_az = deg2rad(0); %azimuth angle wrt +x axis
    R_el = deg2rad(0);

    % generate retroreflector pts
    basis = -R_d/2 : R_d/(R_N-1) : R_d/2;
    [Rx,Ry] = meshgrid(basis,basis);
    Rx=Rx(:);
    Ry=Ry(:);
    F = sqrt(Rx.^2+Ry.^2) >= R_d/2; %outside retro face
    Rx(F) = [];
    Ry(F) = [];

% generate light points
    L_sx = 2*R_d; %side length for light
    L_sy = R_d*2;
    basisX = -R_d : L_sx/(L_Nx-1) : L_sx-R_d;
    basisY = -L_sy/2 : L_sy/(L_Ny-1) : L_sy/2;
    [Lx,Ly] = meshgrid(basisX,basisY);
    Lx=Lx(:);
    Ly=Ly(:);
    F = sqrt(Lx.^2+Ly.^2) <= D_d/2; % remve pts inside PD
    Lx(F) = [];
    Ly(F) = [];
    source_points = [Lx, Ly, h*ones(length(Lx),1)];
    
%% Setup scene
detector = Detector(D_pos, D_norm, D_d);
detected = cell(1,length(v));
detector_hits = zeros(1,length(v));
if debug
    ray_log = cell(N_to_plot_per_trial,6,length(v));
    logged_ray_cnt = 0;
end

for trial = 1:length(v)
     tic;

    dest_points(:,:,trial) = [Rx, Ry, 0*ones(length(Rx),1)] + R_test_pos(:,trial)';
    [refl1, refl2, refl3, cylender, circle] = buildCornerCube(R_d, R_L, R_Ls, R_test_pos(:,trial), R_az, R_el);
    p=1;
    logged_ray_cnt = 0;
    %% raytrace
    disp(['Trial ' num2str(trial) ': v = ' num2str(v(trial)) 'm']);
   
    for i = 1:length(source_points(:,1))
        for j = 1:length(dest_points(:,1,trial))
            start_ray_norm = (dest_points(j,:,trial) - source_points(i,:))/norm(dest_points(j,:,trial) - source_points(i,:));
            ray = Ray(source_points(i,:)',start_ray_norm');
            first_ray = ray;
    
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
                    last_ray_dest = best_old_ray.src_pt + best_old_ray.dir * best_old_ray.tof;
                    detected{trial} = [detected{trial}; [first_ray.src_pt', last_ray_dest']];
                    detector_hits(trial) = detector_hits(trial) + 1;
                else
                    disp("error: bad return ray");
                end
    
                if debug && logged_ray_cnt <= N_to_plot_per_trial
    %                   if best_old_ray.type == "MISSED" || best_old_ray.type == "DETECTED"
                        ray_log{p,k,trial} = best_old_ray;         
                        logged_ray_cnt = logged_ray_cnt + 1;
    %                  end
                end
                k=k+1;
            end
    
        p=p+1;
        end
    end
    toc;
end
%%

% scatter(Lx,Ly, 'y' ,'+');
% axis('equal');
% hold on
% % scatter(Rx,Ry, 'b' ,'o');
% scatter(detected(:,1), detected(:,2) , 'g', '.' );
% % scatter(detected(:,1), detected(:,1) , 'g', '.' );
% hold off
%%
t=tiledlayout('flow','TileSpacing','compact');
for trial = 1:length(v)
    nexttile
    hold on
    viscircles([v(trial),0],R_d/2,'Color','b', 'LineWidth',.1 );
    viscircles(D_pos(1:2)',D_d/2,'Color','b', 'LineWidth',.025 );
%     rectangle('Position',[-R_d -L_sy/2 L_sx L_sy] );
    scatter(dest_points(:,1,trial), dest_points(:,2,trial), [] ,'g','o'); %all dest pts
    scatter(source_points(:,1), source_points(:,2), [] ,[0.9290 0.6940 0.1250],'.'); %all src pts

    if ~isempty(detected{trial})
        scatter(detected{trial}(:,1), detected{trial}(:,2),'b', '.'); %source points that hit det
        scatter(detected{trial}(:,4), detected{trial}(:,5), 'r', '.'); %hits in detector
    end
   
    axis equal
    xlim([-R_d-.01 L_sx-R_d+.01]);
    ylim([-L_sy/2-.01 L_sy/2+.01]);
    title(['v = ' num2str(v(trial))])
    hold off

    if debug && (trial==2)
        nexttile
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
            R_d, R_L, R_Ls, R_test_pos(:, trial), R_az, R_el);

        detector.plot();
        refl1.plot();
        refl2.plot();
        refl3.plot();
        circle.plot();
        cylender.plot();

        colors = ['b','k','k','r','m','c'];
        for i=1:length([ray_log{:,1,trial}])
            for j=1:6%length([ray_log{1,}])
                if ~isempty(ray_log{i,j,trial})% && ray_log{i,j}.type ~= "ABSORBED" && ray_log{i,j}.type ~= "MISSED"
                    ray_log{i,j,trial}.color = colors(j);
                    plot(ray_log{i,j,trial});
                end
            end
        end
        hold off
    end
end
lgd = legend('source pt', 'detected source pt', 'PD hit' );
lgd.Layout.Tile = 'south';
txt = ['Light panel dim: ' num2str(L_sx) 'x' num2str(L_sy) 'm' ...
    '(' num2str(L_Nx) 'x' num2str(L_Ny) ' grid), PD diam = ' num2str(D_d) 'm,' ...
    '(' num2str(R_N) 'x' num2str(R_N) ' grid)'];
title(t,"Effectiv   e Reflecting Area");
subtitle(t,txt);
hold off
% nexttile
% hold on
% x=0:0.025:.4;
% Alpha_theory = calcEffectiveReflectingArea(x,h,R_d,R_ri,R_L,R_Ls);
% plot(x,Alpha_theory)
% Alpha_simulated = 
