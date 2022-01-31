clc; 
while true
clear all; clf; disp(['Simulation started on: ' char(datetime)]);
%% Set simulation Paramters
%parpool(62);
% Meta Params
N=1e6;
debug=false;
% Scene Params
% v = 0:0.05:.8;    % horizontal displacements to simulate
v = 0:0.025:.75;    % horizontal displacements to simulate
v = [v; v];
v = v(:)';

R_test_pos = [-v; zeros(2, length(v))]; % (v,0,0) positions to test retroreflector
H = [1.5];% 3 4.5];         % height of lamp above retrorefltor (m)

% Detector Params
% D_d = 0.004:0.002:0.016;     % detector diameter (1mm to simulate point source)
% D_d = (.1:.05:4)/1000;
D_d = [.1e-3 1e-3 10e-3]  % .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 1 1 1 1 1 1 1 1 1 1]/1e3%.2 .4 .6 .8 1.0 2  4  6  8  10 ...
%12 14 16 18 20 25 30 35 40 45 50 55 60 65 70 75 100 150 200]/1000;


% Retroreflector Params
R_d = .050/5;       % front face diamter (r in the paper, 50mm)
R_ri = 1;        % refractive index of CC material
R_L = .0357/5;     % Length of CC portion of retroreflector (35.7mm)
R_Ls = .0063/5;    % Length of recession from front face to top of CC (6.3 mm)

R_az = deg2rad(0); %azimuth angle wrt +x axis
R_el = deg2rad(0);
% R_d = [ .050/50 .050/5 .050 ];    
% R_L = [ .0357/50 .0357/5 .0357 ];     
% R_Ls = [.0063/50 .0063/5 .0063 ]; 

%light boundaries

elapsed = zeros(length(D_d),length(v));        % tracks run time for each trial 
detector_hits = zeros(length(D_d),length(v));  % tracks number of hits per trial
detected = cell(length(D_d),length(v));
% all_rays = cell(length(D_d),length(v));

for r=1:length(R_d)
for h=H
    D_pos = [0;0;h]; % position of center of detecor, in the center of Light
    D_norm = [0;0;-1]; % unit normal (straight down)D_pos = [0;0;h]; % position of center of detecor, in the center of Light
for d=1:length(D_d)
    disp(['*** PD diameter = ' num2str(D_d(d)*1000) 'mm ***'])
    
    %calculate bounaried of light soruce for monte carlo
    L_ro = R_d(r) + D_d(d);%R_d; %outer radius of light
    L_ri = D_d(d)/2; %inner radios of the light srouce (cuts out PD)
    
    detector = Detector(D_pos, D_norm, D_d(d)); %detector obj for ray tracer
    
    parfor trial = 1:length(v)
        tic;  disp(['Trial ' num2str(trial) ': v = ' num2str(v(trial)) 'm']);

        %generate retroreflector components
        [refl1, refl2, refl3, cylender, circle] = buildCornerCube(R_d(r), R_L(r), R_Ls(r), [0;0;0] + R_test_pos(:,trial), R_az, R_el);
       
        %% raytrace algorithm
        
        for i = 1:N
            %generate the source pt on light
%             r=abs(L_ro*sum(rand(1,10),2)/10-L_ro/2)+L_ri;
%             r = rand(1,1)*(L_ro-L_ri) + L_ri;  %random radius on light source between PD and outer boundary
%             th = rand(1,1)*2*pi;                    %randome angle on light source
%             [x,y,z] = pol2cart(th,r,1.5);            %convert to cartestian
            
            x = (rand(1,1)-0.5)*2*L_ro;
            y = (rand(1,1)-0.5)*2*L_ro;
            while( (x^2+y^2 > L_ro^2) ||  (x^2+y^2 < L_ri^2))
                x = (rand(1,1)-0.5)*2*L_ro;
                y = (rand(1,1)-0.5)*2*L_ro;
            end
            source_pt  = [x y h];
            
            %generate the target point on retroreflector
%             r = rand(1,1)*(R_d/2);                  %random radius on retroreflector
%             th = rand(1,1)*2*pi;                    %rand angle
%             [x,y,z] = pol2cart(th,r,0); %convert to cartesian

            x = (rand(1,1)-0.5)*R_d(r);
            y = (rand(1,1)-0.5)*R_d(r);
            while (x^2+y^2 > (R_d(r)/2)^2)
                x = (rand(1,1)-0.5)*R_d(r);
                y = (rand(1,1)-0.5)*R_d(r);
            end
            dest_pt = [x y 0] + R_test_pos(:,trial)'; %offset by v in the x direction 
%             all_rays{d,trial} = [all_rays{d,trial}; [source_pt, dest_pt]];
            
            
            %calculate the ray normal so ray object can be generated
            source_ray_norm = (dest_pt - source_pt)/norm(dest_pt - source_pt);
            ray = Ray(source_pt',source_ray_norm');
            
            done = false;
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
                    detected{d,trial} = [detected{d,trial}; [source_pt, last_ray_dest']];
                    detector_hits(d,trial) = detector_hits(d,trial) + 1;
                else
                    disp("error: bad return ray");
                end
            end
        end
        toc
        elapsed(d,trial) = toc;
    end
    %%
    suffix=1;
    fn =['monte-carlo-alpha_N' num2str(N) '_PD' num2str(D_d(d)*1000) 'mm_L' num2str(L_ro*1000) 'mm_h' num2str(h) 'm_Rd' num2str(R_d(r)*1000) 'mm(' num2str(suffix) ').mat'];
    while exist(fn, 'file') == 2 
        suffix = suffix+1;
    fn =['monte-carlo-alpha_N' num2str(N) '_PD' num2str(D_d(d)*1000) 'mm_L' num2str(L_ro*1000) 'mm_h' num2str(h) 'm_Rd' num2str(R_d(r)*1000) 'mm(' num2str(suffix) ').mat'];
    end
    
    detected_this_trial = cell(1,length(v));
    for k = 1:length(v)
        detected_this_trial{k} = detected{d,k};
        num_detected_this_trial = detector_hits(d,k);
    end
    rrd=R_d(r);
    save(fn, 'detected_this_trial' ,'num_detected_this_trial','h','L_ro','N','d','rrd','v');
    disp(['--- Total time: ' num2str(sum(elapsed(d,:))/60) ' minutes ---']);
end
end
end
end
%%
%{
t=tiledlayout(length(D_d),length(v),'TileSpacing','compact');
txt = ['Light panel dim: ' num2str(L_ro) 'm OD' ...
    'h = ' num2str(h) ', '...
    'PD diam = ' num2str(D_d(d)*1000) 'mm, ' ...
    'num rays per trial: ' num2str(N) ', '...
    'total time:' num2str(sum(elapsed(:))/60) 'min' ];
title(t,["Effective Reflecting Area" txt]);

for d = 1:length(D_d)
   for trial = 1:length(v)
   nexttile
   hold on
   source_points = all_rays{d,trial}(:,1:3);
    scatter(source_points(:,1), source_points(:,2), [] ,[0.9290 0.6940 0.1250],'.'); %all source pts
    viscircles(D_pos(1:2)',R_d/2,'Color','k', 'LineWidth',.1 );
    viscircles(D_pos(1:2)',D_d(d)/2,'Color','m', 'LineWidth',.025 );
    viscircles(D_pos(1:2)',L_ro,'Color','m', 'LineWidth',.025 );
    if ~isempty(detected{d,trial})
        scatter(detected{d,trial}(:,1), detected{d,trial}(:,2),'b', '.'); %source points that hit det
        scatter(detected{d,trial}(:,4), detected{d,trial}(:,5), 'r', '.'); %hits in detector
    end 
    xlim([-L_ro L_ro])
    ylim([-L_ro L_ro])
    axis equal
    axis square
    hold off
   end
end

%%

t=tiledlayout('flow','TileSpacing','compact');

for trial = 1:length(v)
    nexttile
    hold on
    
    %     rectangle('Position',[-R_d -L_sy/2 L_sx L_sy] );
    scatter(source_points(:,1), source_points(:,2), [] ,[0.9290 0.6940 0.1250],'.'); %all source pts
    viscircles([v(trial),0],R_d/2,'Color','k', 'LineWidth',.1 );
    viscircles(D_pos(1:2)',D_d(d)/2,'Color','k', 'LineWidth',.025 );
    if ~isempty(detected{trial})
        scatter(detected{trial}(:,1), detected{trial}(:,2),'b', '.'); %source points that hit det
        scatter(detected{trial}(:,4), detected{trial}(:,5), 'r', '.'); %hits in detector
    end
    axis equal
    %     xlim([-R_d-.005 L_sx-R_d+.005]);
    ylim([-L_ro/2-.005 L_ro/2+.005]);
    xlim([-L_ro/2-.005 L_ro/2+.005]);
    %         xlim([-.0001-D_d(d)/2  .0001+D_d(d)/2])
    %         ylim([-.0001-D_d(d)/2  .0001+D_d(d)/2])
    
    txt_title=['v = ' num2str(v(trial)) 'm (' num2str(detector_hits(trial)) ' hits)'];
    %     txt_subtitle = ['time = ' num2str(elapsed(trial)/60) ' min'];
    title([convertCharsToStrings(txt_title)]);%, convertCharsToStrings(txt_subtitle )])
    if trial == 1
        lgd = legend('source pt', 'detected source pt', 'PD hit' );
        set(lgd,'location','north outside')
    end
    hold off
end
hold off
nexttile

hold on
% x=0:0.025:.4;
% Alpha_theory = calcEffectiveReflectingArea(x,h,R_d,R_ri,R_L,R_Ls);


% num_detected= zeros(1,length(v));
% for trial=1:length(v)
%     num_detected(trial) = length(detected{trial}(:,1));
% end

% calcEffectiveReflectingArea();

U=zeros(1,length(v));
for u = 1:length(v)
    if ~isempty(detected{u})
        Utemp = unique(detected{u}(:,1:3), 'rows');
        U(u) = length(Utemp(:,1));
    end
end
plot(v,U/U(1))

all = detector_hits/detector_hits(1);
plot(v,all);
legend('calculated (point PD)', ['simulated (' num2str(D_d(d)*1000) 'mm PD)'], 'total hits (include redundant src pts)')
title(['Normalized Effective Reflectiing Area (\alpha) for ' num2str(D_d(d)*1000) 'mm PD'])

%}
disp(['Simulation ended on: ' char(datetime)]);
