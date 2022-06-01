
%%
% elapsed = zeros(length(D_d),length(V));        % tracks run time for each trial 
% detector_hits = zeros(1,length(V));  % tracks number of hits per trial
detected = cell(1,length(V));
% ray_log = cell(N_to_plot,6,length(V));
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
    parfor v = 1:length(V)
        disp(['Thread ' num2str(v) ': v = ' num2str(V(v)) 'm']);

        %generate scene components
        detector = Detector(D_pos, D_norm, D_d(d),"InvertDetection"); %detector obj for ray tracer
        [refl1, refl2, refl3, cylender, circle] = buildCornerCube(R_d(r), R_L(r), R_Ls(r), [0;0;0] + R_test_pos(:,v), R_az, R_el);
       
        %% raytrace algorithm
        k=1;
        for n = 1:N
%             generate the source pt on light ring       
%             x = (rand(1,1)-0.5)*L_do;
%             y = (rand(1,1)-0.5)*L_do;
%             while( (x^2+y^2 > L_ro^2) ||  (x^2+y^2 < L_ri^2))
%                 x = (rand(1,1)-0.5)*L_do;
%                 y = (rand(1,1)-0.5)*L_do;
%             end
%             source_pt  = [x y H(h)];


        % generate source point on the detector (reverse ray tracing)
%         rn=rand(1,1)*2*pi;
%         D_r=D_d(d)/2;
%         x=D_r*cos(rn) ;
%         y=D_r*sin(rn);

            x = (rand(1,1)-0.5)*D_d(d);
            y = (rand(1,1)-0.5)*D_d(d);
            while( (x^2+y^2) > (D_d(d)/2)^2 )
                x = (rand(1,1)-0.5)*D_d(d);
                y = (rand(1,1)-0.5)*D_d(d);
            end
            source_pt  = [x y H(h)-.0001];
            
            % generate point on retroreflector
            x = (rand(1,1)-0.5)*R_d(r);
            y = (rand(1,1)-0.5)*R_d(r);
            while (x^2+y^2 > (R_d(r)/2)^2)
                x = (rand(1,1)-0.5)*R_d(r);
                y = (rand(1,1)-0.5)*R_d(r);
            end
            dest_pt = [x y 0] + R_test_pos(:,v)'; %offset by v in the x direction 
            
            
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
                [new_ray6, old_ray6] = detector.intersect(ray);
                
                hits = [new_ray1, old_ray1; ...
                    new_ray2, old_ray2; ...
                    new_ray3, old_ray3; ...
                    new_ray4, old_ray4; ...
                    new_ray6, old_ray6];
                
                %% find which which object was the one actually hit
                hits = hits([hits(:,2).type] ~= "NULL", :); %filter out any with NULL type and negative tof
                hits = hits([hits(:,2).tof] > 0, :);
                [m, idx] = min([hits(:,2).tof]); % find the row index of the least tof                 
                best_new_ray = hits(idx,1);  %save the new_ray, old_ray pair to best_hit
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
                    detected{v} = [detected{v}; [last_ray_dest', source_pt]];
%                     ray_log{n,i} = best_old_ray;                             
                else
                    disp("error: bad return ray");
                end

                if n <= N_to_plot
%                       if best_old_ray.type == "DETECTED"
%                         ray_log{n,i} = best_old_ray;  
%                      end
                end
                i=i+1;
            end
        end
    end
    %%
    suffix=1;
    fn =['alpha_PD' num2str(D_d(d)*1000) 'mm_H' num2str(H(h)) 'm_Rd' num2str(R_d(r)*1000) 'mm_L' num2str(L_do*1000) 'mm_N' num2str(N), '(' num2str(suffix) ').mat'];
    while exist(fn, 'file') == 2 
        suffix = suffix+1;
    fn =['alpha_PD' num2str(D_d(d)*1000) 'mm_H' num2str(H(h)) 'm_Rd' num2str(R_d(r)*1000) 'mm_L' num2str(L_do*1000) 'mm_N' num2str(N), '(' num2str(suffix) ').mat'];
    end
    
    Rtemp=R_d(r);
    Dtemp=D_d(d);
    Htemp=H(h);
    reversed="reversed";
    save(fn, 'N', 'V', 'detected','Rtemp','Dtemp','Htemp', 'L_do', 'reversed');
    toc;
end
end
end