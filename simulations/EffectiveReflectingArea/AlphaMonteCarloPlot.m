%% Script to do a fancy plot of ERA vs PD size

    N=5e6;
    suffix=1;

    v = 0:0.2:.8;  % horizontal displacements to simulate
    R_test_pos = [-v; zeros(2, length(v))]; % (v,0,0) positions to test retroreflector
    h = 1.5;         % height of lamp above retrorefltor (m)

% Detector Params
    D_d = 0.004:0.002:0.016;     % detector diameter (1mm to simulate point source)
    D_pos = [0;0;h]; % position of center of detecor, in the center of Light
    D_norm = [0;0;-1]; % unit normal (straight down)

% Retroreflector Params
    R_d = .050;       % front face diamter (r in the paper, 50mm)
    R_ri = 1;        % refractive index of CC material
    R_L = .0357;     % Length of CC portion of retroreflector (35.7mm)
    R_Ls = .0063;    % Length of recession from front face to top of CC (6.3 mm)
    R_az = deg2rad(0); %azimuth angle wrt +x axis
    R_el = deg2rad(0);
    
clf

% t=tiledlayout(length(D_d), length(v)+1,'TileSpacing','compact');
t=tiledlayout(length(v)+1, length(D_d) ,'TileSpacing','compact');
 txt = ['Light panel dim: ' num2str(L_ro) 'm,     h=' num2str(h) 'm, N=' num2str(N) ' rays,    total duration: ' num2str(nan) 'hrs'];
title(t,["Effective Reflecting Area vs Photodiode area (Monte Carlo)" txt]);
for d = 1:length(D_d)
     L_ro = R_d+max(D_d); %outer radius of light
    L_ri = D_d(d)/2; %inner radios of the light srouce (cuts out PD)
    
    fn =['monte-carlo-alpha_N' num2str(N) '_PD' num2str(D_d(d)*1000) 'mm_L' num2str(50) 'mm(' num2str(suffix) ').mat'];
    load(fn, 'detected_this_trial');
    detected0=detected_this_trial;
%     
%     fn =['monte-carlo-alpha_N' num2str(N) '_PD' num2str(D_d(d)*1000) 'mm(' num2str(2) ').mat'];
%     load(fn, 'detected_this_trial');
%     detected1=detected_this_trial;
%     
%     fn =['monte-carlo-alpha_N' num2str(2500000) '_PD' num2str(D_d(d)*1000) 'mm(' num2str(1) ').mat'];
%     load(fn, 'detected_this_trial');
%     detected2=detected_this_trial;
%     
%     detected=cell(1,length(v));
%     for trial = 1:length(v)
%         detected{trial} = vertcat(detected0{trial},detected1{trial},detected2{trial});
%     end
    
    
    for trial = 1:length(v)
       nexttile(d+(trial-1)*length(D_d))
       hold on
        viscircles(D_pos(1:2)',R_d/2,'Color','k', 'LineWidth',.1 );
        viscircles(D_pos(1:2)',D_d(d)/2,'Color','m', 'LineWidth',.025 );
        viscircles(D_pos(1:2)',L_ro,'Color','m', 'LineWidth',.025 );
        if ~isempty(detected{trial})
            scatter(detected{trial}(:,1), detected{trial}(:,2),'b', '.'); %source points that hit det
            scatter(detected{trial}(:,4), detected{trial}(:,5), 'r', '.'); %hits in detector
        end 
        xlim([-L_ro L_ro])
        ylim([-L_ro L_ro])
        axis equal
        axis square
        hold off
    end
    hold off
    nexttile(d+trial*length(D_d))

     hold on
    % x=0:0.025:.4;
    % Alpha_theory = calcEffectiveReflectingArea(x,h,R_d,R_ri,R_L,R_Ls);
    % num_detected= zeros(1,length(v));
    % for trial=1:length(v)
    %     num_detected(trial) = length(detected{trial}(:,1));
    % end

    alpha = calcAlphaForPDWithArea(D_d(d),v);
    detector_hits=zeros(1,length(v));
    U=zeros(1,length(v));
    for u = 1:length(v)
        if ~isempty(detected{u})
         detector_hits(u) = length(detected{u}(:,1));
         
         Utemp = unique(detected{u}(:,1:3), 'rows');
         U(u) = length(Utemp(:,1));
        end
    end
    plot(v,alpha,'.')
    plot(v,U/U(1), '.')
    plot(v,U/U(1)-alpha, '-')

       
    all = detector_hits/detector_hits(1);
%     plot(v,all);
    grid on
    lgd=legend(['calculated (' num2str(D_d(d)*1000) 'mm PD)'], ['simulated (' num2str(D_d(d)*1000) 'mm PD)'], "error");%, 'total hits (include redundant src pts)')
    set(lgd,'location','south outside')
    
    title(['\alpha for ' num2str(D_d(d)*1000) 'mm PD'])
hold off
end
%%
figure
colors = [  0 0.4470 0.7410;
            0.8500 0.325 0.0980;
            0.9290 0.6940 0.1250;	
            0.4940 0.184, 0.5560;
            0.4660 0.6740 0.188  ;        	
            0.30 0.7450 0.9330;	          	
            0.6350 0.0780 0.1840;    
            0 0 0;
            0 0 1; 
            0 1 0; 
            1 0 0
                ];
ax.ColorOrder = colors(1:length(D_d),:);
t2=tiledlayout('flow','TileSpacing','normal');
nexttile
hold on
plot(v,alpha)
title('Theoretical');
legend([num2str(D_d'*1000)  repmat(' mm',length(D_d),1)])
hold off
nexttile
hold on
for d=1:length(D_d)
%     fn =['grid-alpha_' num2str(L_Nx) 'x' num2str(L_Ny) 'x' num2str(R_N) '_PD' num2str(D_d(d)*1000) 'mm(' num2str(suffix) ').mat'];
%     fn =['monte-carlo-alpha_N' num2str(N) '_PD' num2str(D_d(d)*1000) 'mm(' num2str(suffix) ').mat'];
      fn =['monte-carlo-alpha_N' num2str(N) '_PD' num2str(D_d(d)*1000) 'mm_L' num2str(50) 'mm(' num2str(suffix) ').mat'];

%     load(fn, 'detected');
        load(fn, 'detected_this_trial');
   U=zeros(1,length(v));
    for u = 1:length(v)
        if ~isempty(detected{u})
         Utemp = unique(detected{u}(:,1:3), 'rows');
         U(u) = length(Utemp(:,1));
        end
    end
    plot(v,U/U(1),'.--','MarkerSize', 15);
    title('Simulated');
    legend([num2str(D_d'*1000)  repmat(' mm',length(D_d),1)])

    
end
legend
hold off

% calcAlphaForPDWithArea(D_d);
