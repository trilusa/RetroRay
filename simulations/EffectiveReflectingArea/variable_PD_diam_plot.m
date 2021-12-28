%% Script to do a fancy plot of ERA vs PD size

% update these numbers to match those of the target data
    L_Nx = 50; %10e6;         % number of rays to cast per trial
    L_Ny = 50;
    R_N = 50;
    v = 0:0.2:.8;  % horizontal displacements to simulate
    R_test_pos = [-v; zeros(2, length(v))]; % (v,0,0) positions to test retroreflector
    h = 1.5;         % height of lamp above retrorefltor (m)

% Detector Params
    D_d = 0.002:0.002:0.020;     % detector diameter (1mm to simulate point source)
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
    dest_points = [Rx, Ry, 0*ones(length(Rx),1)];

    L_sx = 2*R_d+.02; %side length for light
    L_sy = 2*R_d+.02;
    basisX = -L_sx/2 : L_sx/(L_Nx-1) : L_sx/2;
    basisY = -L_sy/2 : L_sy/(L_Ny-1) : L_sy/2;
    [Lx,Ly] = meshgrid(basisX,basisY);
    Lx=Lx(:);
    Ly=Ly(:);
% all_detected = cell(length(D_d),1);
clf

% t=tiledlayout(length(D_d), length(v)+1,'TileSpacing','compact');
t=tiledlayout(length(v)+1, length(D_d) ,'TileSpacing','compact');

txt = ['Light panel dim: ' num2str(L_sx) 'x' num2str(L_sy) 'm' ...
    ' (' num2str(L_Nx) 'x' num2str(L_Ny) ' grid),                PD: ' ...
     num2str(R_N) 'x' num2str(R_N) 'grid,                h=' num2str(h) ',                total duration: ' num2str(59.3) 'hrs'];
title(t,["Effective Reflecting Area vs Photodiode area" txt]);
suffix=1;
for d = 1:length(D_d)
    % generate light points
    
    F = sqrt(Lx.^2+Ly.^2) <= D_d(d)/2; % remve pts inside PD
    Lx(F) = [];
    Ly(F) = [];
    source_points = [Lx, Ly, h*ones(length(Lx),1)];
    
    fn =['grid-alpha_' num2str(L_Nx) 'x' num2str(L_Ny) 'x' num2str(R_N) '_PD' num2str(D_d(d)*1000) 'mm(' num2str(suffix) ').mat'];
    load(fn, 'detected');
    
%     all_detected{d} = detected;
    for trial = 1:length(v)
        nexttile(d+(trial-1)*length(D_d))
        hold on
%         rectangle('Position',[-R_d -L_sy/2 L_sx L_sy] );
        scatter(source_points(:,1), source_points(:,2), [] ,[0.9290 0.6940 0.1250],'.'); %all source pts
        viscircles([v(trial),0],R_d/2,'Color','k', 'LineWidth',.1 );
        viscircles(D_pos(1:2)',D_d(d)/2,'Color','k', 'LineWidth',.025 );
        if ~isempty(detected{trial})
            scatter(detected{trial}(:,1), detected{trial}(:,2),'b', '.'); %source points that hit det
            scatter(detected{trial}(:,4), detected{trial}(:,5), 'r', '.'); %hits in detector
        end
        axis equal
       ylim([-L_sy/2-.005 L_sy/2+.005]);
       xlim([-L_sx/2-.005 L_sx/2+.005]);
    %         xlim([-.0001-D_d(d)/2  .0001+D_d(d)/2])
    %         ylim([-.0001-D_d(d)/2  .0001+D_d(d)/2])

        txt_title=[num2str(v(trial)) 'm (' num2str(length(detected{trial})) ' hits)'];
        if trial == 1
           txt_title=[num2str(D_d(d)*1000) 'mm PD, '  num2str(v(trial)) 'm (' num2str(length(detected{trial})) ' hits)'];
        end
    %     txt_subtitle = ['time = ' num2str(elapsed(trial)/60) ' min'];
        title([convertCharsToStrings(txt_title)]);%, convertCharsToStrings(txt_subtitle )])
%         if (trial == 1) && (d == 1)
%             lgd = legend('source pt', 'detected source pt', 'PD hit' );
%             set(lgd,'location','north outside')
%         end
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
    fn =['grid-alpha_' num2str(L_Nx) 'x' num2str(L_Ny) 'x' num2str(R_N) '_PD' num2str(D_d(d)*1000) 'mm(' num2str(suffix) ').mat'];
    load(fn, 'detected');
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
