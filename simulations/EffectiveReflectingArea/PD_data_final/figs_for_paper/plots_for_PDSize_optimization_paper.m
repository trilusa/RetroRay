clc; clear all;

Vsim = 0:0.025:.75; % matches the simulated data
Vtheo = 0:.01:3.5;  %to use with theoretical calcs

colors=['r';'g';'b'];
t=tiledlayout(1,3);

% scaled retroreflector lookup table (meters)
R_d = [ .050/50 .050/5 .050 ];    
R_L = [ .0357/50 .0357/5 .0357 ];     
R_Ls = [.0063/50 .0063/5 .0063 ];  

H = [1.5 3 4.5]; %meters
D_d = [0.1, 1, 10]*1e-3; %meters
N = 1e6 * 2;

Rsim=[];
Rdata=[];

DSim=[];
Ddata=[];

%% Effect of height on Alpha
% Load the dat for Sim2, calculate alphas, then plot
Hsim=[
    load('monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm(1).mat');
    load('monte-carlo-alpha_N1000000_PD10mm_L20mm_h3m_Rd10mm(1).mat');
    load('monte-carlo-alpha_N1000000_PD10mm_L20mm_h4.5m_Rd10mm(1).mat')
    ];
numHits_H = zeros(length(Hsim),length(Vsim)); 

i=1;
for h=1:length(Hsim)
    for v=1:2:62
        series1 = Hsim(h).detected_this_trial(v);
        series2 = Hsim(h).detected_this_trial(v+1); %doing each horiz disp in two parts,    
        series1 = cell2mat(series1);
        series2 = cell2mat(series2);
        all_rays = [series1; series2]; % Nx6 where cols 1:3 is ray start pos, and 3:6 is ray end pos
        
        numHits_H(h,i) = length(all_rays(:,1)); %number of hit;
        i=i+1;
    end
    i=1;
end
%take alpha to be num_hits (in percent)
alpha_H_sim = (numHits_H ./ numHits_H(:,1)) * 100; 

%calc alpha for sim 2 (PDdiam=10mm, V=[0:25mm:.75m], H=[1.5 3 4.5], RRdiam=10mm)
rrn=3;
alpha_H_theory = [calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(rrn), R_L(rrn), R_Ls(rrn));
                  calcAlphaForPDWithArea(D_d(3), Vtheo, H(2), R_d(rrn), R_L(rrn), R_Ls(rrn));
                  calcAlphaForPDWithArea(D_d(3), Vtheo, H(3), R_d(rrn), R_L(rrn), R_Ls(rrn))] * 100;

% figure(1)
nexttile
hold on 
for i=1:length(Hsim)
    plot(Vsim,alpha_H_sim(i,:), 'Color', colors(i), 'LineStyle' , '--', 'Marker','.');
end
for i=1:length(Hsim)
    plot(Vtheo,alpha_H_theory(i,:),"Color", colors(i),'LineStyle' , '-');
end

lgn_txt = [ repmat('h = ', length(Hsim), 1) num2str([Hsim.h]') repmat(' m',3,1)];
legend(lgn_txt)

title("horiz disp vs alpha for differnt heights")
subtitle("PDdiam=10mm, RRdiam=10mm, N=2e6 rays")
ylabel("\alpha (%)")
xlabel("v (m)")

grid minor
% xlim([0 3.5])
% ylim([0 110])
hold off

%% Effect of PD size on alpha

Dsim=[
    load('monte-carlo-alpha_N1000000_PD0.1mm_L10.1mm_h1.5m_Rd10mm(1).mat'); %.1mm PD
    load('monte-carlo-alpha_N1000000_PD1mm_L11mm_h1.5m_Rd10mm(1).mat'); %1mm
    load('monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm(9).mat'); %10mm
    ];
numHits_D = zeros(length(Dsim),length(Vsim)); 

i=1;
for d=1:length(Dsim)
    for v=1:2:62
        series1 = Dsim(d).detected_this_trial(v);
        series2 = Dsim(d).detected_this_trial(v+1); %doing each horiz disp in two parts,    
        series1 = cell2mat(series1);
        series2 = cell2mat(series2);
        all_rays = [series1; series2]; % Nx6 where cols 1:3 is ray start pos, and 3:6 is ray end pos
        numHits_D(d,i) = length(all_rays(:,1)); %number of hit;
        i=i+1;
    end
    i=1;
end
%take alpha to be num_hits (in percent)
alpha_D_sim = (numHits_D ./ numHits_D(:,1)) * 100; 

%calc alpha for sim 2 (PDdiam=10mm, V=[0:25mm:.75m], H=[1.5 3 4.5], RRdiam=10mm)
alpha_D_theory = [calcAlphaForPDWithArea(D_d(1), Vtheo, H(1), R_d(2), R_L(2), R_Ls(2));
                  calcAlphaForPDWithArea(D_d(2), Vtheo, H(1), R_d(2), R_L(2), R_Ls(2));
                  calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(2), R_L(2), R_Ls(2))] * 100;

% figure(2)
nexttile
hold on 
for i=1:length(Dsim)
    plot(Vsim,alpha_D_sim(i,:), 'Color', colors(i), 'LineStyle' , '--', 'Marker','.');
end
for i=1:length(Dsim)
    plot(Vtheo,alpha_D_theory(i,:),"Color", colors(i),'LineStyle' , '-');
end

lgn_txt = [ repmat('PD diam = ', length(Dsim), 1) num2str(1000*D_d') repmat(' mm',3,1)];
legend(lgn_txt)

title("horiz disp vs alpha for differnt PD sizes")
subtitle("H=1.5m, RRdiam=10mm, N=2e6 rays")
ylabel("\alpha (%)")
xlabel("v (m)")

grid minor
% xlim([0 3.5])
% ylim([0 110])
hold off


%% Effect of CC size on alpha %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsim=[
    load('monte-carlo-alpha_N1000000_PD10mm_L11mm_h1.5m_Rd1mm(1).mat') %RR=1mm
    load('monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm(6).mat') %10mm
    load('monte-carlo-alpha_N1000000_PD10mm_L60mm_h1.5m_Rd50mm(1).mat') %50mm
    ];
numHits_R = zeros(length(Rsim),length(Vsim)); 

i=1;
for r=1:length(Rsim)
    for v=1:2:62
        series1 = Rsim(r).detected_this_trial(v);
        series2 = Rsim(r).detected_this_trial(v+1); %doing each horiz disp in two parts,    
        series1 = cell2mat(series1);
        series2 = cell2mat(series2);
        all_rays = [series1; series2]; % Nx6 where cols 1:3 is ray start pos, and 3:6 is ray end pos
        numHits_R(r,i) = length(all_rays(:,1)); %number of hit;
        i=i+1;
    end
    i=1;
end
%take alpha to be num_hits (in percent)
alpha_R_sim = (numHits_R ./ numHits_R(:,1)) * 100; 

%calc alpha for sim 2 (PDdiam=10mm, V=[0:25mm:.75m], H=[1.5 3 4.5], RRdiam=10mm)
alpha_R_theory = [calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(1), R_L(1), R_Ls(1));
                  calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(2), R_L(2), R_Ls(2));
                  calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(3), R_L(3), R_Ls(3))] * 100;

% figure(3)
nexttile
hold on 
for i=1:length(Rsim)
    plot(Vsim,alpha_R_sim(i,:), 'Color', colors(i), 'LineStyle' , '--', 'Marker','.');
end
for i=1:length(Rsim)
    plot(Vtheo,alpha_R_theory(i,:),"Color", colors(i),'LineStyle' , '-');
end

lgn_txt = [ repmat('CC face diam = ', length(Rsim), 1) num2str(1000*R_d') repmat(' mm',3,1)];
legend(lgn_txt)

title("horiz disp vs alpha for differnt CC diams")
subtitle("H=1.5m, PDdiam=10mm, N=2e6 rays")
ylabel("\alpha (%)")
xlabel("v (m)")

grid minor
% xlim([0 3.5])
% ylim([0 110])
hold off
