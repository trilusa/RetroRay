clc; clear all;

N = 5e6;% load('alpha_PD10mm_H1.5m_Rd1mm_L200mm_N500000(1).mat').N; % matches the simulated data

Vsim = load('alpha_PD10mm_H1.5m_Rd1mm_L200mm_N500000(1).mat').V; % matches the simulated data
Vtheo = 0:.001:100;  %to use with theoretical calcs

colors=['k';'g';'b'];

 t=tiledlayout(3,3);

% scaled retroreflector lookup table (meters)
R_d = [ .050/50 .050/5 .050 ];    
R_L = [ .0357/50 .0357/5 .0357 ];     
R_Ls = [.0063/50 .0063/5 .0063 ];  

H = [1.5 3 4.5]; %meters
D_d = [0.1, 1, 10]*1e-3; %meters
%% Effect of PD size on alpha
tic
% Dsim=[...
% load('alpha_PD0.1mm_H1.5m_Rd10mm_L200mm_N1000000(2).mat')
% load('alpha_PD1mm_H1.5m_Rd10mm_L200mm_N1000000(2).mat')
% load('alpha_PD10mm_H1.5m_Rd10mm_L200mm_N1000000(2).mat')
% 
% ];
% % 
% % Dsim=[...
% % load('alpha_PD0.1mm_H1.5m_Rd10mm_L200mm_N500000(1).mat')
% % load('alpha_PD1mm_H1.5m_Rd10mm_L200mm_N500000(1).mat')
% % load('alpha_PD10mm_H1.5m_Rd10mm_L200mm_N500000(1) (1).mat')
% % ];
Dsim=[
% load('alpha_PD0.1mm_H1.5m_Rd50mm_L200mm_N5000000(1).mat')
%         load('alpha_PD1mm_H1.5m_Rd50mm_L101mm_N5000000(1).mat')
       load('alpha_PD10mm_H1.5m_Rd50mm_L200mm_N5000000(2).mat')];
ERA_D = zeros(length(Dsim),length(Vsim)); 

i=1;
for d=1:length(Dsim)
    for v=1:length(Vsim)
        detected = cell2mat(Dsim(d).detected(v));
       if (isempty(detected) || length(detected(:,1)) < 3)
            ERA_D(d,i) = 0;
        else
            x=detected(:,1); 
            y=detected(:,2);
            [era, ~]= measureERA(x,y,Dsim(d).Dtemp);
            ERA_D(d,i) = era;
        end
        i=i+1;
    end
    i=1;

%     numHits_D(d,v)=length(Dsim(d).detected(:,1)

end
%take alpha to be num_hits (in percent)
alpha_D_sim = (ERA_D ./ ERA_D(:,1)) *100; 

%calc alpha for sim 2 (PDdiam=10mm, V=[0:25mm:.75m], H=[1.5 3 4.5], RRdiam=50mm)
alpha_D_theory = [
%     calcAlphaForPDWithArea(D_d(1), Vtheo, H(1), R_d(3), R_L(3), R_Ls(3));
%                   calcAlphaForPDWithArea(D_d(2), Vtheo, H(1), R_d(3), R_L(3), R_Ls(3));
                  calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(3), R_L(3), R_Ls(3))] * 100;
% figure(1)
% nexttile
hold on 
for i=1:length(Dsim)
    plot(Vsim,alpha_D_sim(i,:), 'Color', colors(i), 'LineStyle' , '--', 'Marker','.');
end
for i=1:length(Dsim)
    plot(Vtheo,alpha_D_theory(i,:),"Color", colors(i),'LineStyle' , '-');
end

lgn_txt = [ repmat('PD diam = ', length(Dsim), 1) num2str(1000*D_d(:,3)') repmat(' mm',length(Dsim),1)];
legend(lgn_txt)

title("horiz disp vs alpha for differnt PD sizes")
subtitle("H=1.5m, RRdiam=50mm, N=5e6 rays")
ylabel("\alpha (%)")
xlabel("v (m)")

grid minor
xlim([0 1])
ylim([0 105])
hold off

toc
%% Effect of CC size on alpha %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rsim=[
%     load('alpha_PD10mm_H1.5m_Rd1mm_L200mm_N1000000(2).mat')
%     load('alpha_PD10mm_H1.5m_Rd10mm_L200mm_N1000000(1).mat')
%      load('alpha_PD10mm_H1.5m_Rd50mm_L200mm_N1000000(1).mat')
% 
% % ];
% Rsim=[
% load('alpha_PD10mm_H1.5m_Rd1mm_L200mm_N500000(1).mat')
% load('alpha_PD10mm_H1.5m_Rd10mm_L200mm_N500000(1).mat')
% load('alpha_PD10mm_H1.5m_Rd50mm_L200mm_N500000(1).mat')
% ];


Rsim=[
%     load('alpha_PD10mm_H1.5m_Rd1mm_L60mm_N10000000(1).mat')
    load('alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000(1).mat')
    load('alpha_PD10mm_H1.5m_Rd10mm_L60mm_N10000000(1).mat')
    load('alpha_PD10mm_H1.5m_Rd50mm_L200mm_N5000000(2).mat')
];
% numHits_R = zeros(length(Rsim),length(Vsim));
ERA_R = zeros(length(Rsim),length(Vsim));

i=1;
now
tic
for r=1:length(Rsim)
    for v=1:length(Vsim)
         detected = cell2mat(Rsim(r).detected(v));
        if (isempty(detected) || length(detected(:,1)) < 3)
%             numHits_R(r,i) = 0;
            ERA_R(r,i) = 0;
        else
%             numHits_R(r,i) = length(detected(:,1)); %number of hit;
            x=detected(:,1); 
            y=detected(:,2);
            [era, ~]= measureERA(x,y,Rsim(r).Dtemp);
            ERA_R(r,i) = era;
        end
        i=i+1;
    end
    i=1;
end
%take alpha to be num_hits (in percent)
% alpha_R_sim = (numHits_R ./ numHits_R(:,1)) *100; 
alpha_R_sim = (ERA_R ./ ERA_R(:,1)) *100; 


% alpha_R_theory = [calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(1), R_L(1), R_Ls(1));
%                   calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(2), R_L(2), R_Ls(2));
%                   calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(3), R_L(3), R_Ls(3))]*100;
% % 

alpha_R_theory = calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d, R_L, R_Ls)*100;

% figure(2)
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
% subtitle("H=1.5m, PDdiam=10mm, N=5e6 rays")
subtitle(['PDdiam=' num2str(D_d(3)*1000) 'mm, H=' num2str(H(1)) 'm, ', num2str(N,'%.1e') ' rays'])

ylabel("\alpha (%)")
xlabel("v (m)")

grid minor
xlim([0 1])
ylim([0 105])
hold off
toc
%% Effect of height on Alpha
% Load the dat for Sim2, calculate alphas, then plot
% Hsim=[
%         load('alpha_PD10mm_H1.5m_Rd1mm_L200mm_N1000000(2).mat')
%         load('alpha_PD10mm_H3m_Rd1mm_L200mm_N1000000(1).mat')
%         load('alpha_PD10mm_H4.5m_Rd1mm_L200mm_N1000000(1).mat')
% ];

Hsim=[  load('alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000(1).mat')
        load('alpha_PD10mm_H3m_Rd1mm_L12mm_N10000000(1).mat')
        load('alpha_PD10mm_H4.5m_Rd1mm_L12mm_N5000000(1).mat')
       ]
     
% 
% Hsim=[
%  load('alpha_PD10mm_H1.5m_Rd50mm_L200mm_N1000000(1).mat')
% load('alpha_PD10mm_H3m_Rd50mm_L200mm_N1000000(1).mat')
% load('alpha_PD10mm_H4.5m_Rd50mm_L200mm_N2.5000(1).mat')
% ];
rrn=1;

ERA_H = zeros(length(Hsim),length(Vsim)); 

i=1;
for h=1:length(Hsim)
    for v=1:length(Vsim)
          detected = cell2mat(Hsim(h).detected(v));
        if (isempty(detected) || length(detected(:,1)) < 3)
            ERA_H(h,i) = 0;
        else
            x=detected(:,1); 
            y=detected(:,2);
            [era, ~]= measureERA(x,y,Hsim(h).Dtemp);
            ERA_H(h,i) = era;
        end
        i=i+1;
    end
    i=1;
end
%take alpha to be num_hits (in percent)
alpha_H_sim =  (ERA_H ./ ERA_H(:,1)) *100; 

%calc alpha for sim 2 (PDdiam=10mm, V=[0:25mm:.75m], H=[1.5 3 4.5], RRdiam=10mm)
% alpha_H_theory = [calcAlphaForPDWithArea(D_d(3), Vtheo, H(1), R_d(1), R_L(1), R_Ls(1));
%                   calcAlphaForPDWithArea(D_d(3), Vtheo, H(2), R_d(1), R_L(1), R_Ls(1));
%                   calcAlphaForPDWithArea(D_d(3), Vtheo, H(3), R_d(1), R_L(1), R_Ls(1))] * 100;

alpha_H_theory = calcAlphaForPDWithArea(D_d(3), Vtheo, H, R_d(rrn), R_L(rrn), R_Ls(rrn))*100;

% figure(3)
nexttile
hold on 
for i=1:length(Hsim)
    plot(Vsim,alpha_H_sim(i,:), 'Color', colors(i), 'LineStyle' , '--', 'Marker','.');
end
for i=1:length(Hsim)
    plot(Vtheo,alpha_H_theory(i,:),"Color", colors(i),'LineStyle' , '-');
end

lgn_txt = [ repmat('h = ', length(Hsim), 1) num2str(H') repmat(' m',3,1)];
legend(lgn_txt)

title("horiz disp vs alpha for differnt heights")
subtitle(['PDdiam=' num2str(D_d(3)*1000) 'mm, RRdiam=' num2str(R_d(rrn)*1000) 'mm, ', num2str(N,'%.1e') ' rays'])
ylabel("\alpha (%)")
xlabel("v (m)")

grid minor
xlim([0 1])
ylim([0 105])
hold off
toc