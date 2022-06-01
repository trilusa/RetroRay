clc; clear all;
load('final_data_large.mat');
load('alpha_PD1mm_H1.5m_Rd50mm_L101mm_N10000000(1).mat', 'V')
D_d = {[0.1,1,10]*1e-3, 10e-3, 1e-3};
R_d = {10e-3, [1 10 50]*1e-3, 50e-3};
H   = {1.5, 1.5, [1.3 3 4.5]};
params = {D_d{1}, R_d{1},H{1};
    D_d{2},R_d{2},H{2};
    D_d{3},R_d{3},H{3}};

figure(100001)
t=tiledlayout(3,3,"Padding","compact", "TileSpacing","normal");
c = ['r' 'g' 'b'];
% 'Color'order({c(1) C(2) C(3); C(1) C(2) C(3)});


lgntxt = {[repmat('PD diam = ',3, 1) num2str(1.*[.1 1 10 ]') repmat(' mm',3,1)];
    [repmat('CC diam = ', 3, 1) num2str([1 10 50]') repmat(' mm',3,1)];
    [repmat('h = ', 3, 1) num2str([1.5 3.0 4.5]') repmat(' m',3,1)]};
titletxt = {
    ['Sim ' num2str(1) ': N=10e7 rays, CCdiam=' num2str(params{1,2}*1e3) 'mm, Height='  num2str(num2str(params{1,3})) 'm']
    ['Sim ' num2str(2) ': N=10e7 rays, PDdiam=' num2str(params{2,1}*1e3) 'mm, Height='  num2str(num2str(params{2,3})) 'm']
    ['Sim ' num2str(3) ': N=10e7 rays, PDdiam=' num2str(params{3,1}*1e3) 'mm, CCdiam='  num2str(num2str(params{3,2}*1e3)) 'mm']
    };

C={ [0 0.4470 0.7410]
    [0.8500 0.3250 0.0980]
    [0.9290 0.6940 0.1250]};
transparency = .8;

i=1;
rlim=.75;
fidx = find(V==rlim);
V=V(1:fidx);
for s=1:3:9

    %% alpha
    nexttile
    yyaxis left
    hold on
%      [alpha_th, ERA_th ] = calcAlphaForPDWithArea(D(i),V,H(i),R(i),RL,RLs)
    plot(V, 100*alpha_th(s,  1:fidx),'-', 'Color',C{1})
    plot(V, 100*alpha_th(s+1,1:fidx),'-' ,'Color',C{2})
    plot(V, 100*alpha_th(s+2,1:fidx), '-','Color',C{3})
    plot(V, 100*alpha(s,1:fidx),'--.', 'Color',C{1})
    plot(V, 100*alpha(s+1,1:fidx),'--.','Color',C{2})
    plot(V, 100*alpha(s+2,1:fidx),'--.','Color',C{3})
    ylabel('\alpha (%)')
    ylim([0 105])
    xlim([0 rlim])

    yyaxis right
    plot(V,abs(100*(alpha_th(s,1:fidx)-alpha(s,1:fidx))),'--', 'Color',[C{1} transparency])
    plot(V,abs(100*(alpha_th(s+1,1:fidx)-alpha(s+1,1:fidx))),'--','Color',[C{2} transparency])
    plot(V,abs(100*(alpha_th(s+2,1:fidx)-alpha(s+2,1:fidx))),'--','Color',[C{3} transparency])
    ylabel("|\alpha error|")
%     ylim([-.1 .1])
%     title(titletxt(i))
    hold off
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    grid minor

%%      ERA
    nexttile
    yyaxis left
    hold on
    p(1)=plot(V, 1e6*ERA_th(s,  1:fidx), '-','Color',C{1});
    p(2)=plot(V, 1e6*ERA_th(s+1,1:fidx),'-', 'Color',C{2});
    p(3)=plot(V, 1e6*ERA_th(s+2,1:fidx), '-','Color',C{3});

    plot(V,  1e6*ERA(s,1:fidx),'--.', 'Color',C{1})
    plot(V,  1e6*ERA(s+1,1:fidx),'--.','Color',C{2})
    plot(V,  1e6*ERA(s+2,1:fidx),'--.','Color',C{3})
    ylabel('ERA (mm^2)')
    xlim([0 rlim])

    yyaxis right
    plot(V,abs(1e6*(ERA_th(s,1:fidx)-ERA(s,1:fidx))),'--','Color',[C{1} transparency] )
    plot(V,abs(1e6*(ERA_th(s+1,1:fidx)-ERA(s+1,1:fidx))),'--', 'Color',[C{2} transparency])
    plot(V,abs(1e6*(ERA_th(s+2,1:fidx)-ERA(s+2,1:fidx))),'--','Color',[C{3} transparency])


    ylabel("|ERA error| (mm^2)")
%     ylim([- 5e-4])
    title(titletxt(i))

    hold off
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    grid minor


    %% PD Area
    nexttile
    yyaxis left
    hold on
    p(1)=plot(V, 1e6*(Apd_th(s, 1:fidx)), '-','Color',C{1});
    p(2)= plot(V,  1e6*(Apd_th(s+1,1:fidx)),'-', 'Color',C{2});
    p(3)=plot(V,  1e6*(Apd_th(s+2,1:fidx)), '-','Color',C{3});

    plot(V,  1e6*Apd(s,1:fidx),'--.', 'Color',C{1})
    plot(V,  1e6*Apd(s+1,1:fidx),'--.','Color',C{2})
    plot(V,  1e6*Apd(s+2,1:fidx),'--.','Color',C{3})
    ylabel('A_{pd} (mm^2)')
     xlim([0 rlim])

    yyaxis right
    plot(V,abs(1e6*( (Apd_th(s,1:fidx)-Apd(s,1:fidx) ))),'--','Color',[C{1} transparency] )
    plot(V,abs(1e6*( (Apd_th(s+1,1:fidx)-Apd(s+1,1:fidx)))),'--', 'Color',[C{2} transparency])
    plot(V,abs(1e6*( (Apd_th(s+2,1:fidx)-Apd(s+2,1:fidx)))),'--','Color',[C{3} transparency])

    ylabel("|A_{pd} Error| (mm^2)")
    legend(p(:),lgntxt{i}, 'Location','northeast', 'Orientation','Vertical' );
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    grid minor

    
    i=i+1;
end
title(t,'Theory vs Simulated for \alpha, ERA, and A_{PD}')
xlabel(t,'v (m)')


