%  clear all; clc; sim=[];
 sim=[];
 %%
 %Comment out any line to not plot it -------------------------
% Sim 1
% sim = [sim; load('alpha_PD1mm_H1.5m_Rd50mm_L101mm_N5000000(1).mat')];
%  sim = [sim; load('alpha_PD10mm_H1.5m_Rd50mm_L200mm_N5000000(2).mat')];
%sim 2
% sim=[sim; load('alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000(1).mat')];
% sim=[sim; load('alpha_PD1mm_H1.5m_Rd10mm_L21mm_N10000000(1).mat')];
% sim = [sim; load('alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000(1).mat')];
% sim = [sim; load('alpha_PD10mm_H1.5m_Rd10mm_L60mm_N10000000(1).mat')];
% %sim 3
% sim = [sim; load('alpha_PD10mm_H3m_Rd1mm_L200mm_N5000000(1).mat')];
% sim = [sim; load('alpha_PD10mm_H4.5m_Rd1mm_L200mm_N5000000(1).mat')];

% sim = [sim; load('alpha_PD10mm_H1.5m_Rd50mm_L110mm_N50000(1).mat')];
% sim = [sim; load('alpha_PD10mm_H1.5m_Rd50mm_L110mm_N20000(1).mat')];
% sim = [sim; load('alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000(6).mat')];

% sim = [sim;       load('0.1PD_alldet.mat')];
% sim = [sim;       load('0.1PD_alldet.mat')];
% sim = [sim; load('alpha_PD1mm_H1.5m_Rd10mm_L21mm_N10000000(1).mat')];
% sim = [sim; load('alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000(2).mat')];

%  sim = [sim; load('alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000(2).mat')];
%  sim = [sim; load('alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000(2).mat')];
 sim = [sim; load('alpha_PD10mmfi_H1.5m_Rd50mm_L110mm_N10000000(1).mat')];

% sim = [sim;load('alpha_PD1mm_H1.5m_Rd50mm_L101mm_N10000000(1).mat')];
% sim = [sim;load('alpha_PD1mm_H3m_Rd50mm_L101mm_N10000000(1).mat')];
% sim = [sim;load('alpha_PD1mm_H4.5m_Rd50mm_L101mm_N10000000(1).mat')];

% sim = [sim; load('alpha_PD1mm_H3m_Rd50mm_L101mm_N10000(1).mat')];
% sim = [sim; load('alpha_PD1mm_H4.5m_Rd50mm_L101mm_N10000(1).mat')];
%%
%----------------------------------------------------------------------
V= sim(1).V;
%
% figure(123123)

t=tiledlayout('flow','TileSpacing','tight','Padding','tight');
% t=tiledlayout(length(sim), 8,'TileSpacing','loose','Padding','none');
% t=tiledlayout(2,2,'TileSpacing','none','Padding','compact');
alpha = zeros(length(sim),length(V));
ERA = zeros(length(sim),length(V));
ERAth = zeros(length(sim),length(V));
alphath = zeros(length(sim),length(V));

for s=1:length(sim) 
    disp(s);
    detected=sim(s).detected;
    N= sim(s).N;
    
    R= sim(s).Rtemp;
    D= sim(s).Dtemp;
    H= sim(s).Htemp;
    L_ro = sim(s).L_do/2;  
    idx=zeros(length(V),1);
    idx(30:5:length(V)-1)=1;
%     V = V(1:end-1);
    cnt=1;
    detv1 = unique(detected{1}, 'rows');
    %simulated era
    [ERA0, ~] = measureERA(detv1(:,1), detv1(:,2) , D);  
    for v = 1:5:length(V)-1
        detv=unique(detected{v}, 'rows');
        tic
       nexttile
       hold on
        viscircles([V(v),0], R,'Color','k', 'LineWidth',.25 );
%         viscircles([0, 0], D/2,'Color','b', 'LineWidth',.25 );
%         viscircles([0, 0], L_ro,'Color','m', 'LineWidth',.5 );
        if ~isempty(detv)
            id=randperm(length(detv(:,1)),floor(length(detv(:,1))/1));
            detv=detv(id,:);
            xto=detv(:,1);
            yto=detv(:,2);
            xti=detv(:,4);
            yti=detv(:,5);
    
            xo=xto(((xto.^2+yto.^2)<3*R));
            yo=yto((xto.^2+yto.^2)<3*R);
            xi=xti(((xto.^2+yto.^2)<3*R));
            yi=yti((xto.^2+yto.^2)<3*R);

            if ~(length(xo) < 3)
            scatter(xo*1e3, yo*1e3,10, [0,.4,1],'o','filled','MarkerFaceAlpha', .25); %source points that hit det
            scatter(xi*1e3, yi*1e3,10, [.75,0,0],'o','filled','MarkerFaceAlpha', 1); %hits in detector
            [era, xhull,yhull] = measureERA(xo, yo, D);
            ERA(s,v)=era;
            plot(xhull*1e3,yhull*1e3,'k',LineWidth=1);
            alpha(s,v)=era/ERA0;
            end
        end         

        side=((D/2 + R)*2)*1e3;
        xlim([-side/2 side/2])
        ylim([-side/2 side/2])
        
        if(cnt==7)
            nexttile(s*7-4)
%              title(...%side*.25,(side/2)+(side*.05), ...
%                     ['Dia_{cc}= ' num2str(R*1000) ...
%                     'mm, Dia_{PD}= ' num2str(D*1000) ...
%                     'mm, H = ' num2str(H) ...
%                     'm,  Dia_{L}= ' num2str(sim(s).L_do*1000) ...
%                     'mm, N = ' num2str(N,'%.1e') ' rays' ],...
%                     'FontSize', 11 ...
%                     );
            nexttile(s*cnt)
         end   
        if(cnt==1)
            ylabel('mm')
        end
%         if(cnt~=1)
            xticklabels([])
            yticklabels([])
%         end
        xlabel(['x = ' num2str(V(v)) 'm']);%, \alpha=' num2str(alpha)])
        axis square
%         grid minor
        cnt=cnt+1;
        toc
    end
end
%%
% 
% 
% t=tiledlayout(3,3);
% theory = cell(3,1);
% for s=1:3
%     L = .0357*(params{s,Ridx}/50e-3);
%     Ls = .0063*(params{s,Ridx}/50e-3);
% %     measureERA(x,y,sim_params{s,PDidx}
%     [alpha_th, era, Apd] = calcAlphaForPDWithArea(params{s,PDidx},V,params{s,Hidx},params{s,Ridx},L,Ls);
%     theory{s} = {alpha_th*1e2, era, Apd};
% 
%     
%     nexttile
%     plot(V,alpha_th)
%     plot(V, alpha);
%     ylabel('\alpha (%)')
%     grid minor
%     axis square
%     title(['Sim ' num2str(s) ':'])
%    
%     nexttile
%     plot(V,era)
%     ylabel('ERA (mm^2)')
%     grid minor
%     axis square
% 
%     nexttile
%     plot(V,Apd);
%     ylabel('A_{PD} (mm^2)')
%     legend(lgntxt{s},"Location", 'eastoutside');
%     grid minor
%     axis square
% end
% nexttile(2)
%  title("Theory vs Simulated for \alpha, ERA, and A_{PD}")
% xlabel(t,'v (m)')
