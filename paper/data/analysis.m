clc; clear all;
%%

% data = [load('0.1PD_alldet.mat'),...
%         load('alpha_PD1mm_H1.5m_Rd10mm_L21mm_N10000000(1).mat'),...
%         load('alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000(2).mat'); %replace this one with combunation
%        
%         load('alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000(2).mat'),...
%         load('alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000(2).mat'),...
%         load('alpha_PD10mm_H1.5m_Rd50mm_L110mm_N10000000(1).mat');
% 
%         load('alpha_PD1mm_H1.5m_Rd50mm_L101mm_N10000000(1).mat'),...
%         load('alpha_PD1mm_H3m_Rd50mm_L101mm_N10000000(1).mat'),...
%         load('alpha_PD1mm_H4.5m_Rd50mm_L101mm_N10000000(1).mat')];

data = [load('0.1PD_alldet.mat'),...
        load('alpha_PD1mm_H1.5m_Rd10mm_L21mm_N10000000(1).mat'),...
        load('alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000(2).mat'); %replace this one with combunation
       
        load('alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000(2).mat'),...
        load('alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000(2).mat'),...
        load('alpha_PD10mm_H1.5m_Rd50mm_L110mm_N10000000(1).mat');

        load('alpha_PD1mm_H1.5m_Rd50mm_L101mm_N10000000(1).mat'),...
        load('alpha_PD1mm_H3m_Rd50mm_L101mm_N10000000(1).mat'),...
        load('alpha_PD1mm_H4.5m_Rd50mm_L101mm_N10000000(1).mat')];
data=data';
sim=data(:);
%%
alpha=ones(length(sim),length(sim(1).V));
ERA=ones(length(sim),length(sim(1).V));
Apd=ones(length(sim),length(sim(1).V));

alpha_th=ones(length(sim),length(sim(1).V));
ERA_th=ones(length(sim),length(sim(1).V));
Apd_th=ones(length(sim),length(sim(1).V));

tic
for s=1:length(sim) 
    detected=sim(s).detected;
    N= sim(s).N;
    V= sim(s).V;
    R= sim(s).Rtemp;
    RL = .0357*(R/50e-3);
    RLs = .0063*(R/50e-3);

    D= sim(s).Dtemp;
    H= sim(s).Htemp;
    L_do=sim(s).L_do;
    L_ro = L_do/2;  
    idx=zeros(length(V),1);
    idx(1:5:length(V)-1)=1;
    cnt=1;
    detv1 = unique(detected{1}, 'rows');
    %simulated era
    [ERA0, ~] = measureERA(detv1(:,1), detv1(:,2) , D);
    for v = 1:length(sim(1).V)
        disp(num2str([s v toc]));
        detv=unique(detected{v},'rows');
        if ~isempty(detv)
            xt=detv(:,1);
            yt=detv(:,2);
            x=xt(((xt.^2+yt.^2)<3*R));
            y=yt((xt.^2+yt.^2)<3*R);
            if ~(length(x) < 3)
                [ERA(s,v), xhull,yhull] = measureERA(x, y, D);
                alpha(s,v)=ERA(s,v)/ERA0;
            end
            Apd(s,v) = measureApd(detv(:,4),detv(:,5));
        end
    end
    %%
   % theoretical era 
   [atemp, eratemp, APDtemp] = calcAlphaForPDWithArea(D,V,H,R,RL,RLs);
   alpha_th(s,:) = atemp;
   ERA_th(s,:) = eratemp;
   Apd_th(s,:) = real(APDtemp);
%    save('processed_alpha_era_theory_sim(3).mat',"Apd_th",'-append')

end
toc
 save('final_data_small.mat','alpha','alpha_th',"ERA","ERA_th","Apd" ,"Apd_th");