sim=[
    load('alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000(1).mat')
load('alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000(3).mat')
load('alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000(4).mat')
load('alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000(5).mat')
load('alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000(6).mat')
load('alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000(8).mat')]


alldet=cell(1,41);
for s=1:6
    for v=1:41 
    alldet{v} = [alldet{v};  sim(s).detected{v}];
    end
end
 N= sim(s).N;
    V= sim(s).V;
    L_do = sim(s).L_do;
    Rtemp= sim(s).Rtemp;
    Dtemp= sim(s).Dtemp;
    Htemp= sim(s).Htemp;
detected=alldet;
save('0.1PD_alldet.mat',"detected",'N','V',"Rtemp","Htemp","Dtemp","L_do")