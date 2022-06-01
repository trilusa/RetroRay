
%%
tic
simdata={{
%     case1
load('01-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('02-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('03-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('04-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('05-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('06-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('07-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('08-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('09-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('10-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('11-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('12-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('13-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('14-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('15-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('16-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('17-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('18-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N10000000.mat')
load('19-alpha_PD0.1mm_H1.5m_Rd10mm_L20.1mm_N100000000.mat')
load('20-alpha_PD0.1mm_H1.5m_Rd10mm_L200mm_N1000000.mat')
load('21-alpha_PD0.1mm_H1.5m_Rd10mm_L200mm_N1000000.mat')
load('monte-carlo-alpha_N1000000_PD0.1mm_L10.1mm_h1.5m_Rd10mm(1).mat')
load('monte-carlo-alpha_N1000000_PD0.1mm_L10.1mm_h1.5m_Rd10mm(2).mat')
load('monte-carlo-alpha_N1000000_PD0.1mm_L10.1mm_h1.5m_Rd10mm(3).mat')};

% case2
{load('00-alpha_PD1mm_H1.5m_Rd10mm_L21mm_N100000000.mat')
load('01-alpha_PD1mm_H1.5m_Rd10mm_L200mm_N1000000.mat')
load('02-alpha_PD1mm_H1.5m_Rd10mm_L200mm_N1000000.mat')
load('03-alpha_PD1mm_H1.5m_Rd10mm_L21mm_N10000000.mat')
load('04-alpha_PD1mm_H1.5m_Rd10mm_L21mm_N10000000.mat')
load('05-alpha_PD1mm_H1.5m_Rd10mm_L21mm_N10000000.mat')
load('monte-carlo-alpha_N1000000_PD1mm_L11mm_h1.5m_Rd10mm(1).mat')
load('monte-carlo-alpha_N1000000_PD1mm_L11mm_h1.5m_Rd10mm(2).mat')
load('monte-carlo-alpha_N1000000_PD1mm_L11mm_h1.5m_Rd10mm(3).mat')};
% case 3
{ load('00-alpha_PD10mm_H1.5m_Rd10mm_L30mm_N15000000.mat')
load('00-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('01-alpha_PD10mm_H1.5m_Rd10mm_L200mm_N1000000.mat')
load('01-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm(.mat')
load('02-alpha_PD10mm_H1.5m_Rd10mm_L200mm_N1000000.mat')
load('02-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm(.mat')
load('03-alpha_PD10mm_H1.5m_Rd10mm_L200mm_N500000.mat')
load('03-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('04-alpha_PD10mm_H1.5m_Rd10mm_L200mm_N5000000.mat')
load('04-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('05-alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000.mat')
load('05-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('06-alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000.mat')
load('06-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('07-alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000.mat')
load('07-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('08-alpha_PD10mm_H1.5m_Rd10mm_L60mm_N10000000.mat')
load('08-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('09-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('10-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')};
% case 4
{load('00-alpha_PD10mm_H1.5m_Rd1mm_L12mm_N15000000.mat')
load('01-alpha_PD10mm_H1.5m_Rd1mm_L12mm_N25000000.mat')
load('02-alpha_PD10mm_H1.5m_Rd1mm_L12mm_N1000000.mat')
load('03-alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000.mat')
load('04-alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000.mat')
load('05-alpha_PD10mm_H1.5m_Rd1mm_L12mm_N10000000.mat')
load('06-alpha_PD10mm_H1.5m_Rd1mm_L200mm_N1000000.mat')
load('07-alpha_PD10mm_H1.5m_Rd1mm_L200mm_N1000000.mat')
load('08-alpha_PD10mm_H1.5m_Rd1mm_L200mm_N500000.mat')
load('09-alpha_PD10mm_H1.5m_Rd1mm_L200mm_N500000.mat')
load('10-alpha_PD10mm_H1.5m_Rd1mm_L200mm_N5000000.mat')
load('11-alpha_PD10mm_H1.5m_Rd1mm_L200mm_N5000000.mat')
load('12-alpha_PD10mm_H1.5m_Rd1mm_L60mm_N10000000.mat')
load('monte-carlo-alpha_N1000000_PD10mm_L11mm_h1.5m_Rd1mm(1).mat')
load('monte-carlo-alpha_N1000000_PD10mm_L11mm_h1.5m_Rd1mm(2).mat')
load('monte-carlo-alpha_N1000000_PD10mm_L11mm_h1.5m_Rd1mm(3).mat')
load('monte-carlo-alpha_N1000000_PD10mm_L11mm_h1.5m_Rd1mm(4).mat')};
% 5

{load('00-alpha_PD10mm_H1.5m_Rd10mm_L30mm_N15000000.mat')
load('00-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('01-alpha_PD10mm_H1.5m_Rd10mm_L200mm_N1000000.mat')
load('01-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm(.mat')
load('02-alpha_PD10mm_H1.5m_Rd10mm_L200mm_N1000000.mat')
load('02-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm(.mat')
load('03-alpha_PD10mm_H1.5m_Rd10mm_L200mm_N500000.mat')
load('03-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('04-alpha_PD10mm_H1.5m_Rd10mm_L200mm_N5000000.mat')
load('04-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('05-alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000.mat')
load('05-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('06-alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000.mat')
load('06-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('07-alpha_PD10mm_H1.5m_Rd10mm_L30mm_N10000000.mat')
load('07-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('08-alpha_PD10mm_H1.5m_Rd10mm_L60mm_N10000000.mat')
load('08-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('09-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')
load('10-monte-carlo-alpha_N1000000_PD10mm_L20mm_h1.5m_Rd10mm.mat')};
% 6
{load('00-alpha_PD10mm_H1.5m_Rd50mm_L110mm_N10000000.mat')
load('01-alpha_PD10mm_H1.5m_Rd50mm_L110mm_N10000000.mat')
load('02-alpha_PD10mm_H1.5m_Rd50mm_L110mm_N5000000.mat')
load('03-alpha_PD10mm_H1.5m_Rd50mm_L200mm_N1000000.mat')
load('04-alpha_PD10mm_H1.5m_Rd50mm_L200mm_N1000000.mat')
load('05-alpha_PD10mm_H1.5m_Rd50mm_L200mm_N500000.mat')
load('06-alpha_PD10mm_H1.5m_Rd50mm_L200mm_N500000.mat')
load('07-alpha_PD10mm_H1.5m_Rd50mm_L200mm_N5000000.mat')
load('08-alpha_PD10mm_H1.5m_Rd50mm_L200mm_N5000000.mat')
load('monte-carlo-alpha_N1000000_PD10mm_L60mm_h1.5m_Rd50mm(1).mat')
load('monte-carlo-alpha_N1000000_PD10mm_L60mm_h1.5m_Rd50mm(2).mat')
load('monte-carlo-alpha_N1000000_PD10mm_L60mm_h1.5m_Rd50mm(3).mat')};
% 7
{ load('00-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N10000000.mat')
load('01-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N10000000.mat')
load('02-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N100000000.mat')
load('03-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N5000000.mat')
load('04-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N5000000.mat')
load('05-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N5000000.mat')
load('06-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N5000000.mat')
load('07-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N5000000.mat')
load('08-alpha_PD1mm_H1.5m_Rd50mm_L101mm_N5000000.mat')
load('09-alpha_PD1mm_H1.5m_Rd50mm_L200mm_N1000000.mat')
load('10-alpha_PD1mm_H1.5m_Rd50mm_L200mm_N1000000.mat')
load('11-alpha_PD1mm_H1.5m_Rd50mm_L200mm_N500000.mat')
load('12-alpha_PD1mm_H1.5m_Rd50mm_L200mm_N5000000.mat')};
% 8
{load('00-alpha_PD1mm_H3m_Rd50mm_L101mm_N10000000.mat')
load('01-alpha_PD1mm_H3m_Rd50mm_L101mm_N10000000.mat')
load('02-alpha_PD1mm_H3m_Rd50mm_L101mm_N100000000.mat')
load('03-alpha_PD1mm_H3m_Rd50mm_L200mm_N1000000.mat')};
% 9
{ load('00-alpha_PD1mm_H4.5m_Rd50mm_L101mm_N10000000.mat')
load('01-alpha_PD1mm_H4.5m_Rd50mm_L101mm_N10000000.mat')
load('02-alpha_PD1mm_H4.5m_Rd50mm_L101mm_N100000000.mat')
load('03-alpha_PD1mm_H4.5m_Rd50mm_L200mm_N1000000.mat')}};
toc
%%
%remove the h
data=cell(1,9);

clc
for i = 1:length(simdata)
    tic
    for j=length(simdata{i})
        if(isfield(simdata{i}{j},'detected_this_trial'))
            simdata{i}{j}.detected=simdata{i}{j}.detected_this_trial;
        end
        for k = 1:length(simdata{i}{j}.detected)
            if(~isempty(simdata{i}{j}.detected{k}))
                srcXY=simdata{i}{j}.detected{k}(:,1:2);
                destXY=simdata{i}{j}.detected{k}(:,4:5);
                data{i}=[data{i} ; srcXY destXY-1];
            end
        end
    end
    toc
end