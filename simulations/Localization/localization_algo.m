% calcAlphaForPDWithArea(10e-3,0:);
clc, clear all

PDsquare_spacing = .1;
H=4.5;
D_d = 10e-3;
R_d = 50e-3;
R_L = .0357;     % Length of CC portion of retroreflector (35.7mm)
R_Ls = .0063; 
NumPD=5; %num PDs

%% generate fake retrorefector position
Rpos = [.45551 -.1 0];
PDpos = [ 0 0 H;
        PDsquare_spacing  PDsquare_spacing H;
       -PDsquare_spacing  PDsquare_spacing H;
        PDsquare_spacing -PDsquare_spacing H;
       -PDsquare_spacing -PDsquare_spacing H];
alpha_meas = calcAlphaBetweenPDandRR(Rpos,PDpos,D_d, R_d ,R_L,R_Ls)+ 0.005*randn(NumPD,1)

%%
calcAlphaErr = @(s) alpha_meas-calcAlphaBetweenPDandRR(s, PDpos,D_d, R_d ,R_L,R_Ls)
s=[0 0 0];
sfit=lsqnonlin(calcAlphaErr,s)
Rpos
%%
function alpha = calcAlphaBetweenPDandRR(s,PDpos, D_d, R_d ,R_L,R_Ls)
    N=length(PDpos);
    V=zeros(N,1);
    H=zeros(N,1);
    alpha = zeros(N,1);

    for i = 1:N
        V(i) =  sqrt(sum((s(1:2)- PDpos(i,1:2)).^2)); %distance on xy-plane
        H(i) = PDpos(i,3) - s(3); %vertical distance
        % calculate predicted alphas with the model
        alpha(i) = calcAlphaForPDWithArea(D_d, V(i),H(i),R_d ,R_L ,R_Ls);
    end
end








