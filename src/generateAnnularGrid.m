%% generate light ray grid on annulus
clc, clear all; tic
h=1.5;
N = 1000; % rays per side
L_s = .25;
n = -L_s/2 : L_s/(N-1) : L_s/2;
D_d = 0.001;     % detector diameter (1mm to simulate point source)

[Lx,Ly] = meshgrid(n,n);
Lx=Lx(:);
Ly=Ly(:);
F = sqrt(Lx.^2+Ly.^2) <= D_d/2; %is inside PD
Lx(F) = [];
Ly(F) = [];
sources = [Lx, Ly, h*ones(length(Lx),1)];
toc
scatter(Lx,Ly, 'y' ,'.');
axis('equal');
hold on

N = 500; % rays per side
R_d = .050;       % front face diamter (r in the paper, 50mm)
n = -R_d/2 : R_d/(N-1) : R_d/2;
[Rx,Ry] = meshgrid(n,n);
Rx=Rx(:);
Ry=Ry(:);
F = sqrt(Rx.^2+Ry.^2) >= R_d/2; %outside retro face
Rx(F) = [];
Ry(F) = [];
dests = [Rx, Ry, 0*ones(length(Rx),1)];
scatter(Rx,Ry, 'b' ,'o');
hold off
% 
% for i = 1:N
%     for j = 1:N
%         if (X(i,j)) <= s_len
%             F(i,j) = 1;
%         else
%             F(i,j) = 0;
%         end
%     end
% % end
% surf(X,Y,F)
