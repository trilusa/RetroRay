clear all;close all;clc;

L = 35.7e-3; %length of retroreflector
Ls = 6.3e-3; %length of recession
rr = 50e-3;  %diameter of the front face of retroreflector

PDSize = [10e-3 100e-3]; %m
h = 1.5; %meter
HdistanceRange = [0:0.001:0.8]; %meter

i = 1;
for Hdistance = HdistanceRange  %radiance angle = incidence angle

    theta = atand(Hdistance/h);
    D(i) = 2*L*tand(theta);
    Ds(i) = 2*Ls*tand(theta);

    j = 1;
    for PD_r = PDSize/2
        r = rr+PD_r;
        ERA(i,j) = 2*r^2*(acos((D(i)+Ds(i))/r)-(D(i)+Ds(i))/r*sin(acos((D(i)+Ds(i))/r)))-pi*PD_r^2;
        Alpha(i,j) = ERA(i,j)/(pi*r^2-pi*PD_r^2)*100;
        j = j+1;
    end
    
i = i+1;
end

% figure
% hold on
% plot([0.1:0.01:500],ERA(1,:));
% plot([0.1:0.01:500],ERA(2,:));
% plot([0.1:0.01:500],ERA(3,:));
% plot([0.1:0.01:500],ERA(4,:));
% hold off
% legend('\theta=0','\theta=10','\theta=20','\theta=30')
% xlabel('PD radius (mm)')
% ylabel('Effective Reflecting Area (mm^2)')

figure
hold on
plot(HdistanceRange,Alpha(:,1));
plot(HdistanceRange,Alpha(:,2));
hold off
legend('10 mm','100 mm')
xlabel('Horizontal Distance (m)')
ylabel('\alpha (%)')
