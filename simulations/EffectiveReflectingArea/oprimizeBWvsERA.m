%% optimize BW and ERA as function of PD diameter
D_d = .0005:.0001:.05;  %.2cm to 2cm
v=0:0.2:0.8;
alpha = calcAlphaForPDWithArea(D_d,v);
clf
hold on
plot(D_d*1000,alpha,'--')

BW = calcPDBandwidth(D_d);

plot(D_d*1000, BW/max(BW))
title('Bndwidth and ERA (\alpha) vs Active Area Diameter')
xlabel('PD diam (mm)')
% ylabel('F_{3db} (Hz)')
lgn_txt = [repmat("\alpha (", length(v),1) + num2str(v') + repmat(" m)",length(v),1)];
 lgn_txt = [lgn_txt; "f_{-3db}"];
legend(lgn_txt);
xlim([0 50])


hold off
grid on