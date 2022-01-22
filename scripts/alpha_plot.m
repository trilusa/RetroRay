D_d = (0.1:0.01:200)*1e-3
plot(D_d*1000,calcAlphaForPDWithArea(D_d,v));
grid minor
xlabel 'PD diameter (mm)'
ylabel '\alpha'
title '\alpha vs PD diameter (.1 mm  to 4mm)'
ylim( [0 1.01])
legend({'v = 0m', 'v = .2m', 'v = .4m', 'v = .6m', 'v = .8m', })

% f3db = calcPDBandwidth(D_d);
% plot(D_d*1000,f3db)
% grid minor
% xlabel 'PD diameter (mm)'
% ylabel 'BW (freq_{3db})'
% title 'BW vs PD diameter (.1 mm  to 4mm)'
