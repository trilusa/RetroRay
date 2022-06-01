N=10e4;
Xl = zeros(1,N);
Yl = zeros(1,N);
Xr = zeros(1,N);
Yr = zeros(1,N);


D_d=10e-3;
R_d=10e-3;
L_do = D_d+2*R_d;
L_ro = L_do/2;
L_ri = D_d/2;

for i = 1:N
%             generate the source pt on light ring       
            x = (rand(1,1)-0.5)*L_do;
            y = (rand(1,1)-0.5)*L_do;
            while( (x^2+y^2 > L_ro^2) ||  (x^2+y^2 < L_ri^2))
                x = (rand(1,1)-0.5)*L_do;
                y = (rand(1,1)-0.5)*L_do;
            end
            Xl(i)=x;
            Yl(i)=y;


        % generate source point on the detector (reverse ray tracing)
%             x = (rand(1,1)-0.5)*D_d(d);
%             y = (rand(1,1)-0.5)*D_d(d);
%             while( (x^2+y^2) > (D_d(d)/2)^2 )
%                 x = (rand(1,1)-0.5)*D_d(d);
%                 y = (rand(1,1)-0.5)*D_d(d);
%             end
%             source_pt  = [x y H(h)-.001];
            
            % generate point on retroreflector
            x = (rand(1,1)-0.5)*R_d;
            y = (rand(1,1)-0.5)*R_d;
            while (x^2+y^2 > (R_d/2)^2)
                x = (rand(1,1)-0.5)*R_d;
                y = (rand(1,1)-0.5)*R_d;
            end
            
            Xr(i)=x;
            Yr(i)=y;

end
%%
% figure(2)
t=tiledlayout(1,2);


nexttile
plot(Xl,Yl,'b.','MarkerSize',.0001)
title([num2str(N, '%.1e') ' samples of The light source'])
ax=gca;
axis square
grid on

nexttile
plot(Xr,Yr,'r.','MarkerSize',.0001)
title([num2str(N, '%.1e') ' samples of The retroreflector foace'])
xlim(ax.XLim);
ylim(ax.YLim)
axis square
grid on
