v = 0.0:.02:1;   %horizonal positions
D_d = .004:.004:.016;  %diode diameters to test
Pt=1;  % Tx power, 1 watt
h=1.5; %m
m=1; %lambertian source with half angle = 60 degrees
H = 1.5;
Rd=50e-3;
L = .0357 * (Rd/50e-3);
Ls = .0063 * (Rd/50e-3);


%initialize variables to hold results
Pr = zeros(length(D_d), length(v)); % holds absolute power values
Pr_db = zeros(length(D_d), length(v)); % holds power in dB wrt power recived at v=0

%for every diode and position, calculate the receiever power
for d=1:length(D_d)
As = pi*(D_d(d)/2)^2; %calc sensing area
for i = 1:length(v)
    alpha = calcAlphaForPDWithArea(D_d(d),v(i),H,Rd,L,Ls); %calc alpha
    dist = sqrt(v(i)^2+h^2); %distance from light to corcnercube
    theta = atan(v(i)/h); %incidence angle wrt CC normal

    %This is eqn 7 in info com paper
    Pr(d,i)= Pt * (((m+1)*As) / (8*pi*dist^2)) * (cos(theta)^(m+1)) * alpha;
    Pr_db(d,i) = 20*log10(Pr(d,i)/Pr(d,1));

   %power calc goes negative when ERA = 0, so log10() become imaginary
   if imag(Pr_db(d,i)) ~= 0
       Pr_db(d,i) = -inf;
   end

end
end

%%plotting
clf
t=tiledlayout('flow');
nexttile
plot(v,Pr(1:4:end,:),'-')
lgn_txt = repmat("PD diam = ", length(D_d(1:4:end)), 1) + num2str(D_d(1:4:end)'*1000) + repmat(" mm",length(D_d(1:4:end)),1);
lgd=legend(lgn_txt);
title("Horizontal Dist vs Rx Power")
% ylabel("P (dB_{FS})")
ylabel("P (W)")

xlabel("Horizontal Displacement (m)")
set(lgd,'location','northeast')
grid on


nexttile
plot(D_d*1000,Pr(:,1:10:end))
lgn_txt = repmat("v = ", length(v(1:10:end)), 1) + num2str(v(1:10:end)'*1000) + repmat(" mm",length(v(1:10:end)),1);
lgd=legend(lgn_txt);
title("PD diam vs Rx Power")
% ylabel("P (dB_{FS})")
ylabel("P (W)")

set(lgd,'location','northwest')
xlabel("PD Diameter (mm)")
grid on