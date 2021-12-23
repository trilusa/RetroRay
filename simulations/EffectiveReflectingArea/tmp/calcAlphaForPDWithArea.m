 function alpha =  calcAlphaForPDWithArea(PD_d,v)
%     v= 0:0.01:5;
    h=1.5;
    R_d=.05;
%     PD_d = .01;%0.002:0.002:0.020;
    PD_r = PD_d/2;
    
    n=1;
    L = .0357;
    Ls = .0063;
    
    alpha=zeros(length(PD_r),length(v));
    for d=1:length(PD_r)
        r = R_d + PD_r(d);
        for i=1:length(v) %horizontal displacements to test
            theta = atan(v(i)/h);
            phi = asin(sin(theta)/n);
            D = 2*L*tan(phi);
            Ds = 2*Ls*tan(phi);
            era = (2*r^2) * ( acos((D+Ds)/r) - ((D+Ds)/r)*sin(acos((D+Ds)/r)));
            era = era - pi*PD_r(d)^2;
            alpha(d,i) = real(  era / (pi*(r^2 - PD_r(d)^2)));
        end
%           plot(v,real(alpha(d,:)));
    end   
end