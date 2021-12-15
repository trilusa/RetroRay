 function alpha =  calcAlphaForPDWithArea(PD_d)
    v= 0:0.01:1;
    h=1.5;
    R_d=.05;
%     PD_d = 0.002:0.002:0.020;
    PD_r = PD_d/2;
    
    n=1;
    L = .0357;
    Ls = .0063;
    
    for d=1:length(PD_r)
        alpha=zeros(length(d),length(v));
        r = R_d + PD_r(d);
        for i=1:length(v) %horizontal displacements to test
            theta = atan(v(i)/h);
            phi = asin(sin(theta)/n);
            D = 2*L*tan(phi);
            Ds = 2*Ls*tan(phi);
            era = (2*r^2) * ( acos((D+Ds)/r) - ((D+Ds)/r)*sin(acos((D+Ds)/r)));
            era = era - pi*PD_r(d)^2;
            alpha(d,i) =  era / (pi*(r^2 - PD_r(d)^2));
        end
         plot(v,real(alpha(d,:)));
    end   
end