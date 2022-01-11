%  function alpha =  calcEffectiveReflectingArea(v,h,r,n,L,Ls)
    v= 0:0.01:1
    h=1.5
    r=.05
    n=1
    L = .0357
    Ls = .0063

    alpha=zeros(length(v),1);
    for i=1:length(v) %horizontal displacements to test
        theta = atan(v(i)/h)
        phi = asin(sin(theta)/n)
        D = 2*L*tan(phi)
        Ds = 2*Ls*tan(phi)
        alpha(i) = (2*r^2) * ( acos((D+Ds)/r) - ((D+Ds)/r)*sin(acos((D+Ds)/r))) / (pi*r^2)
    end
%  end