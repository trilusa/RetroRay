function [alpha, era, Apd] =  calcAlphaForPDWithArea(Dia_PD,V,H,Dia_CC,L,Ls)
% make sure only one paramter is changing, for simplicity
if(length(Dia_PD) > 1)
    if(length(H)==1 && length(Dia_CC)==1)
        era=ones(length(Dia_PD),length(V));
        alpha=ones(length(Dia_PD),length(V));
        Apd = ones(length(Dia_PD),length(V));
    else
        error('one independent var please')
    end
end
if(length(H) > 1)
    if(length(Dia_PD)==1 && length(Dia_CC)==1)
        era=ones(length(H),length(V));
        alpha=ones(length(H),length(V));
        Apd=ones(length(H),length(V));

    else
        error('one independent var please')
    end
end
if(length(Dia_CC) > 1)
    if(length(H)==1 && length(Dia_PD)==1)
        era=ones(length(Dia_CC),length(V));
        alpha=ones(length(Dia_CC),length(V));
        Apd=ones(length(Dia_CC),length(V));

    else
        error('one independent var please')
    end
else
     era=ones(1,length(V));
     alpha=ones(1,length(V));
      Apd=ones(1,length(V));
end


%check for consistency of retroreflector dimensions
if((length(Dia_CC) ~= length(L)) || (length(Dia_CC) ~= length(Ls)))
    error('inconsistent retroreflector dimensions');
end

R_PD = Dia_PD./2;
n=1; %refractive index
%     L = .0357;
%     Ls = .0063;
i=1;
for r = 1:length(Dia_CC)
    for h=1:length(H)
        for d=1:length(Dia_PD)
            RdDr = Dia_CC(r) + R_PD(d);
            for v=1:length(V) %horizontal displacements to test
                V(v) = max( .5*(Dia_CC(r)+Dia_PD(d)), V(v)); %adjust for points directly under light

                %calculate ERA
                theta = atan(V(v)/H(h));
                phi = asin(sin(theta)/n);
                D = 2*L(r)*tan(phi);
                Ds = 2*Ls(r)*tan(phi);
                era(i,v) = real(((2*RdDr^2) * ( acos((D+Ds)/RdDr) - ((D+Ds)/RdDr)*sin(acos((D+Ds)/RdDr)))) - pi*R_PD(d)^2);
                
                %Calculate alpha
                alpha(i,v) =  era(i,v) / (pi*(RdDr^2 - R_PD(d)^2));

                %Calculate A_pd according to derivation in paper
                 %Case 1
                if(Dia_CC(r)<.5*Dia_PD(d))
                    phi_a = acos((D+Ds)/Dia_CC(r));
                    r_sec = .5*Dia_PD(d)-Dia_CC(r);
                    L_seg = 2*(D+Ds)+2*r_sec*cos(phi_a);
                    phi_b = acos(1-(2*L_seg^2)/(Dia_PD(d)^2));

                    A_sector = phi_a*r_sec^2;
                    A_trapezoid = .5*r_sec*sin(phi_a)*(2*(D+Ds)+L_seg);
                    A_segment = (1/8)*Dia_PD(d)^2*(phi_b-sin(phi_b));

                    Apd(i,v) = .25*pi*Dia_PD(d)^2 - 2*(A_sector+A_trapezoid+A_segment);
                %Case 2    
                elseif( (Dia_CC(r)>=.5*Dia_PD(d)) && sqrt(Dia_CC(r)^2 - ((D+Ds)^2)) < .5*Dia_PD(d))
                    phi_b=2*acos((2*sqrt(Dia_CC(r)^2-(D+Ds)^2)) / Dia_PD(d));
                    A_segment= (1/8)*Dia_PD(d)^2*(phi_b-sin(phi_b));

                    Apd(i,v)=.25*pi*Dia_PD(d)^2 - 2*A_segment;
                %Case 3
                elseif(sqrt(Dia_CC(r)^2 - ((D+Ds)^2)) >= .5*Dia_PD(d))
                    Apd(i,v)=.25*pi*Dia_PD(d)^2;
                end
            end
            i=i+1;
        end
    end
end
end

