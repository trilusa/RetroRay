function Apd_area = measureApd(x,y)
    %split the points into two halves
     %take the boundary of each
    %find area of boundaties
    %add areas together, return

    if ~iscolumn(x) 
        error("Bad input, needs x and y need to be column vects")
    end
    pts = [x y];
    upper=pts(pts(:,2)>=0,:); %points with y greather than zero
    lower=pts(pts(:,2)<0,:);

    lowerB = lower(boundary(lower),:);
    upperB = upper(boundary(upper),:);

%     plot(lowerB(:,1),lowerB(:,2))
%     hold on
%     plot(upperB(:,1),upperB(:,2))

    lowerA = abs(trapz(lowerB(:,1),lowerB(:,2)));
    upperA = abs(trapz(upperB(:,1),upperB(:,2)));

%     scatter(lower(:,1),lower(:,2),'.')
%     scatter(upper(:,1),upper(:,2),'.')
    Apd_area = lowerA + upperA;

   
end