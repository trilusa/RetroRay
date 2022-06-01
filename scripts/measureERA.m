function [ERA, xhull, yhull]  = measureERA(x,y,D_d)
    DT=delaunayTriangulation(x,y);
    k=convexHull(DT);
    xhull = DT.Points(k,1);
    yhull = DT.Points(k,2);
    area=abs(trapz(x(k),y(k)));
    ERA = area - pi*(D_d/2)^2;
end

