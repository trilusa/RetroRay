function [refl1, refl2, refl3, cylender, circle] = buildCornerCube(R_d, R_L, R_Ls, R_pos,  R_az, R_el )
    Rx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
    Ry = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
    Rz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
    
    %% calculate side length from params
     s = .068; %This is for a 50mm cube
     %so I will scale it based on what R_d is
     s = .068 * (R_d / 50e-3);

%   initial  veritices along bases
    vertex=[0; 0; 0];
    face=s*eye(3);

    face = Rx(pi/5)*Ry(-pi/4)*face;     %orient facing up
    face = face - [0; 0; R_L + R_Ls];          %shift down
    vertex = vertex - [0; 0; R_L + R_Ls];
    face = face + R_pos;    %shift to test position
    vertex = vertex + R_pos;
     
    face = Rz(R_az)*Ry(R_el)*face;      %rotate to specified elev and az

%     centroid = (face(:,1)+face(:,2)+face(:,3))/3;

    [nx, ny, nz] = sph2cart(R_az, R_el+(pi/2), 1);
    n = [nx; ny; nz];

    %  return reflectors/absorbers
    cylender = CylinderAperature(R_d/2, R_Ls, R_pos-[0;0;R_Ls], n);
    circle=CircleAperature(R_d/2, R_pos, n);
    refl1 = TriangleReflector(vertex, face(:,1), face(:,2));
    refl2 = TriangleReflector(vertex, face(:,3), face(:,1));
    refl3 = TriangleReflector(vertex, face(:,2), face(:,3));
end