function [thetaInRad, phiInRad] = equiDisPointsOnUpperSphericalCube(maxTiltAngle)
    %% Cube Projection
    
    % Tangent wrap function with parameter Theta.
    % This function distorts the rectangular grid on the cube face.
    TRF = @(Point, Theta) tan(Point * Theta) / tan(Theta);
    
    % Optimal value as given in literature
    Theta = 0.8687;
    % Max. value of Theta angle in rad.
    ThetaMax = maxTiltAngle * pi / 180;
    
    % Grid no. over a and b axes
%     % For 65 Points
%     Na = 5.5;
%     Nb = 5.5;
    % For 32 Points
    Na = 4.0;
    Nb = 4.0;
    
    % Griding a and b axes
    n = 0 : 1 : Na - 1;
    m = 0 : 1 : Nb - 1;
    
    % Grids over a and b axes
    Sa = 1 - (2 * n + 1) / Na;
    Sb = 1 - (2 * m + 1) / Nb;
    
    % find griding intersection points
    Mesh = [];
    for OutCounter = 1 : 1 : length(Sa)
        for InCounter = 1 : 1 : length(Sb)
            Mesh = [Mesh; Sa(OutCounter) Sb(InCounter)];
        end
    end
    
    % Find points over the cube
    x = 1;
    y = TRF(Mesh(:, 1), Theta);
    z = TRF(Mesh(:, 2), Theta);
    r = sqrt(1 + power(y, 2) + power(z, 2));
    
    % Normalize values regarding assumed circumferenced spherical by cube
    X = 1 ./ r;
    Y = y ./ r;
    Z = z ./ r;
    
    % Sampled points over one of the cubes faces
    Face = [X Y Z];
    
    % This section performs rotations corresponding to the 6 faces of a cube.
    % Input: a list of points on the cube's face centered at x = 1, y = z = 0.
    % Output: a list of points on all faces.
    spherePoints = [ Face; ...
                   [ Face(:, 2) -Face(:, 1)  Face(:, 3)]; ...
                   [-Face(:, 1) -Face(:, 2)  Face(:, 3)]; ...
                   [ Face(:, 2)  Face(:, 1)  Face(:, 3)]; ...
                   [-Face(:, 3)  Face(:, 2)  Face(:, 1)]; ...
                   [ Face(:, 3)  Face(:, 2) -Face(:, 1)]];
    
    % Filter points based on a given condition.
    % Finds points over the upper part of spherical coordinates which their Theta angle are equal or less than 75 degree
    [pointsRow, pointsColumn] = size(spherePoints);
    FilteredPoints = [];
    for rowCounter = 1 : 1 : pointsRow
        if spherePoints(rowCounter, 3) > cos(ThetaMax)
            FilteredPoints = [FilteredPoints; spherePoints(rowCounter, :)];
        end
    end

    % Convert Cartesian to Spherical coordinate
    [phiInRad, thetaInRad, r] = cart2sph(FilteredPoints(:, 1), FilteredPoints(:, 2), FilteredPoints(:, 3));
    % Correct elevation angle (change its range from [(-pi / 2) (pi / 2)] to [0 pi])
    thetaInRad = transpose(pi / 2 - thetaInRad);
    phiInRad = transpose(phiInRad);
    
end
