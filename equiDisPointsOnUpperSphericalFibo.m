function [Theta, Phi] = equiDisPointsOnUpperSphericalFibo(N, maxTiltAngle)
    %% Fibo Method
    % This function distributes N points over the upper part of spherical
    % N is the number of samples which is distributed equidistantly over the upper part of spherical
    
    % Golden ratio in radians
    GR = (sqrt(5) + 1) / 2;
    
    % ???
    Fun = @(X, Y) power((2 * Y + 1) / (2 * N), X);
    
    % Max. value of Theta angle in rad.
    Theta = maxTiltAngle * pi / 180;
    determineMaxThetaAngle = cos(Theta);
    % Calculate z axis values
    z = zeros(1, N);
    for Counter = 1 : 1 : N
        z(Counter) = 1 - (1 - determineMaxThetaAngle) * Fun(1, Counter - 1);
    end
    
    % Calculate Phi values in rad.
    Phi = zeros(1, N);
    for Counter = 1 : 1 : N
        Phi(Counter) = mod((2 * pi * N / GR) * Fun(1, Counter - 1) + 0, 2 * pi);
    end
    % Calculate in plane radius
%     Rz = sqrt(1 - power(z, 2));
    
    % Calculate theta in rad.
    Theta = acos(z);
    
    % Required x, y and z values
%     Points = [Rz .* cos(Phi); Rz .* sin(Phi); z];
end
