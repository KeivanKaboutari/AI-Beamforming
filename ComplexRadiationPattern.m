function RetComplexRadiationPattern = ComplexRadiationPattern(UValues, VValues)
    %% Complex Radiation Pattern, Complex value of radiation pattern's samples at different values of u and v (u and v are in Spherical Coordinate)
    
    global d Lambda sx sy u0 v0;
    
    k0 = 2 * pi / Lambda;

    % Pattern of the antenna (EXAMPLE). u and v are scalar values.
    Pattern = @(u, v) sum(sinc(sx * (u - u0) / (2 * pi)) .* sinc(sy * (v - v0) / (2 * pi)));
    
    % Complex value of radiation pattern's samples at different values of u and v
    % k and l are number of elements in each Row and Column, respectively.
    RetComplexRadiationPattern = zeros(numel(UValues), numel(VValues));
    for kCounter = 1 : 1 : numel(UValues)
        for lCounter = 1 : 1 : numel(VValues)
%             RetComplexRadiationPattern(kCounter, lCounter) = Pattern(UValues(kCounter), VValues(lCounter));
            RetComplexRadiationPattern(kCounter, lCounter) = sqrt(1 - power(UValues(kCounter) / (k0 * d), 2)) * Pattern(UValues(kCounter), VValues(lCounter));
        end
    end
end