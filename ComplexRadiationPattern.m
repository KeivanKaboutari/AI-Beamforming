function RetComplexRadiationPattern = ComplexRadiationPattern(UValues, VValues, amp, phase)
    %% Complex Radiation Pattern, Complex value of radiation pattern's samples at different values of u and v (u and v are in Spherical Coordinate)
    
    global sx sy u0 v0 k0 d;
    
    % Pattern of the antenna (EXAMPLE). u and v are scalar values.
    Pattern = @(u, v) sum(amp .* exp(1j .* phase) .* sinc(sx * (u - u0) / (2 * pi)) .* sinc(sy * (v - v0) / (2 * pi)));
    
    % Complex value of radiation pattern's samples at different values of u and v
    % k and l are number of elements in each Column and Row, respectively.
    RetComplexRadiationPattern = zeros(numel(VValues), numel(UValues));
    for lCounter = 1 : 1 : numel(VValues)
        for kCounter = 1 : 1 : numel(UValues)
            RetComplexRadiationPattern(lCounter, kCounter) = sqrt(1 - power(UValues(kCounter) / (k0 * d), 2)) * Pattern(UValues(kCounter), VValues(lCounter));
        end
    end
end