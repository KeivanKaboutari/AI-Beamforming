function RetComplexRadiationPattern = ComplexRadiationPatternUV(UValue, VValue)
    %% Complex Radiation Pattern, Complex value of radiation pattern's samples at different values of u and v (u and v are in Spherical Coordinate)
    global sx sy u0 v0;
    
    % Pattern of the antenna (EXAMPLE). u and v are scalar values.
    Pattern = @(u, v) sum(sinc(sx * (u - u0) / (2 * pi)) .* sinc(sy * (v - v0) / (2 * pi)));
    
    % Complex value of radiation pattern's samples at different values of u and v
    RetComplexRadiationPattern = Pattern(UValue, VValue);
end