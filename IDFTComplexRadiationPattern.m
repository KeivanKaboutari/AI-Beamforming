function RetSincPatternPhase = IDFTComplexRadiationPattern(ComplexRadiationPattern, Lx, Ly)
    %% Inverse Discrete Fourier Transform from complex values and pattern of the antenna
    
    global Rf d k0;
    [N, M] = size(ComplexRadiationPattern);
    
    Rad = @ (n, m) sqrt(Rf^2 + ((n - (N - 1) / 2 + Lx) ^ 2 + (m - (M - 1) / 2) ^ 2 + Ly) * d ^ 2);
    
    mValues = linspace(0, M - 1, M);
    nValues = linspace(0, N - 1, N);
    kValues = linspace(0, M - 1, M);
    lValues = linspace(0, N - 1, N);
    
    RetSincPatternPhase = zeros(N, M);
    for mCounter = 1 : 1 : M
        for nCounter = 1 : 1 : N
            for kCounter = 1 : 1: M
                for lCounter = 1 : 1 : N
%                     RetSincPatternPhase(mCounter, nCounter) = RetSincPatternPhase(mCounter, nCounter) + ComplexRadiationPattern(kCounter, lCounter) * exp(-2 * pi * 1i * ((mValues(mCounter) * (kValues(kCounter) - (M - 1) / 2) / M) + (nValues(nCounter) * (lValues(lCounter) - (N - 1) / 2) / N)));
                    RetSincPatternPhase(nCounter, mCounter) = RetSincPatternPhase(nCounter, mCounter) + ComplexRadiationPattern(lCounter, kCounter) * exp(-2 * pi * 1i * ((nValues(nCounter) - (N - 1) / 2) * (lValues(lCounter) - (N - 1) / 2) / N + (mValues(mCounter) - (M - 1) / 2) * (kValues(kCounter) - (M - 1) / 2) / M));
                end
            end
            % Just pick up the phase information (k0 * Rad(m, n) is coming due to the additional phase
            % delay that happens during propagation (all plane waves are not normal to the MS plane))
            RetSincPatternPhase(nCounter, mCounter) = angle(RetSincPatternPhase(nCounter, mCounter)) + rem(k0 * Rad(nValues(nCounter), mValues(mCounter)), 2 * pi);
        end
    end
end