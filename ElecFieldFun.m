function RetAbsElecField = ElecFieldFun(ComplexPhaseFactor, UVPoints)
    %% Electric Field Function
    
    global Rf d Lambda;
    [M, N] = size(ComplexPhaseFactor);
    k0 = 2 * pi / Lambda;
    u = UVPoints(:, 1);
    v = UVPoints(:, 2);
    SamplingNumber = numel(u(:, 1));
    
    mValues = linspace(0, M - 1, M);
    nValues = linspace(0, N - 1, N);
    
    RetComElecField = zeros(1, SamplingNumber);
    RetAbsElecField = zeros(1, SamplingNumber);
    Rad = @ (m, n) sqrt(Rf^2 + ((m - (M - 1) / 2) ^ 2 + (n - (N - 1) / 2) ^ 2) * d ^ 2);
    for SampleCounter = 1 : 1 : SamplingNumber
        for mCounter = 1 : 1 : M
            for nCounter = 1 : 1 : N
                RetComElecField(SampleCounter) = RetComElecField(SampleCounter) + sqrt(1 - power(u(SampleCounter, 1) / (k0 * d), 2)) * (ComplexPhaseFactor(mCounter, nCounter) * exp(1i * (mValues(mCounter) - (M - 1) / 2) * u(SampleCounter, 1) + 1i * (nValues(nCounter) - (N - 1) / 2) * v(SampleCounter, 1) - 1i * rem(k0 * Rad(mValues(mCounter), nValues(nCounter)), 2 * pi))) / Rad(mValues(mCounter), nValues(nCounter));
            end
        end
        RetAbsElecField(SampleCounter) = abs(RetComElecField(1, SampleCounter));
    end
end