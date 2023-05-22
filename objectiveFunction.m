function ReturnCost = objectiveFunction(randomPhase, Lx, Ly) %#codegen
%% Cost function to minimize Sidelobe Levels (SLL)
%     global NFE PhaseValuesPerIteration PatternE PatternENormalized PatternESmooth PatternESmoothedNormalized CostValuePerIter;
    global NFE PhaseValuesPerIteration PatternESmooth PatternESmoothedNormalized CostValuePerIter;
    global M N;
%     global kValues lValues;
    global SmoothingFactor kValuesSmoothed lValuesSmoothed
    global FixedOutsideSamplingPointsValues FixedInsideSamplingPointsValues OutsideSamplingPoints InsideSamplingPoints;
    global L ControlFactor AbsElectricFieldSmoothedMax CorrectionFactorCPR Coeff1 Coeff2;
%     global Criterion;
    
    NFE = NFE + 1;
    
%     % Finds phases which are close to the boundaries
%     if ~isempty(find(randomPhase > Criterion, M * N))
%         randomPhase(1, find(randomPhase > Criterion, M * N)) = 0;
%     end
%     if ~isempty(find(randomPhase < -Criterion, M * N))
%         randomPhase(1, find(randomPhase < -Criterion, M * N)) = 0;
%     end

    stochasticPhase = reshape(randomPhase, N, M);
    
    % Confine Function is confining the phase within the [-pi, pi] interval.
    Confine = @(x) pi - mod(pi - x, 2 * pi);
%     Confine = @(x) pi - rem(pi - x, 2 * pi);
    
    confinedStochasticPhase = Confine(stochasticPhase);
    
    complexPhaseFactor = 1 * exp(1i * confinedStochasticPhase);
    
    % Save new phases
    PhaseValuesPerIteration(:, :, NFE) = confinedStochasticPhase;
    
    % Calculate Electrical Field
%     AbsElectricField = zeros(N, M);
%     for kCounter = 1 : 1 : M
%         u = -pi + kValues(kCounter) * 2 * pi / (M - 1);
%         for lCounter = 1 : 1 : N
%             v = -pi + lValues(lCounter) * 2 * pi / (N - 1);
%             AbsElectricField(lCounter, kCounter) = ElecFieldFun(complexPhaseFactor, [u, v], Lx, Ly);
%         end
%     end
%     % Element-base Electrical Field
%     ElemBaseElectricFieldMax = max(AbsElectricField, [], 'all');
%     % Save new electrical patterns
%     PatternE(:, :, NFE) = AbsElectricField(:, :);
%     PatternENormalized(:, :, NFE) = AbsElectricField(:, :) / ElemBaseElectricFieldMax;
    
    % Calculate Smoothed Electrical Field
    AbsElectricFieldSmoothed = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
    for kCounter = 1 : 1 : round(SmoothingFactor * M)
        u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
        for lCounter = 1 : 1 : round(SmoothingFactor * N)
            v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
            AbsElectricFieldSmoothed(lCounter, kCounter) = ElecFieldFun(complexPhaseFactor, [u, v], Lx, Ly);
        end
    end
    % Smoothed Electrical Field
    PatternESmoothedNormalizedMax = max(AbsElectricFieldSmoothed, [], 'all');
    % Save new smoothed electrical patterns
    PatternESmooth(:, :, NFE) = AbsElectricFieldSmoothed(:, :);
    PatternESmoothedNormalized(:, :, NFE) = AbsElectricFieldSmoothed(:, :) / PatternESmoothedNormalizedMax;
    
    % Determine cost function and save costs per iteration
    % AbsElecFieldOutsideProArea = ElecFieldFun(complexPhaseFactor, [FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, 2), FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, 1)], Lx, Ly);
    AbsElecFieldOutsideProAreaFiltered = ElecFieldFun(complexPhaseFactor, [FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, 1), FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, 2)], Lx, Ly) .* transpose(modifiedGaussianFilter(FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, 1), FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, 2), 0));
    % AbsElecFieldOutsideProAreaFiltered = ElecFieldFun(complexPhaseFactor, [FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, 1), FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, 2)], Lx, Ly);
    EngMaxEF = ceil(AbsElectricFieldSmoothedMax + 5);
%     AbsElecFieldInsideProArea = zeros(1, InsideSamplingPoints, L);
    AbsElecFieldInsideProAreaFiltered = zeros(1, InsideSamplingPoints, L);
    CRPInsideProArea = zeros(1, InsideSamplingPoints, L);
%     InsideProAreaCostFun = zeros(1, L);
    InsideProAreaCostFunFiltered = zeros(1, L);
    for BeamCounter = 1 : 1 : L
        % AbsElecFieldInsideProArea(:, :, BeamCounter) = ElecFieldFun(complexPhaseFactor, [FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, 1, BeamCounter), FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, 2, BeamCounter)], Lx, Ly);
        AbsElecFieldInsideProAreaFiltered(:, :, BeamCounter) = ElecFieldFun(complexPhaseFactor, [FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, 1, BeamCounter), FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, 2, BeamCounter)], Lx, Ly) .* transpose(modifiedGaussianFilter(FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, 1, BeamCounter), FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, 2, BeamCounter), BeamCounter));
        % AbsElecFieldInsideProAreaFiltered(:, :, BeamCounter) = ElecFieldFun(complexPhaseFactor, [FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, 1, BeamCounter), FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, 2, BeamCounter)], Lx, Ly);
        for PointCounter = 1 : 1 : InsideSamplingPoints
            CRPInsideProArea(1, PointCounter, BeamCounter) = ComplexRadiationPatternUV(FixedInsideSamplingPointsValues(PointCounter, 1, BeamCounter), FixedInsideSamplingPointsValues(PointCounter, 2, BeamCounter));
        end
%         InsideProAreaCostFun(1, BeamCounter) = sum(power(abs(AbsElecFieldInsideProArea(:, :, BeamCounter) - ControlFactor(BeamCounter) .* CorrectionFactorCPR(BeamCounter) .* BeamMaxEF(BeamCounter) .* CRPInsideProArea(1, :, BeamCounter)), 1));
%         InsideProAreaCostFun(1, BeamCounter) = sum(power(abs(AbsElecFieldInsideProArea(:, :, BeamCounter) - ControlFactor(BeamCounter) .* CorrectionFactorCPR(BeamCounter) .* EngMaxEF .* CRPInsideProArea(1, :, BeamCounter)), 1.3));
        InsideProAreaCostFunFiltered(1, BeamCounter) = sum(power(abs(AbsElecFieldInsideProAreaFiltered(:, :, BeamCounter) - ControlFactor(BeamCounter) .* CorrectionFactorCPR(BeamCounter) .* EngMaxEF .* CRPInsideProArea(1, :, BeamCounter)), 1.3));
    end
%     OutsideProAreaCostFun = @() sum(power(AbsElecFieldOutsideProArea, 2));
    OutsideProAreaCostFun = @() sum(power(AbsElecFieldOutsideProAreaFiltered, 2));
    
    % Save new cost function evaluation value
    CostValuePerIter(:, NFE) = Coeff1 * OutsideProAreaCostFun() + Coeff2 * sum(InsideProAreaCostFunFiltered);
    
    % Return cost function evaluation value
    ReturnCost = CostValuePerIter(:, NFE);
end