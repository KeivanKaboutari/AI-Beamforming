function optimizationFunction(stochasticPhase, evaluationNumber, Band)
    global NFE PhaseValuesPerIteration optFlag Key Criterion M N;
    solution = [];
    objectiveValue = 0;
%     ccons = [];
%     ceqcons = [];
    
    % Pass fixed parameters to objfun
    ObjectiveFunction = @(StochasticPhase) objectiveFunction(StochasticPhase);
%     ConstraintFunction = @(x, y) constraintFunction(x);
    
    % Set nondefault solver options
%     Options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'MaxFunctionEvaluations', Iteration, 'StepTolerance', 1e-50, 'ConstraintTolerance', 1e-50, 'OptimalityTolerance', 1e-50, 'PlotFcn', {'optimplotfvalconstr', 'optimplotstepsize', 'optimplotfirstorderopt'});
%     Options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', evaluationNumber, 'StepTolerance', 1e-6, 'ConstraintTolerance', 1e-50, 'OptimalityTolerance', 1e-50, 'FunctionTolerance', 1e-12, 'PlotFcn', {'optimplotfvalconstr', 'optimplotstepsize', 'optimplotfirstorderopt'});
    Options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', evaluationNumber, 'StepTolerance', 1e-3, 'ConstraintTolerance', 1e-3, 'OptimalityTolerance', 1e-3, 'FunctionTolerance', 1e-3, 'PlotFcn', {'optimplotfvalconstr', 'optimplotstepsize', 'optimplotfirstorderopt'});
    % Solve
    %[solution, objectiveValue] = fmincon(ObjectiveFunction, stochasticPhase, [], [], [], [], LowerBand, UpperBand, @constraintFcn, Options);
    [solution, objectiveValue] = fmincon(ObjectiveFunction, stochasticPhase, [], [], [], [], Band.Lower, Band.Upper, [], Options);
    
    % Update phase
    stochasticPhase = PhaseValuesPerIteration(:, :, NFE);
    
    % Finds phases which are close to the boundaries
    if ~isempty(find(stochasticPhase > Criterion, M * N)) || ~isempty(find(stochasticPhase < -Criterion, M * N))
        stochasticPhase(stochasticPhase > Criterion) = 0;
        stochasticPhase(stochasticPhase < -Criterion) = 0;
        optFlag = 1;
    end
    
    if (optFlag == 1 && Key ~= 2)
        Key = Key + 1;
        NFE = 0;
        %% Call Optimization Function
        OptimizationFunction = @(optimizedPhase, evaluationNumber, Band) optimizationFunction(optimizedPhase, evaluationNumber, Band);
        OptimizationFunction(stochasticPhase, evaluationNumber, Band);
        optFlag = 0;
    end

    % Clear variables
    clearvars options
    
%     [ccons,ceqcons] = ConstraintFunction(solution);
end