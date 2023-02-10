clear all;
close all;
clc;

%% Define global varibales

% Number of cost function evaluation (NFE)
global NFE;

global u0 v0 ControlFactor;
global CorrectionFactorCPR;
global ComplexPhaseFactor;
global AbsElectricFieldSmoothedMax;
global OutsideSamplingPoints InsideSamplingPoints;
global FixedOutsideSamplingPointsValues FixedInsideSamplingPointsValues;
global InsideProtectedAreaScatteredPointsProperties;
global PhaseValuesPerIteration;
global optFlag Key;
global Coeff1 Coeff2;
global CostValuePerIter;
global PatternESmooth PatternESmoothedNormalized;

%% Control functionality and output of the program
% When the key is equal to 1, it will represent the results
% When the key is equal to 0, it will not represent the results

% Plot not smoothed Complex Radiation Pattern (not smoothed)
NCRPKey = 0;

% Plot smoothed Complex Radiation Pattern
SCRPKey = 0;

% Plot phase distribution
PhaseKey = 0;

% Plot not smoothed electric field
NEFieldKey = 0;

% Plot smoothed electric field
SEFieldKey = 0;

% Plot sampling area
SAreaKey = 0;

% Run optmization
OptRunKey = 1;

% Plot optimized results
OptReskey = 0;

%% In this example calculation, the time dependence is exp(+i * ω * t), as is generally assumed in electrical engineering.
% Therefore, when an EM wave propagates along a distance Rad in free space it acquires a phase delay given by the phase factor exp(-i * k0 * Rad),
% where k0 = ω / c = 2 * pi / Lambda (wave number or angular wavenumber).

% Below, Rf is the focal distance, which is the distance from the source antenna to the MS.
% The MS is formed by N * M elements (reflecting or transmitting, which does not matter for this calculation).
% Let us define the MS parameters: Rf, M, N, d
global Rf;
Rf = 50.000000;
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Rf.mat', 'Rf');
% Number of elements in each Columns
global M;
M = 10;
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\M.mat', 'M');
% Number of elements in each Rows
global N;
N = 3;
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\N.mat', 'N');

% Operational frequency
Frequency = 5e9;

% Determine radius of Protected Areas (PAs)
% Increasing or decreasing Q's value increases or decreases the radiuses of PAs.
Q = 1 / 4;
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Q.mat', 'Q');

% Let us define the sinc(x) function as a pattern
% Let us now define some radiation pattern we want to realize. In particular, it makes sense to consider patterns
% given by a superposition of sinc(sx * (u - u0) / 2) * sinc(sy * (v - v0) / 2) terms, as these are the basic patterns
% of a rectangular aperture with the size sx * d by sy * d. Thus for the radiation pattern with L beams we have
% Let us define the radiation pattern parameters with L beams:
global L;
L = 1;
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\L.mat', 'L');

% Rad(n, m) is the radial distance from the source antenna to the (m, n)-th element of the MS.
% d [cm] is the period of the MS (the distance between the centers of its elements), which is the same along x and y axes.
global d;
d = 0.68;
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\d.mat', 'd');

% m and n are scalar values which presents mth and nth element of the metasurface array
% Rad = @ (n, m) sqrt(Rf^2 + ((n - (N - 1) / 2) ^ 2 + (m - (M - 1) / 2) ^ 2) * d ^ 2);

% Here, F(u, v) is the (arbitrary) radiation pattern function expressed in the reciprocal space with the coordinates
% u = k0 * d * sin(θ) * cos(φ), v = k0 * d * sin(θ) * sin(φ). Because the far-field radiation pattern
% of a rectangular aperture in such coordinates is equivalent to the Fourier transform of the aperture wave field,
% the necessary amplitude and phase distribution on the MS can be found with the inverse Fourier transform
% (we use the discrete version of it). In the following, we neglect the magnitude and just consider the phase.
% Besides the phase distribution dictated by the radiation pattern, we must compensate for the phase delay
% due to the propagation from the source antenna to a given point on the MS (the last term).

% Aperture size along X and Y axes (Sampling number per axes)
% Make beam narrower ==> Increase sx and sy
% Make beam wider    ==> decrease sx and sy
global sx sy;
sx = M / 1;
sy = N / 1;

% Let us define the free space wavelength (Lambda [cm]) and wavenumber:
global Lambda;
Lambda = physconst('LightSpeed') / Frequency * 100;
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Lambda.mat', 'Lambda');
global k0;
k0 = 2 * pi / Lambda;

% MS number of elements
elementNo = N * M;

% Generate Theta and Phi samples which is distributed equidistantly over the upper part of spherical
% Select point distribution over upper half hemisphere equidistantly
% 1 for selecting Fibo method and 2 for selecting Cube method
PointDistributionMethod = 1;
% Determine max tilt angle for beams [Degree]
maxTiltAngle = 75;
if PointDistributionMethod == 1
    angleSamplingNo = 32;
    [Theta, Phi] = equiDisPointsOnUpperSphericalFibo(angleSamplingNo, maxTiltAngle);
elseif PointDistributionMethod == 2
    [Theta, Phi] = equiDisPointsOnUpperSphericalCube(maxTiltAngle);
    angleSamplingNo = length(Theta);
end
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Theta.mat', 'Theta');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Phi.mat', 'Phi');

% Number of desired data sets
noDesData = 100;

% Smoothing Factor to make patterns more smoother
global SmoothingFactor;
SmoothingFactor = 5;

% Get values around the u0 and v0 to find main beams amplitude
DistanceFromCenter = 4;

% Criterian of the optimization algorithm to change some of phase values which avoids criteriaon
global Criterion;
Criterion = 3;

%% Define vectors to save data
% Theta, Phi, N * M Phases (without noise)
thetaPhiPhaseWoN = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% Theta, Phi, N * M Phases (with noise)
thetaPhiPhaseWN = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% Theta, Phi, N * M Phases (optimized)
thetaPhiPhaseOpt = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% N * M Phases (optimized and modified)
modifiedPhaseOpt = zeros(angleSamplingNo, elementNo, noDesData);
% Theta, Phi, N * M CRP
% thetaPhiPatCR = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% Theta, Phi, N * M CRP (normalized)
% thetaPhiPatCRNormalized = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% (M * SmoothingFactor) * (N * SmoothingFactor) CRP
thetaPhiPatCRSmoothed = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% (M * SmoothingFactor) * (N * SmoothingFactor) CRP (normalized)
thetaPhiPatCRSmoothedNormalized = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% Theta, Phi, N * M E-pattern
% thetaPhiPatE = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% Theta, Phi, N * M E-pattern (normalized)
% thetaPhiPatENormalized = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% Theta, Phi, N * M E-patter (optimized)
% thetaPhiPatEOpt = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% Theta, Phi, N * M E-patter (normalized and optimized)
% thetaPhiPatEOptNormalized = zeros(angleSamplingNo, 2 * L + elementNo, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) E-pattern
thetaPhiPatESmoothed = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) E-pattern (normalized)
thetaPhiPatESmoothedNormalized = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) E-pattern (optimized)
thetaPhiPatESmoothedOpt = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) E-pattern (normalized and optimized)
thetaPhiPatESmoothedOptNormalized = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);

% Initialize Theta and Phi values into the related arrays
thetaPhiPhaseWoN(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);
thetaPhiPhaseWN(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);
thetaPhiPhaseOpt(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);

% thetaPhiPatCR(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);
% thetaPhiPatCRNormalized(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);

% thetaPhiPatE(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);
% thetaPhiPatENormalized(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);
% thetaPhiPatEOpt(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);
% thetaPhiPatEOptNormalized(:, 1 : 2, :) = repmat([transpose(Theta), transpose(Phi)], [1 1 noDesData]);

%% Generate meshes and grids for 2D and 3D plots

% Mesh regarding Massive MIMO Antenna rows and columns indices (Element Base)
% Column
ElemBaseMeshX = linspace(-floor(M / 2), floor(M / 2), M);
ElementsX = linspace(1, M, M);
% Row
ElemBaseMeshY = linspace(-floor(N / 2), floor(N / 2), N);
ElementsY = linspace(1, N, N);
% Generate the mesh grid
[ElemBaseX, ElemBaseY] = meshgrid(ElemBaseMeshX, ElemBaseMeshY);
[gridElemBaseX, gridElemBaseY] = meshgrid(ElementsX, ElementsY);

% k and l are number of elements in each Column and Row, respectively.
global kValues lValues;
kValues = linspace(0, M - 1, M);
lValues = linspace(0, N - 1, N);
% Convert coordinate to UV space
% Column
UVBaseMeshX = 2 * k0 * d * (kValues - (M - 1) / 2) / (M - 1);
% Row
UVBaseMeshY = 2 * k0 * d * (lValues - (N - 1) / 2) / (N - 1);
% Generate mesh grid
[gridUVBaseX, gridUVBaseY] = meshgrid(UVBaseMeshX, UVBaseMeshY);

% Smoothed mesh regarding Massive MIMO Antenna
global kValuesSmoothed lValuesSmoothed;
kValuesSmoothed = linspace(0, M - 1, round(SmoothingFactor * M));
lValuesSmoothed = linspace(0, N - 1, round(SmoothingFactor * N));
% Convert coordinate to UV space
global UVBaseMeshXSmoothed UVBaseMeshYSmoothed gridUVBaseXSmoothed gridUVBaseYSmoothed;
% Column
UVBaseMeshXSmoothed = 2 * k0 * d * (kValuesSmoothed - (M - 1) / 2) / (M - 1);
% Row
UVBaseMeshYSmoothed = 2 * k0 * d * (lValuesSmoothed - (N - 1) / 2) / (N - 1);
% Generate the mesh grid
[gridUVBaseXSmoothed, gridUVBaseYSmoothed] = meshgrid(UVBaseMeshXSmoothed, UVBaseMeshYSmoothed);

% For finding exact loacation of beams and plotting E-field in uv-space
MaxLocUVBaseMeshX = 2 * pi * (kValues - (M - 1) / 2) / (M - 1);
MaxLocUVBaseMeshY = 2 * pi * (lValues - (N - 1) / 2) / (N - 1);
% Generate the mesh grid for finding exact loacation of beams and plotting E-field in uv-space
% [gridMaxLocUVBaseX, gridMaxLocUVBaseY] = meshgrid(MaxLocUVBaseMeshX, MaxLocUVBaseMeshY);

% For plotting E-field in uv-space and finding exact loacation of beams
MaxLocUVBaseMeshXSmoothed = 2 * pi * (kValuesSmoothed - (M - 1) / 2) / (M - 1);
MaxLocUVBaseMeshYSmoothed = 2 * pi * (lValuesSmoothed - (N - 1) / 2) / (N - 1);
% Generate mesh grid for plotting E-field in uv-space and finding exact loacation of beams
global gridMaxLocUVBaseXSmoothed gridMaxLocUVBaseYSmoothed;
[gridMaxLocUVBaseXSmoothed, gridMaxLocUVBaseYSmoothed] = meshgrid(MaxLocUVBaseMeshXSmoothed, MaxLocUVBaseMeshYSmoothed);

% Calculate distance from center of the visible area (indeed, MS)
centerDistance = sqrt(power(gridMaxLocUVBaseXSmoothed, 2) + power(gridMaxLocUVBaseYSmoothed, 2));

% Generate grid to plot in real dimension on the MS
SizeEleX = linspace(-(M - 1) * d / 2, (M - 1) * d / 2, M);
SizeEleY = linspace(-(N - 1) * d / 2, (N - 1) * d / 2, N);

%% Define visible area
% Center of the visible area
U0 = 0;
V0 = 0;
% Radius of visible area in UV-space
visibleAreaRadius = k0 * d;
% Define Theta to plot circle for visible area in UV-space
global visibleAreaTheta;
visibleAreaTheta = 0 : pi / 100 : 2 * pi;
% Generate grid for visible area in UV-space
global visibleAreaX visibleAreaY;
visibleAreaX = visibleAreaRadius .* cos(visibleAreaTheta) + V0;
visibleAreaY = visibleAreaRadius .* sin(visibleAreaTheta) + U0;

%% Define elliptic protected area
% Protected area (ovals)
global ProtectedAreaPhi;
ProtectedAreaPhi = 0 : pi / 100 : 2 * pi;

 % Calculate the HPBW (Q = 1 / 10 is equal to -20 dB for (u, v))
syms uvSymbole;
% uEqu = sinc(sx * (uvSymbole - u0) / 2) - Q == 0;
uEqu = sinc(uvSymbole) - Q == 0;
uHPBWatOrigin = double(abs(vpasolve(uEqu, uvSymbole)));

% vEqu = sinc(sy * (uvSymbole - v0) / 2) - Q == 0;
vEqu = sinc(uvSymbole) - Q == 0;
vHPBWatOrigin = double(abs(vpasolve(vEqu, uvSymbole)));

global EFXunit EFYunit;
EFXunit = zeros(1, length(ProtectedAreaPhi), L);
EFYunit = zeros(1, length(ProtectedAreaPhi), L);

%% Set axes color for figures
global left_color right_color;
left_color = [0 0 0];
right_color = [0 0 0];

%% Beamforming and optimization
for dataSetCounter = 1 : 1 : noDesData
    
    % Add Zero Mean Unifrom distributed random noise to phase For reproducibility
    % rng default;  % To generate reproducable noise
    % rng(1);  % To generate reproducable noise
    % Lower band of the noise
    Lb = -1 / 2;
    % Upper band of the noise
    Ub = 1 / 2;
    % Generate noise
    Noise = Lb + (Ub - Lb) .* rand(N, M);
    Text = append('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\dataSet', num2str(dataSetCounter), 'Noise', '.mat');
    save(Text, 'Noise');
    
    for angleCounter = 1 : 1 : angleSamplingNo
        
        NFE = 0;
        
        % Determine a new path for saving data
        PathText = append('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\noDesData', num2str(dataSetCounter), 'Theta', num2str(angleCounter), ', ', num2str(Theta(angleCounter)), ' and Phi', num2str(angleCounter), ', ', num2str(Phi(angleCounter)));
        mkdir(PathText);
        
        % Beam directions and states (u = k0 * d * sin(θ) * cos(φ), v = k0 * d * sin(θ) * sin(φ)) (Example: (0, 0) is at the center and perpendicular)
        % Column (M) --> U --> X (Axis) is x label
        % Row (N) --> V --> Y (Axis) is y label
        u0 = zeros(1, L);
        v0 = zeros(1, L);
        ControlFactor = zeros(1, L);
        for Counter = 1 : 1 : L
            u0(Counter) = k0 * d * sin(Theta(angleCounter)) * cos(Phi(angleCounter));
            v0(Counter) = k0 * d * sin(Theta(angleCounter)) * sin(Phi(angleCounter));
            ControlFactor(Counter) = 1;
        end
        % Custimize the ControlFactor to control beam's gain
        % ControlFactor(2) = 0.5;
        % ControlFactor = ones(1, L);
        
        % Save u0 and v0
        Text = append(PathText, '\u0Theta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
        save(Text, 'u0');
        Text = append(PathText, '\v0Theta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
        save(Text, 'v0');
        
        %% =============================================================================================================================================================================================
        %% Complex value of radiation pattern's samples at different values of u and v (u and v are in uv-space)
        % Pattern of the antenna (EXAMPLE). u and v are scalar values.
        % Pattern = @(u, v) sum(sinc(sx * (u - u0) / (2 * pi)) .* sinc(sy * (v - v0) / (2 * pi)));
        
        % Complex value of radiation pattern's samples at different values of u and v
%         ComplexRP = ComplexRadiationPatternOld(kValues, lValues);
%         ComplexRPSmoothed = ComplexRadiationPatternOld(kValuesSmoothed, lValuesSmoothed);
        ComplexRP = abs(ComplexRadiationPattern(2 * pi * (kValues - (M - 1) / 2) / (M - 1), 2 * pi * (lValues - (N - 1) / 2) / (N - 1)));
        ComplexRPSmoothed = abs(ComplexRadiationPattern(2 * pi * (kValuesSmoothed - (M - 1) / 2) / (M - 1), 2 * pi * (lValuesSmoothed - (N - 1) / 2) / (N - 1)));
        
        if NCRPKey == 1
            ComplexRPMax = max(ComplexRP, [], 'all');
            
            % Plot 2D and 3D complex value of radiation pattern in xy
            % Arguments are (FigureName, Xaxis, Xlabel, Yaxis, Ylabel, Zvalue, Zlabel)
            Plot2Dand3D('Sample Pattern per Element', ElemBaseX, 'Columns', ElemBaseY, 'Rows', ComplexRP, '|C(u,v)|', ComplexRPMax, [1, 1], [0, 0], 0);
            
            % Print the max values of complex value of radiation pattern
            for Counter = 1 : 1 : L
                [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshX - u0(Counter)));
                [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshY - v0(Counter)));
                disp('=====================================================================');
                disp(['Project pattern Max. direction for ' num2str(Counter) ' beam is (u, x, Column (M)): ' num2str(SizeEleX(MaxXLoc(Counter)))]);
                disp(['Project pattern Max. direction for ' num2str(Counter) ' beam is (v, y, Row (N)): ' num2str(SizeEleY(MaxYLoc(Counter)))]);
                disp(['Project pattern Max. value for ' num2str(Counter) ' beam is: ' num2str(ComplexRP(MaxYLoc(Counter), MaxXLoc(Counter)))]);
            end
        end
        
        % Print the max values of smoothed complex value of radiation pattern
        ComplexRPMaxSmoothed = 0;
        BeamMaxCRP = zeros(1, L);
        MaxXLocRange = zeros(L, DistanceFromCenter * 2 + 1);
        MaxYLocRange = zeros(L, DistanceFromCenter * 2 + 1);
        BeamMatCPR = zeros(DistanceFromCenter * 2 + 1, DistanceFromCenter * 2 + 1, L);
        CorrectionFactorCPR = zeros(1, L);
        CorrectionMatrixCRP = zeros(DistanceFromCenter * 2 + 1, DistanceFromCenter * 2 + 1, L);
        for Counter = 1 : 1 : L
            [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshXSmoothed - u0(Counter)));
            [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshYSmoothed - v0(Counter)));
            MaxXLocRange(Counter, :) = [MaxXLoc(Counter) - DistanceFromCenter : 1 : MaxXLoc(Counter) + DistanceFromCenter];
            MaxYLocRange(Counter, :) = [MaxYLoc(Counter) - DistanceFromCenter : 1 : MaxYLoc(Counter) + DistanceFromCenter];
            BeamMatCPR(:, :, Counter) = ComplexRPSmoothed(MaxYLocRange(Counter, :), MaxXLocRange(Counter, :));
            BeamMaxCRP(Counter) = max(BeamMatCPR(:, :, Counter), [], 'all');
            CorrectionFactorCPR(Counter) = 1 ./ BeamMaxCRP(Counter);
            CorrectionMatrixCRP(:, :, Counter) = CorrectionFactorCPR(Counter) .* BeamMatCPR(:, :, Counter);
            
            Text = append(PathText, '\ComplexRPMaxSmoothedTheta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
%             CRPSmoothedMAX = ComplexRPSmoothed(MaxYLoc(Counter), MaxXLoc(Counter));
            CRPSmoothedMAX = BeamMaxCRP(Counter);
            save(Text, 'CRPSmoothedMAX');
%             if CRPSmoothedMAX > ComplexRPMaxSmoothed
%                 ComplexRPMaxSmoothed = CRPSmoothedMAX;
%             end
            ComplexRPMaxSmoothed = max(ComplexRPSmoothed, [], 'all');
            
            if SCRPKey == 1
                disp('=====================================================================');
                disp(['Smoothed Project pattern Max. direction for ' num2str(Counter) ' beam is (u, x, Column (M)): ' num2str(UVBaseMeshXSmoothed(MaxXLoc(Counter)))]);
                disp(['Smoothed Project pattern Max. direction for ' num2str(Counter) ' beam is (v, y, Row (N)): ' num2str(UVBaseMeshYSmoothed(MaxYLoc(Counter)))]);
                disp(['Smoothed Project pattern Max. value for ' num2str(Counter) ' beam is: ' num2str(CRPSmoothedMAX)]);
            end
        end

%         thetaPhiPatCRNormalized(angleCounter, 2 * L + 1 : end, dataSetCounter) = reshape(ComplexRP ./ ComplexRPMax, 1, elementNo);
        thetaPhiPatCRSmoothed(angleCounter, :, dataSetCounter) = reshape(ComplexRPSmoothed, 1, elementNo * SmoothingFactor * SmoothingFactor);
        thetaPhiPatCRSmoothedNormalized(angleCounter, :, dataSetCounter) = reshape(ComplexRPSmoothed ./ ComplexRPMaxSmoothed, 1, elementNo * SmoothingFactor * SmoothingFactor);
        
        if SCRPKey == 1
            % Plot smoothed 2D and 3D complex value of radiation pattern in xy and uv space
            % Arguments are (FigureName, Xaxis, Xlabel, Yaxis, Ylabel, Zvalue, Zlabel)
            Plot2Dand3D('Sample Pattern in uv-space', gridUVBaseXSmoothed, 'u-axis', gridUVBaseYSmoothed, 'v-axis', ComplexRPSmoothed, '|C(u,v)|', ComplexRPMaxSmoothed, [1, 1], [0, 0], 0);
            
            % Smoothed complex radiation pattern in spherical coordinate
            SphericalPlot3D('Complex radiation pattern in spherical space', 'Normalized |C(x,y)|', 'CRP');
            
            % Plot of smoothed CRP in the xy-plane (Theta = 90 deg) calculated in polar system
            PolarPlot2D('Project pattern in the xy-plane (Theta = 90 deg)', 'xy-plane', 'CRP');
            
            % Plot of smoothed CRP in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
            PolarPlot2D('Project pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'CRP');
            
            % Plot of smoothed CRP in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
            PolarPlot2D('Project pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'CRP');
        end
        
        %% =============================================================================================================================================================================================
        %% Inverse Discrete Fourier Transform from complex values and pattern of the antenna
        SincPatternPhase = IDFTComplexRadiationPattern(ComplexRP);
        
        % Confine Function is confining the phase within the [-pi, pi] interval.
        Confine = @(x) pi - mod(pi - x, 2 * pi);
%         Confine = @(x) pi - rem(pi - x, 2 * pi);
        
        % Save phase without noise
        ConfinedSincPatternPhaseWithoutNoise = Confine(SincPatternPhase);
        thetaPhiPhaseWoN(angleCounter, 2 * L + 1 : end, dataSetCounter) = reshape(ConfinedSincPatternPhaseWithoutNoise, 1, elementNo);
        % Noisy phase
        NoisySincPatternPhase = SincPatternPhase + Noise;
%         NoisySincPatternPhase = SincPatternPhase + 0;
        ConfinedSincPatternPhaseNoisy = Confine(NoisySincPatternPhase);
        stochasticPhase = reshape(ConfinedSincPatternPhaseNoisy, 1, elementNo);
        % Save phase without noise
        thetaPhiPhaseWN(angleCounter, 2 * L + 1 : end, dataSetCounter) = reshape(ConfinedSincPatternPhaseNoisy, 1, elementNo);
        
        % Finds phases which are close to the boundaries
        if ~isempty(find(stochasticPhase > Criterion, elementNo))
            stochasticPhase(1, find(stochasticPhase > Criterion, elementNo)) = 0;
        end
        if ~isempty(find(stochasticPhase < -Criterion, elementNo))
            stochasticPhase(1, find(stochasticPhase < -Criterion, elementNo)) = 0;
        end
        
        % Thus, the complex phase factor (e.g., local reflection coefficient) at each of the MS elements can be written as (magnitude is 1)
        ComplexPhaseFactor = 1 * exp(1i * ConfinedSincPatternPhaseNoisy);
        
        %% =============================================================================================================================================================================================
        %% Plot phase distribution on the MS elements
        % For the ease of plotting, let us define the function phase(x, y) as a piecewise constant function
        % on the MS surface with the values of the phase at the MS elements determined by Phase(m, n),
        % actually, it converts from discrete form to continuous form.
        if PhaseKey == 1
            ElemBasePhase = @(y, x) Confine(ConfinedSincPatternPhaseNoisy(1 + round(y / d + (N - 1) / 2), 1 + round(x / d + (M - 1) / 2)));
            
            % Let us caclulate and plot the phase distribution on the MS that realises the given pattern
            ConfinedElemBasePhase2D = zeros(N, M);
            for Row = 1 : 1 : N
                for Column = 1 : 1 : M
                    ConfinedElemBasePhase2D(Row, Column) = ElemBasePhase(SizeEleY(Row), SizeEleX(Column));
                end
            end
            
            % Plot phase distribution on MS
            PlotPhase('2D Phase Distributaion', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', ConfinedElemBasePhase2D, 'Element based', ...
                                                gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
                                                SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
                                                                                                                                       [1, 1, 1]);
            
        end
        %% =============================================================================================================================================================================================
        %% Calculate Electrical Field
        % Let us calculate and plot 2D far-field radiation pattern (3D electrical field) produced by the MS. Here, we take
        % into account both phase and amplitude variation of the source antenna field on the MS, as for a spherical wave.
        % The Rad(n, m) term at denominator is presenting attenuation due to the propagation
        
        if NEFieldKey == 1
            AbsElectricField = zeros(N, M);
            for kCounter = 1 : 1 : M
                u = -pi + kValues(kCounter) * 2 * pi / (M - 1);
                % u = 2 * k0 * d * (kValues(kCounter) - (M - 1) / 2) / (M - 1);
                for lCounter = 1 : 1 : N
                    v = -pi + lValues(lCounter) * 2 * pi / (N - 1);
                    % v = 2 * k0 * d * (lValues(kCounter) - (N - 1) / 2) / (N - 1);
                    AbsElectricField(lCounter, kCounter) = ElecFieldFun(ComplexPhaseFactor, [u, v]);
                end
            end
            
            % Element-base Electrical Field
            ElemBaseElectricFieldMax = max(AbsElectricField, [], 'all');
    
            % Save E pattern
%             thetaPhiPatE(angleCounter, 2 * L + 1 : 2 * L + elementNo, dataSetCounter) = reshape(AbsElectricField, 1, elementNo);
%             thetaPhiPatENormalized(angleCounter, 2 * L + 1 : 2 * L + elementNo, dataSetCounter) = reshape(AbsElectricField ./ ElemBaseElectricFieldMax, 1, elementNo);
            
            % Plot electric field in xy coordinate
            Plot2Dand3D('2D and 3D Far-Field Radiation Pattern (Element Base)', ElemBaseX, 'Columns', ElemBaseY, 'Rows', AbsElectricField, '|E(u,v)|', ElemBaseElectricFieldMax, [1, 1], [0, 0], 0);
            
            % Print the max values of Electric Field
            for Counter = 1 : 1 : L
                [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshX - u0(Counter)));
                [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshY - v0(Counter)));
                disp('=====================================================================');
                disp(['Electric Field Max. direction for ' num2str(Counter) ' beam is (u, x, Column (M)): ' num2str(SizeEleX(MaxXLoc(Counter)))]);
                disp(['Electric Field Max. direction for ' num2str(Counter) ' beam is (v, y, Row (N)): ' num2str(SizeEleY(MaxYLoc(Counter)))]);
                disp(['Electric Field Max. value for ' num2str(Counter) ' beam is: ' num2str(ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshX(MaxXLoc(Counter)), MaxLocUVBaseMeshY(MaxYLoc(Counter))]))]);
            end
        end
        
        % Calculate Smoothed Electrical Field
        AbsElectricFieldSmoothed = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
        for kCounter = 1 : 1 : round(SmoothingFactor * M)
            u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
            % u = 2 * k0 * d * (kValuesSmoothed(kCounter) - (M - 1) / 2) / (M - 1);
            for lCounter = 1 : 1 : round(SmoothingFactor * N)
                v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
                % v = 2 * k0 * d * (lValuesSmoothed(lCounter) - (N - 1) / 2) / (N - 1);
                AbsElectricFieldSmoothed(lCounter, kCounter) = ElecFieldFun(ComplexPhaseFactor, [u, v]);
            end
        end
        
        BeamMatEF = zeros(DistanceFromCenter * 2 + 1, DistanceFromCenter * 2 + 1, L);
        BeamMaxEF = zeros(1, L);
        % Print the max values of smoothed complex value of radiation pattern
        AbsElectricFieldSmoothedMax = 0;
        for Counter = 1 : 1 : L
            [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshXSmoothed - u0(Counter)));
            [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshYSmoothed - v0(Counter)));
            BeamMatEF(:, :, Counter) = AbsElectricFieldSmoothed(MaxYLocRange(Counter, :), MaxXLocRange(Counter, :));
            BeamMaxEF(Counter) = max(BeamMatEF(:, :, Counter), [], 'all');
            
            Text = append(PathText, '\EFMaxSmoothedTheta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
%             EFSmoothedMAX = ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshYSmoothed(MaxYLoc(Counter)), MaxLocUVBaseMeshXSmoothed(MaxXLoc(Counter))]);
            EFSmoothedMAX = BeamMaxEF(Counter);
            save(Text, 'EFSmoothedMAX');
%             if EFSmoothedMAX > AbsElectricFieldSmoothedMax
%                 AbsElectricFieldSmoothedMax = EFSmoothedMAX;
%             end
            % Element-base Electrical Field
            AbsElectricFieldSmoothedMax = max(AbsElectricFieldSmoothed, [], 'all');
            
            if SEFieldKey == 1
                disp('=====================================================================');
                disp(['Smoothed Electric Field Max. direction for ' num2str(Counter) ' beam is (u, x, Column (M)): ' num2str(MaxLocUVBaseMeshXSmoothed(MaxXLoc(Counter)))]);
                disp(['Smoothed Electric Field Max. direction for ' num2str(Counter) ' beam is (v, y, Row (N)) :' num2str(MaxLocUVBaseMeshYSmoothed(MaxYLoc(Counter)))]);
                disp(['Smoothed Electric Field Max. value for ' num2str(Counter) ' beam is: ' num2str(EFSmoothedMAX)]);
            end
        end
        
        % Save E pattern smoothed
        thetaPhiPatESmoothed(angleCounter, :, dataSetCounter) = reshape(AbsElectricFieldSmoothed, 1, elementNo * SmoothingFactor * SmoothingFactor);
        thetaPhiPatESmoothedNormalized(angleCounter, :, dataSetCounter) = reshape(AbsElectricFieldSmoothed ./ AbsElectricFieldSmoothedMax, 1, elementNo * SmoothingFactor * SmoothingFactor);
        
        %% =============================================================================================================================================================================================
        %% Determine sampling and protecting regions and get samples from the input arguments
        
        % % Calculate the Jacobian transformation from (Theta, Phi) to (u, v) and determinant jacobian 
        % syms THETA PHI;
        % Vector = [THETA, PHI];
        % Fu0 = [k0 * d * sin(THETA * pi / 180) * cos(PHI * pi / 180)];
        % Fv0 = [k0 * d * sin(THETA * pi / 180) * sin(PHI * pi / 180)];
        % Func = [Fu0, Fv0];
        % Jac = jacobian(Func, Vector);
        % DetJac = det(Jac);
        % 
        % % For example for the case of (Theta, Phi) = (30, 45) in degree
        % Theta = 30;
        % Phi = 45;
        % ThetaPhiValues = [Theta, Phi];
        % JacValue = double(subs(Jac, Vector, ThetaPhiValues));
        % DetJacValue = double(subs(DetJac, Vector, ThetaPhiValues));
        
        % Apply scalling to the sinc function (Defining protected area)
        ScaleduHPBWatOrigin = uHPBWatOrigin * 2 / sx;
        ScaledvHPBWatOrigin = vHPBWatOrigin * 2 / sy;

        % Determine bigger raduis to generate sampling points inside the protected areas (Elleptic)
        if ScaleduHPBWatOrigin >= ScaledvHPBWatOrigin
            samplingRadius = ScaleduHPBWatOrigin;
        else
            samplingRadius = ScaledvHPBWatOrigin;
        end
        
        % Calculate the preotected area radius for circle scenario
%         ProtectedRadius = sqrt(power(ScaleduHPBWatOrigin, 2) + power(ScaledvHPBWatOrigin, 2));
        
        % Protection area (sqrt((u - u0) ^ 2 + (v - v0) ^ 2) < ProtectedRadius) on the UV space (circle)
        % For protection from main beams, inequality (sqrt((u - u0) ^ 2 + (v - v0) ^ 2) - ProtectedRadius) should be positive.
        % PA = @(u, v) sqrt((u - u0) ^ 2 + (v - v0) ^ 2) - ProtectedRadius; (PA = @(u, v) sqrt((u - u0) ^ 2 + (v - v0) ^ 2) - ProtectedRadius;)
        for Counter = 1 : 1 : L
            % Circle protected area
%             PA{Counter} = @(u, v) sqrt((u - u0(Counter)) ^ 2 + (v - v0(Counter)) ^ 2) - ProtectedRadius;
            % Ellips protected area
            PA{Counter} = @(u, v) sqrt((u - u0(Counter)) .^ 2 / power(ScaleduHPBWatOrigin, 2) + (v - v0(Counter)) .^ 2 / power(ScaledvHPBWatOrigin, 2)) - 1;
        end
        
        for Counter = 1 : 1 : L
        %     Xunit(1, :, Counter) = ProtectedRadius(Counter) * cos(ProtectedAreaPhi) + u0(Counter);
        %     Yunit(1, :, Counter) = ProtectedRadius(Counter) * sin(ProtectedAreaPhi) + v0(Counter);
        %     Xunit(1, :, Counter) = ProtectedRadius * sin(ProtectedAreaPhi) + v0(Counter);
        %     Yunit(1, :, Counter) = ProtectedRadius * cos(ProtectedAreaPhi) + u0(Counter);
            EFXunit(1, :, Counter) = ScaleduHPBWatOrigin * cos(ProtectedAreaPhi) + u0(Counter);
            EFYunit(1, :, Counter) = ScaledvHPBWatOrigin * sin(ProtectedAreaPhi) + v0(Counter);
        end
        
        if SEFieldKey == 1
            % Plot electric field in uv coordinate
            Plot2Dand3D('2D and 3D Far-Field Radiation Pattern (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', AbsElectricFieldSmoothed, '|E(u,v)|', AbsElectricFieldSmoothedMax, [1, 1], [1, 1], 0);
            
            % Smoothed electrical field in spherical coordinate
            SphericalPlot3D('Electrical field in spherical space', 'Normalized |E(x,y)|', 'EP');
            
            % Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
            PolarPlot2D('Electric pattern in the xy-plane (Theta = 90 deg)', 'xy-plane', 'EP');
            
            % Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
            PolarPlot2D('Electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'EP');
            
            % Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
            PolarPlot2D('Electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'EP');
        end
        
        %% Equidistant distribution of sampling points
        % Equidistant diameter [cm]
        OutsidePAEDDiameter = 0.1;
        InsidePAEDDiameter = 0.1;
%         UVBaseMeshXYInsidePA = linspace(-ProtectedRadius, ProtectedRadius, 2);
        UVBaseMeshXYInsidePA = linspace(-samplingRadius, samplingRadius, 2);
        [OutsidePAPotentialUPoints, OutsidePAPotentialVPoints] = EquiDistant2DPoints(MaxLocUVBaseMeshX, MaxLocUVBaseMeshY, OutsidePAEDDiameter);
        [InsidePAPotentialUPoints, InsidePAPotentialVPoints] = EquiDistant2DPoints(UVBaseMeshXYInsidePA, UVBaseMeshXYInsidePA, InsidePAEDDiameter);
        
        % Scattering point number outside and inside of protecting area
        ApproximateScatteringPointNumberOutsidePA = numel(OutsidePAPotentialUPoints);
        ApproximateScatteringPointNumberInsidePA = numel(InsidePAPotentialUPoints);
        
        % Number of scattered point outsice of protecting area
        OutsideProtectedAreaScatteredPointsCounter = 0;
        InsideProtectedAreaScatteredPointsCounter = zeros(1, L);
        % Location of scttered points (u, v)
        OutsideProtectedAreaScatteredPointsProperties = NaN(ApproximateScatteringPointNumberOutsidePA, 2);
        % Location of scttered points (u, v, z)
        InsideProtectedAreaScatteredPointsPropertiesWithNaN = NaN(ApproximateScatteringPointNumberInsidePA, 3, L);
        
        % Determine sampling steps
        PointCounterSamplingStep = 1;
        
        % Outside sampling of protected area
        for PointCounter = 1 : PointCounterSamplingStep : ApproximateScatteringPointNumberOutsidePA
            % Flag is used to check the scattered point are not inside of the protected area, default: the scattered point is out of PA
            Flag = 1;
            for Counter = 1 : 1 : L
                if PA{Counter}(OutsidePAPotentialUPoints(PointCounter), OutsidePAPotentialVPoints(PointCounter)) < 0
                    Flag = 0;
        %             InsideProtectedAreaScatteredPointsCounter = InsideProtectedAreaScatteredPointsCounter + 1;
        %             InsideProtectedAreaScatteredPointsProperties(InsideProtectedAreaScatteredPointsCounter, :) = [InsidePAPotentialUPoints(PointCounter), InsidePAPotentialVPoints(PointCounter), ElecFieldFun(ComplexPhaseFactor, [InsidePAPotentialUPoints(PointCounter), InsidePAPotentialVPoints(PointCounter)])];
                    break;
                end
            end
            
            if Flag == 1
                OutsideProtectedAreaScatteredPointsCounter = OutsideProtectedAreaScatteredPointsCounter + 1;
                % Position order (X, Y, Z)
%                 OutsideProtectedAreaScatteredPointsProperties(OutsideProtectedAreaScatteredPointsCounter, :) = [OutsidePAPotentialVPoints(PointCounter), OutsidePAPotentialUPoints(PointCounter), ElecFieldFun(ComplexPhaseFactor, [OutsidePAPotentialUPoints(PointCounter), OutsidePAPotentialVPoints(PointCounter)])];
                OutsideProtectedAreaScatteredPointsProperties(OutsideProtectedAreaScatteredPointsCounter, :) = [OutsidePAPotentialUPoints(PointCounter), OutsidePAPotentialVPoints(PointCounter)];
            end
        end
        
        % Getting samples from edges of the metasurface for outside of protected area
        for Ucounter = 1 : 1 : length(MaxLocUVBaseMeshXSmoothed) - 2
            UandVNeg(Ucounter, :) = [MaxLocUVBaseMeshXSmoothed(Ucounter + 1), MaxLocUVBaseMeshYSmoothed(1), ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshXSmoothed(Ucounter + 1), MaxLocUVBaseMeshYSmoothed(1)])];
            UandVPos(Ucounter, :) = [MaxLocUVBaseMeshXSmoothed(Ucounter + 1), MaxLocUVBaseMeshYSmoothed(length(MaxLocUVBaseMeshYSmoothed)), ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshXSmoothed(Ucounter + 1), MaxLocUVBaseMeshYSmoothed(length(MaxLocUVBaseMeshYSmoothed))])];
        end
        
        for Vcounter = 1 : 1 : length(MaxLocUVBaseMeshYSmoothed)
            VandUNeg(Vcounter, :) = [MaxLocUVBaseMeshXSmoothed(1), MaxLocUVBaseMeshYSmoothed(Vcounter), ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshXSmoothed(1), MaxLocUVBaseMeshYSmoothed(Vcounter)])];
            VandUPos(Vcounter, :) = [MaxLocUVBaseMeshXSmoothed(length(MaxLocUVBaseMeshXSmoothed)), MaxLocUVBaseMeshYSmoothed(Vcounter), ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshXSmoothed(length(MaxLocUVBaseMeshXSmoothed)), MaxLocUVBaseMeshYSmoothed(Vcounter)])];
        end
        
        % Sampling points values (E values) regarding fixed sampling points (location) (% Omit the NaN values from outside sampling points)
        % FixedOutsideSamplingPointsValues = OutsideProtectedAreaScatteredPointsProperties((all((~isnan(OutsideProtectedAreaScatteredPointsProperties(:, 1 : 2))), 2)), 1 : 2);
        % FixedOutsideSamplingPointsValues= rmmissing(OutsideProtectedAreaScatteredPointsProperties(:, 1 : 2));
        OutsideSamplingPoints = length([[[OutsideProtectedAreaScatteredPointsProperties((all((~isnan(OutsideProtectedAreaScatteredPointsProperties(:, 1 : 2))), 2)), 1 : 2)]; [UandVNeg(:, 1 : 2)]; [UandVPos(:, 1 : 2)]; [VandUNeg(:, 1 : 2)]; [VandUPos(:, 1 : 2)]]]);
        FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, :) = [[[OutsideProtectedAreaScatteredPointsProperties((all((~isnan(OutsideProtectedAreaScatteredPointsProperties(:, 1 : 2))), 2)), 1 : 2)]; [UandVNeg(:, 1 : 2)]; [UandVPos(:, 1 : 2)]; [VandUNeg(:, 1 : 2)]; [VandUPos(:, 1 : 2)]]];
        Text = append(PathText, '\FixedOutsideSamplingPointsValuesTheta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
        save(Text, 'FixedOutsideSamplingPointsValues');
        
        % For the scattered points inside of the protected area, a rectangular grid with equidistant points are produced and shifted to the beams directions
        InsidePAPotentialUPointsEachBeam = zeros(ApproximateScatteringPointNumberInsidePA, L);
        InsidePAPotentialVPointsEachBeam = zeros(ApproximateScatteringPointNumberInsidePA, L);
        for Counter = 1 : 1 : ApproximateScatteringPointNumberInsidePA
            for BeamCounter = 1 : 1 : L
                InsidePAPotentialUPointsEachBeam(Counter, BeamCounter) = InsidePAPotentialUPoints(Counter) + u0(BeamCounter);
                InsidePAPotentialVPointsEachBeam(Counter, BeamCounter) = InsidePAPotentialVPoints(Counter) + v0(BeamCounter);
            end
        end
        
        % Inside sampling points
        % Beam Counter
        for Counter = 1 : 1 : L
            % Inside sampling points counter for the protected area 
            for PointCounter = 1 : PointCounterSamplingStep : ApproximateScatteringPointNumberInsidePA
                if PA{Counter}(InsidePAPotentialUPointsEachBeam(PointCounter, Counter), InsidePAPotentialVPointsEachBeam(PointCounter, Counter)) <= 0
                    InsideProtectedAreaScatteredPointsCounter(1, Counter) = InsideProtectedAreaScatteredPointsCounter(1, Counter) + 1;
                    InsideProtectedAreaScatteredPointsPropertiesWithNaN(InsideProtectedAreaScatteredPointsCounter(1, Counter), :, Counter) = [InsidePAPotentialUPointsEachBeam(PointCounter, Counter), InsidePAPotentialVPointsEachBeam(PointCounter, Counter), ElecFieldFun(ComplexPhaseFactor, [InsidePAPotentialUPointsEachBeam(PointCounter, Counter), InsidePAPotentialVPointsEachBeam(PointCounter, Counter)])];
                end
            end
        end
        
        % Omit the NaN values from inside sampling points
        InsideSamplingPoints = length(rmmissing(InsideProtectedAreaScatteredPointsPropertiesWithNaN(:, :, 1)));
        InsideProtectedAreaScatteredPointsProperties = zeros(InsideSamplingPoints, 3, L);
        for Counter = 1 : 1 : L
            InsideProtectedAreaScatteredPointsProperties(:, :, Counter) = rmmissing(InsideProtectedAreaScatteredPointsPropertiesWithNaN(:, :, Counter));
        end
        clear InsideProtectedAreaScatteredPointsPropertiesWithNaN;
        
        for BeamCounter = 1 : 1 : L
            FixedInsideSamplingPointsValues(1 : InsideSamplingPoints, :, BeamCounter) = InsideProtectedAreaScatteredPointsProperties(:, 1 : 2, BeamCounter);
        end
        Text = append(PathText, '\FixedInsideSamplingPointsValuesTheta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
        save(Text, 'FixedInsideSamplingPointsValues');
        
        if SAreaKey == 1
            % Plot the protected areas and sampling points
            % Let us calculate and plot 2D and 3D far-field radiation patterns produced by the MS.
            Plot2Dand3D('2D and 3D Far-Field Radiation Pattern (UV base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', AbsElectricFieldSmoothed, '|E(u,v)|', AbsElectricFieldSmoothedMax, [1, 1], [1, 1], 1);
        end
        
        %% =============================================================================================================================================================================================
        % Optimization
        
        if OptRunKey == 1
            % Number of iteration
            evaluationNumber = elementNo * 167 + 1;
            
            % Set the initial phase values as an input of an optimization algorithm
            PhaseValuesPerIteration = zeros(N, M, evaluationNumber);
            
            Coeff1 = 1;
            Coeff2 = 10;
            CostValuePerIter = zeros(1, evaluationNumber);
            
            % Determine Boundaries (Declare Band Properties)
            % Bands Order are Lower and Upper
            Property.Lower = [];
            % Velocities Order are VelWi, VelHis and VelWis
            Property.Upper = [];
            Band = repmat(Property, 1, 1);
            clear Property;
            Band.Lower = -pi * ones(1, elementNo);
            Band.Upper = pi * ones(1, elementNo);
            
            %% Call Optimization Function
            optFlag = 0;
            Key = 0;
            OptimizationFunction = @(optimizedPhase, evaluationNumber, Band) optimizationFunction(optimizedPhase, evaluationNumber, Band);
            
            OptimizationFunction(stochasticPhase, evaluationNumber, Band);
            
            OptPhase = PhaseValuesPerIteration(:, :, NFE);
            
            % The optimized complex phase factor (e.g., local reflection coefficient) at each of the MS elements (magnitude is 1)
            global optComplexPhaseFactor;
            optComplexPhaseFactor = 1 * exp(1i * OptPhase);
            
            % Print the max values of smoothed complex value of radiation pattern
            optAbsElectricFieldSmoothedMax = 0;
            for Counter = 1 : 1 : L
                [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshXSmoothed - u0(Counter)));
                [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshYSmoothed - v0(Counter)));
                
                Text = append(PathText, '\OptEFMaxSmoothedTheta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
                OptEFSmoothedMAX = ElecFieldFun(optComplexPhaseFactor, [MaxLocUVBaseMeshXSmoothed(MaxXLoc), MaxLocUVBaseMeshYSmoothed(MaxYLoc)]);
                save(Text, 'OptEFSmoothedMAX');
                if OptEFSmoothedMAX > optAbsElectricFieldSmoothedMax
                    optAbsElectricFieldSmoothedMax = OptEFSmoothedMAX;
                end
                
                disp('=====================================================================');
                disp(['Smoothed Opt. Electric Field Max. direction for ' num2str(Counter) ' beam is (u, x, Column (M)): ' num2str(MaxLocUVBaseMeshXSmoothed(MaxXLoc(Counter)))]);
                disp(['Smoothed Opt. Electric Field Max. direction for ' num2str(Counter) ' beam is (v, y, Row (N)) :' num2str(MaxLocUVBaseMeshYSmoothed(MaxYLoc(Counter)))]);
                disp(['Smoothed Opt. Electric Field Max. value for ' num2str(Counter) ' beam is: ' num2str(OptEFSmoothedMAX)]);
            end
            
            if OptReskey == 1
                %% Plot optimized electric field in uv coordinate
                
                ElemBaseOptPhase = @(y, x) Confine(OptPhase(1 + round(y / d + (N - 1) / 2), 1 + round(x / d + (M - 1) / 2)));
                
                % Let us caclulate and plot the phase distribution on the MS that realises the given pattern
                ConfinedElemBaseOptPhase2D = zeros(N, M);
                for Row = 1 : 1 : N
                    for Column = 1 : 1 : M
                        ConfinedElemBaseOptPhase2D(Row, Column) = ElemBaseOptPhase(SizeEleY(Row), SizeEleX(Column));
                    end
                end
                
                % Determine number of iteration
                phaseIter = 100;
                
                % Assign the original phases
                modifiedPhase = ConfinedElemBaseOptPhase2D;
                
                % Post-processing on the phase distribution
                for Counter = 1 : 1 : phaseIter
                    % Take the average from phases and subtract from original ones
                    modifiedPhase = modifiedPhase - mean(modifiedPhase, 'all');
                    ConfinedPhaseModification = Confine(modifiedPhase);
                    if modifiedPhase == ConfinedPhaseModification
                        break;
                    end
                    modifiedPhase = ConfinedPhaseModification;
                end
                modifiedPhaseOpt(angleCounter, :, dataSetCounter) = reshape(modifiedPhase, 1, elementNo);
                
                % Plot phase distribution on MS
                PlotPhase('2D Optimized Phase Distributaion (Modified)', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', modifiedPhase, 'Element based', ...
                                                              gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
                                                              SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
                                                                                                                                                     [1, 1, 1]);
                
                % Plot phase distribution on MS
                PlotPhase('2D Optimized Phase Distributaion', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', ConfinedElemBaseOptPhase2D, 'Element based', ...
                                                              gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
                                                              SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
                                                                                                                                                     [1, 1, 1]);
                
                % Calculate Smoothed Electrical Field
                global optAbsElectricFieldSmoothed;
                optAbsElectricFieldSmoothed = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
                for kCounter = 1 : 1 : round(SmoothingFactor * M)
                    u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
                    %  u = 2 * k0 * d * (kValuesSmoothed(kCounter) - (M - 1) / 2) / (M - 1);
                    for lCounter = 1 : 1 : round(SmoothingFactor * N)
                        v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
                        %  v = 2 * k0 * d * (lValuesSmoothed(lCounter) - (N - 1) / 2) / (N - 1);
                        optAbsElectricFieldSmoothed(lCounter, kCounter) = ElecFieldFun(optComplexPhaseFactor, [u, v]);
                    end
                end
                
                % Element-base Electrical Field
%                 optAbsElectricFieldSmoothedMax = max(optAbsElectricFieldSmoothed, [], 'all');
                optAbsElectricFieldSmoothedMax = max(optAbsElectricFieldSmoothed(centerDistance < visibleAreaRadius));
                
                % Plot 2D and 3D Optimized Far-Field Radiation Pattern in uv-Base
                Plot2Dand3D('2D and 3D Optimized Far-Field Radiation Pattern (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', optAbsElectricFieldSmoothed, '|E_{Opt.}(u,v)|', optAbsElectricFieldSmoothedMax, [1, 1], [1, 1], 0);
                
                % Smoothed optimized electrical field in spherical coordinate
                SphericalPlot3D('Optimized electrical field in spherical space', 'Normalized |E_{Opt.}(x,y)|', 'OptEP');
                
                % Plot of optimized electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
                PolarPlot2D('Optimized electric pattern in the xy-plane (Theta = 90 deg)', 'xy-plane', 'OptEP');
                
                % Plot of optimized electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
                PolarPlot2D('Optimized electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'OptEP');
                
                % Plot of optimized electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
                PolarPlot2D('Optimized electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'OptEP');
                
                % Plot of smoothed CRP, EF, and OptEF in the xy-plane (Theta = 90 deg) calculated in polar system
                PolarPlot2D('Project pattern, EF, and OptEF in the xy-plane (Theta = 90 deg)', 'xy-plane', 'All');
                
                % Plot of smoothed CRP, EF, and OptEF in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
                PolarPlot2D('Optimized electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'All');
                
                % Plot of smoothed CRP, EF, and OptEF in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
                PolarPlot2D('Optimized electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'All');
            end
            
            % Save data
            Text = append(PathText, '\NFETheta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
            save(Text, 'NFE');
            Text = append(PathText, '\PhaseValuesPerIterationTheta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
            save(Text, 'PhaseValuesPerIteration');
            Text = append(PathText, '\CostValuePerIterTheta', num2str(angleCounter), 'Phi', num2str(angleCounter), '.mat');
            save(Text, 'CostValuePerIter');
            thetaPhiPhaseOpt(angleCounter, 2 * L + 1 : end, dataSetCounter) = reshape(OptPhase, 1, elementNo);
            %  thetaPhiPatEOpt(angleCounter, 2 * L + 1 : end, dataSetCounter) = reshape(PatternE(:, :, NFE), 1, elementNo);
            %  thetaPhiPatEOptNormalized(angleCounter, 2 * L + 1 : end, dataSetCounter) = reshape(PatternENormalized(:, :, NFE), 1, elementNo);
            thetaPhiPatESmoothedOpt(angleCounter, :, dataSetCounter) = reshape(PatternESmooth(:, :, NFE), 1, elementNo * SmoothingFactor * SmoothingFactor);
            thetaPhiPatESmoothedOptNormalized(angleCounter, :, dataSetCounter) = reshape(PatternESmoothedNormalized(:, :, NFE), 1, elementNo * SmoothingFactor * SmoothingFactor);
        end
    end
end

save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPhaseWoN.mat', 'thetaPhiPhaseWoN');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPhaseWN.mat', 'thetaPhiPhaseWN');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPhaseOpt.mat', 'thetaPhiPhaseOpt');
if OptReskey == 1
    save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\modifiedPhaseOpt.mat', 'modifiedPhaseOpt');
end
% save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatCR.mat', 'thetaPhiPatCR');
% save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatCRNormalized.mat', 'thetaPhiPatCRNormalized');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatCRSmoothed.mat', 'thetaPhiPatCRSmoothed');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatCRSmoothedNormalized.mat', 'thetaPhiPatCRSmoothedNormalized');


% save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatE.mat', 'thetaPhiPatE');
% save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatENormalized.mat', 'thetaPhiPatENormalized');
% save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatEOpt.mat', 'thetaPhiPatEOpt');
% save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatEOptNormalized.mat', 'thetaPhiPatEOptNormalized');

save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatESmoothed.mat', 'thetaPhiPatESmoothed');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatESmoothedNormalized.mat', 'thetaPhiPatESmoothedNormalized');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatESmoothedOpt.mat', 'thetaPhiPatESmoothedOpt');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\thetaPhiPatESmoothedOptNormalized.mat', 'thetaPhiPatESmoothedOptNormalized');

save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Coeff1.mat', 'Coeff1');
save('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Coeff2.mat', 'Coeff2');
