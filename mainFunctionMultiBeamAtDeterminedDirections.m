clear all;
close all;
clc;

%% Enable and presents the mentioned sections
% When the key is equal to 1, it will represent the results
% When the key is equal to 0, it will not represent the results

% Plot not smoothed Complex Radiation Pattern (not smoothed)
NCRPKey = 0;

% Plot smoothed Complex Radiation Pattern
SCRPKey = 1;

% Plot phase distribution
PhaseKey = 1;

% Plot not smoothed electric field
NEFieldKey = 1;

% Plot smoothed electric field
SEFieldKey = 1;

% Plot sampling area
SAreaKey = 1;

% Run optmization
OptRunKey = 0;

% Plot optimized results
OptReskey = 1;

%% In this example calculation, the time dependence is exp(+i * ω * t), as is generally assumed in electrical engineering.
% Therefore, when an EM wave propagates along a distance Rad in free space it acquires a phase delay given by the phase factor exp(-i * k0 * Rad),
% where k0 = ω / c = 2 * pi / Lambda (wave number or angular wavenumber).

% Below, Rf is the focal distance, which is the distance from the source antenna to the MS.
% The MS is formed by M * N elements (reflecting or transmitting, which does not matter for this calculation).
% Let us define the MS parameters: Rf, M, N, d
global Rf;
Rf = 50.000000;
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\Rf.mat', 'Rf');
% Number of elements in each Rows
global M;
M = 10;
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\M.mat', 'M');
% Number of elements in each Columns
global N;
N = 3;
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\N.mat', 'N');

% Operational frequency
Frequency = 10e9;

% Determine radius of Protected Areas (PAs)
% Increasing or decreasing Q's value increases or decreases the radiuses of PAs.
Q = 1 / 4;

% Let us define the sinc(x) function as a pattern
% Let us now define some radiation pattern we want to realize. In particular, it makes sense to consider patterns
% given by a superposition of sinc(sx * (u - u0) / 2) * sinc(sy * (v - v0) / 2) terms, as these are the basic patterns
% of a rectangular aperture with the size sx * d by sy * d. Thus for the radiation pattern with L beams we have
% Let us define the radiation pattern parameters with L beams:
global L;
L = 2;
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\L.mat', 'L');

% Number of angles for constructing beams in arbitrary directions
% Such as combination of Theta = [50, 25, 0, 25, 50] and Phi = [90, 90, 0, -90, -90] ==> {[Theta1, Phi1], [Theta2, Phi2], [Theta3, Phi3], [Theta4, Phi4], [Theta5, Phi5]}
% ThetaPhiBeams = {[50, 90], [25, 90], [0, 0], [25, -90], [50, -90]};
ThetaPhiBeams = {[50, 0], [25, 0], [0, 0], [25, 180], [50, 180]};

% Rad(m, n) is the radial distance from the source antenna to the (m, n)-th element of the MS.
% d [cm] is the period of the MS (the distance between the centers of its elements), which is the same along x and y axes.
global d;
d = 0.68;
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\d.mat', 'd');

% m and n are scalar values which presents mth and nth element of the metasurface array
% Rad = @ (m, n) sqrt(Rf^2 + ((m - (M - 1) / 2) ^ 2 + (n - (N - 1) / 2) ^ 2) * d ^ 2);

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
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\Lambda.mat', 'Lambda');
k0 = 2 * pi / Lambda;

% Number of cost function evaluation (NFE)
global NFE optFlag Key;

global u0 v0 ControlFactor;

global SmoothingFactor kValues lValues kValuesSmoothed lValuesSmoothed UVBaseMeshYSmoothed UVBaseMeshXSmoothed UVBaseXSmoothed UVBaseYSmoothed;
% global ComplexRPSmoothed;
global CorrectionFactorCPR;
global AbsElectricFieldSmoothedMax;
% global BeamMaxEF;
% global PatternE PatternENormalized;
global PatternESmooth PatternESmoothedNormalized;

global PhaseValuesPerIteration;
% global AveBeamMaxEF;
global Coeff1 Coeff2;
global CostValuePerIter;

global FixedOutsideSamplingPointsValues FixedInsideSamplingPointsValues;
global OutsideSamplingPoints InsideSamplingPoints;
% FixedOutsideSamplingPointsValues = zeros(5000, 2);
% FixedInsideSamplingPointsValues = zeros(200, 2, L);

% Generate Theta and Phi samples which is distributed equidistantly over the upper part of spherical
elementNo = M * N;

% Number of angle sets (nCr = n! / r! * (n - r)!)
angleSamplingNo = factorial(length(ThetaPhiBeams)) / (factorial(L) * factorial(length(ThetaPhiBeams) - L));
% Determine beam combinations, each row includes Theta and Phi for each beam. For example for L = 2, we have [[Theta1, Phi1], [Theta2, Phi2]] in each row
BeamCombinations = cell2mat(nchoosek(ThetaPhiBeams, L)) * pi / 180;
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\BeamCombinations.mat', 'BeamCombinations');

Theta = zeros(L * angleSamplingNo, 1);
Phi = zeros(L * angleSamplingNo, 1);
ColumnCounter = 1;
for Counter = 1 : 1 : L
    Theta((Counter - 1) * angleSamplingNo + 1 : Counter * angleSamplingNo) = BeamCombinations(:, ColumnCounter);
    Phi((Counter - 1) * angleSamplingNo + 1 : Counter * angleSamplingNo) = BeamCombinations(:, ColumnCounter + 1);
    ColumnCounter = ColumnCounter + 2;
end
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\Theta.mat', 'Theta');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\Phi.mat', 'Phi');

% Number of desired data sets
noDesData = 10;

% Smoothing Factor to make patterns more smoother
SmoothingFactor = 5;

% Data
% Theta, Phi, M * N Phases (without noise)
thetaPhiPhaseWoN = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, M * N Phases (with noise)
thetaPhiPhaseWN = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, M * N Phases (optimized)
thetaPhiPhaseOpt = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, M * N CRP
% thetaPhiPatCR = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, M * N CRP (normalized)
% thetaPhiPatCRNormalized = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) CRP
thetaPhiPatCRSmoothed = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) CRP (normalized)
thetaPhiPatCRSmoothedNormalized = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);

% Theta, Phi, M * N E-pattern
% thetaPhiPatE = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, M * N E-pattern (normalized)
% thetaPhiPatENormalized = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, M * N E-patter (optimized)
% thetaPhiPatEOpt = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, M * N E-patter (normalized and optimized)
% thetaPhiPatEOptNormalized = zeros(angleSamplingNo, L * 2 + elementNo, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) E-pattern
thetaPhiPatESmoothed = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) E-pattern (normalized)
thetaPhiPatESmoothedNormalized = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) E-pattern (optimized)
thetaPhiPatESmoothedOpt = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);
% Theta, Phi, (M * SmoothingFactor) * (N * SmoothingFactor) E-pattern (normalized and optimized)
thetaPhiPatESmoothedOptNormalized = zeros(angleSamplingNo, elementNo * SmoothingFactor * SmoothingFactor, noDesData);

% Get values around the u0 and v0 to from main beams amplitude
DistanceFromCenter = 7;

% Mesh regarding Massive MIMO Antenna rows and columns indices (Element Base)
% Column
ElemBaseMeshX = linspace(-floor(N / 2), floor(N / 2), N);
ElementsX = linspace(1, N, N);
% Row
ElemBaseMeshY = linspace(-floor(M / 2), floor(M / 2), M);
ElementsY = linspace(1, M, M);
% Generate the mesh grid
[ElemBaseX, ElemBaseY] = meshgrid(ElemBaseMeshX, ElemBaseMeshY);
[gridElemBaseX, gridElemBaseY] = meshgrid(ElementsX, ElementsY);

% k and l are number of elements in each Row and Column, respectively.
kValues = linspace(0, M - 1, M);
lValues = linspace(0, N - 1, N);
% Convert coordinate to UV space
% Column
UVBaseMeshX = 2 * k0 * d * (lValues - (N - 1) / 2) / (N - 1);
% Row
UVBaseMeshY = 2 * k0 * d * (kValues - (M - 1) / 2) / (M - 1);
% Generate mesh grid
[UVBaseX, UVBaseY] = meshgrid(UVBaseMeshX, UVBaseMeshY);

% Smoothed mesh regarding Massive MIMO Antenna
kValuesSmoothed = linspace(0, M - 1, round(SmoothingFactor * M));
lValuesSmoothed = linspace(0, N - 1, round(SmoothingFactor * N));

% Convert coordinate to UV space
% Column
UVBaseMeshXSmoothed = 2 * k0 * d * (lValuesSmoothed - (N - 1) / 2) / (N - 1);
% Row
UVBaseMeshYSmoothed = 2 * k0 * d * (kValuesSmoothed - (M - 1) / 2) / (M - 1);
% Generate the mesh grid
[UVBaseXSmoothed, UVBaseYSmoothed] = meshgrid(UVBaseMeshXSmoothed, UVBaseMeshYSmoothed);

% For finding exact loacation of beams and plotting E-field in uv-space
MaxLocUVBaseMeshX = 2 * pi * (lValues - (N - 1) / 2) / (N - 1);
MaxLocUVBaseMeshY = 2 * pi * (kValues - (M - 1) / 2) / (M - 1);
% Generate the mesh grid for finding exact loacation of beams and plotting E-field in uv-space
[MaxLocUVBaseX, MaxLocUVBaseY] = meshgrid(MaxLocUVBaseMeshX, MaxLocUVBaseMeshY);

% For plotting E-field in uv-space and finding exact loacation of beams
MaxLocUVBaseMeshXSmoothed = 2 * pi * (lValuesSmoothed - (N - 1) / 2) / (N - 1);
MaxLocUVBaseMeshYSmoothed = 2 * pi * (kValuesSmoothed - (M - 1) / 2) / (M - 1);
% Generate mesh grid for plotting E-field in uv-space and finding exact loacation of beams
[MaxLocUVBaseXSmoothed, MaxLocUVBaseYSmoothed] = meshgrid(MaxLocUVBaseMeshXSmoothed, MaxLocUVBaseMeshYSmoothed);

% Generate mesh grid for plotting in spherical coordinate using Phi and Theta values
gridThetaSmoothed = linspace(0 , pi / 2, SmoothingFactor * M);
gridPhiSmoothed = linspace(0 , 2 * pi, SmoothingFactor * N);
plotUVBaseMeshXSmoothed = zeros(SmoothingFactor * M, SmoothingFactor * N);
plotUVBaseMeshYSmoothed = zeros(SmoothingFactor * M, SmoothingFactor * N);
for thetaCounter = 1 : 1 : SmoothingFactor * M
    for phiCounter = 1 : 1 : SmoothingFactor * N
        plotUVBaseMeshXSmoothed(thetaCounter, phiCounter) = k0 * d * sin(gridThetaSmoothed(thetaCounter)) * cos(gridPhiSmoothed(phiCounter));
        plotUVBaseMeshYSmoothed(thetaCounter, phiCounter) = k0 * d * sin(gridThetaSmoothed(thetaCounter)) * sin(gridPhiSmoothed(phiCounter));
    end
end

% Generate grid to plot in real dimension on the MS
SizeEleX = linspace(-(N - 1) * d / 2, (N - 1) * d / 2, N);
SizeEleY = linspace(-(M - 1) * d / 2, (M - 1) * d / 2, M);

% Criterian of the optimization algorithm to change some of phase values
global Criterion;
Criterion = 3;

% Center of the visible area
U0 = 0;
V0 = 0;
% Radius of visible area in UV-space
visibleAreaRadius = k0 * d;
% Define Theta to plot circle for visible area in UV-space
visibleAreaTheta = 0 : pi / 100 : 2 * pi;
% Generate grid for visible area in UV-space
visibleAreaX = visibleAreaRadius .* cos(visibleAreaTheta) + V0;
visibleAreaY = visibleAreaRadius .* sin(visibleAreaTheta) + U0;

% Distance matrix for determining visible area
Distance = sqrt(power(MaxLocUVBaseXSmoothed, 2) + power(MaxLocUVBaseYSmoothed, 2));

% Determine which elements are inside visible area
[InVisAreaX, InVisAreaY] = find(Distance <= visibleAreaRadius);

% Protected area (ovals)
ProtectedAreaPhi = 0 : pi / 100 : 2 * pi;

% Spherical plot of the Complex Radiation Pattern
% Generate two vectors for Theta and Phi
SmoothedMeshThetaSphCoor = linspace(0, pi / 2, SmoothingFactor * M);
SmoothedMeshPhiSpeCoor = linspace(0, 2 * pi, SmoothingFactor * N);

for kCounter = 1 : 1 : round(SmoothingFactor * M)
    for lCounter = 1 : 1 : round(SmoothingFactor * N)
        % Generating mesh in spherical coordinate
        % Number of rows
        Ymn(kCounter, lCounter) = sin(SmoothedMeshThetaSphCoor(kCounter)) * cos(SmoothedMeshPhiSpeCoor(lCounter));
        % Number of columns
        Xmn(kCounter, lCounter) = sin(SmoothedMeshThetaSphCoor(kCounter)) * sin(SmoothedMeshPhiSpeCoor(lCounter));
    end
end
Zmn = real(sqrt(1 - power(Xmn, 2) - power(Ymn, 2)));

% Generate vector of Theta and Phi for electric pattern in the planes
SmoothedMeshThetaPlane = linspace(-pi / 2, pi / 2, SmoothingFactor * elementNo);
SmoothedMeshPhiPlane = linspace(0, 2 * pi, SmoothingFactor * elementNo);

% xz-plane (phi = 90 or 270 deg, therefore u = 0)
phiValueXZ = 90 * pi / 180;
% xy-plane (theta = 90 deg)
ThetaValue = 90 * pi / 180;
% yz-plane (phi = 0 or 180 deg, therefore v = 0)
phiValueYZ = 0 * pi / 180;

xzSphCorAbsElectricFieldSmoothed = zeros(1, SmoothingFactor * elementNo);
xySphCorAbsElectricFieldSmoothed = zeros(1, SmoothingFactor * elementNo);
yzSphCorAbsElectricFieldSmoothed = zeros(1, SmoothingFactor * elementNo);

for Counter = 1 : 1 : SmoothingFactor * elementNo
    % Number of rows (XZ)
    xzYmn(Counter) = sin(SmoothedMeshThetaPlane(Counter)) * cos(phiValueXZ);
    % Number of columns (XZ)
    xzXmn(Counter) = sin(SmoothedMeshThetaPlane(Counter)) * sin(phiValueXZ);

    % Number of rows (XY)
    xyYmn(Counter) = sin(ThetaValue) * cos(SmoothedMeshPhiPlane(Counter));
    % Number of columns (XY)
    xyXmn(Counter) = sin(ThetaValue) * sin(SmoothedMeshPhiPlane(Counter));

    % Number of rows (YZ)
    yzYmn(Counter) = sin(SmoothedMeshThetaPlane(Counter)) * cos(phiValueYZ);
    % Number of columns (YZ)
    yzXmn(Counter) = sin(SmoothedMeshThetaPlane(Counter)) * sin(phiValueYZ);
end
xzZmn = real(sqrt(1 - power(xzXmn, 2) - power(xzYmn, 2)));
xyZmn = real(sqrt(1 - power(xyXmn, 2) - power(xyYmn, 2)));
yzZmn = real(sqrt(1 - power(yzXmn, 2) - power(yzYmn, 2)));

% Set axes color for figures
left_color = [0 0 0];
right_color = [0 0 0];

% Save Theta and Phi
thetaPhiPhaseWoN(:, 1 : 2 * L, :) = repmat(BeamCombinations, [1 1 noDesData]);
thetaPhiPhaseWN(:, 1 : 2 * L, :) = repmat(BeamCombinations, [1 1 noDesData]);
thetaPhiPhaseOpt(:, 1 : 2 * L, :) = repmat(BeamCombinations, [1 1 noDesData]);

for dataSetCounter = 1 : 1 : noDesData
    
    % Add Zero Mean Unifrom distributed random noise to phase For reproducibility
    % rng default;  % To generate reproducable noise
    % rng(1);  % To generate reproducable noise
    % Lower band of the noise
    Lb = -1 / 2;
    % Upper band of the noise
    Ub = 1 / 2;
    % Generate noise
    Noise = Lb + (Ub - Lb) .* rand(M, N);
    Text = append('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\dataSet', num2str(dataSetCounter), 'Noise', '.mat');
    save(Text, 'Noise');
    
    for angleCounter = 1 : 1 : angleSamplingNo
        
        NFE = 0;
        
        % Determine a new path for saving data
        PathText = append('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\noDesData', num2str(dataSetCounter), 'ThetaAndPhiPairNo', num2str(angleCounter));
        mkdir(PathText);
        
        % Beam directions and states (u = k0 * d * sin(θ) * cos(φ), v = k0 * d * sin(θ) * sin(φ)) (Example: (0, 0) is at the center and perpendicular)
        % Column (N) --> V --> X (Axis) is x label and ylabel('Row (M) --> U --> Y (Axis) is y label
        u0 = zeros(1, L);
        v0 = zeros(1, L);
        ControlFactor = zeros(1, L);
        ColumnCounter = 1;
        for Counter = 1 : 1 : L
            u0(Counter) = k0 * d * sin(BeamCombinations(angleCounter, ColumnCounter)) * cos(BeamCombinations(angleCounter, ColumnCounter + 1));
            v0(Counter) = k0 * d * sin(BeamCombinations(angleCounter, ColumnCounter)) * sin(BeamCombinations(angleCounter, ColumnCounter + 1));
            ControlFactor(Counter) = 1;
            ColumnCounter = ColumnCounter + 2;
        end
        % Custimize the ControlFactor to control beam's gain
        % ControlFactor(2) = 0.5;
        
        % Save u0 and v0
        Text = append(PathText, '\u0ThetaAndPhiPairNo', num2str(angleCounter), '.mat');
        save(Text, 'u0');
        Text = append(PathText, '\v0ThetaAndPhiPairNo', num2str(angleCounter), '.mat');
        save(Text, 'v0');
        
        %% =============================================================================================================================================================================================
        %% Complex value of radiation pattern's samples at different values of u and v (u and v are in Spherical Coordinate)
        % Pattern of the antenna (EXAMPLE). u and v are scalar values.
        % Pattern = @(u, v) sum(sinc(sx * (u - u0) / (2 * pi)) .* sinc(sy * (v - v0) / (2 * pi)));
        
        % Complex value of radiation pattern's samples at different values of u and v
%         ComplexRP = ComplexRadiationPatternOld(kValues, lValues);
%         ComplexRPSmoothed = ComplexRadiationPatternOld(kValuesSmoothed, lValuesSmoothed);
        ComplexRP = abs(ComplexRadiationPattern(2 * pi * (kValues - (M - 1) / 2) / (M - 1), 2 * pi * (lValues - (N - 1) / 2) / (N - 1)));
        ComplexRPSmoothed = abs(ComplexRadiationPattern(2 * pi * (kValuesSmoothed - (M - 1) / 2) / (M - 1), 2 * pi * (lValuesSmoothed - (N - 1) / 2) / (N - 1)));
        % Save ComplexRPMax and Save ComplexRPMaxSmoothed
%         thetaPhiPatCR(angleCounter, L * 2 + 1 : end, dataSetCounter) = reshape(ComplexRP, 1, elementNo);
        
        if NCRPKey == 1
            ComplexRPMax = max(ComplexRP, [], 'all');

            % Plot 2D and 3D complex value of radiation pattern in xy
            Fig = figure('Name', 'Sample Pattern in xy coordinate', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            subplot(1, 2, 1);
            ElemBaseSurf2DC = surf(ElemBaseX, ElemBaseY, ComplexRP);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('Columns', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('Rows', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            title('Sample Pattern per Element', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            set(ElemBaseSurf2DC, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            view(2);
            subplot(1, 2, 2);
            ElemBaseSurf3DC = surf(ElemBaseX, ElemBaseY, ComplexRP / ComplexRPMax);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('Columns', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('Rows', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            title('Normalized Sample Pattern per Element', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            set(ElemBaseSurf3DC, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            
            % Print the max values of complex value of radiation pattern
            for Counter = 1 : 1 : L
                [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshX - v0(Counter)));
                [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshY - u0(Counter)));
                disp('=====================================================================');
                disp(['CRP Max. direction for ' num2str(Counter) ' beam is (u, y, Row (M)): ' num2str(SizeEleY(MaxYLoc(Counter)))]);
                disp(['CRP Max. direction for ' num2str(Counter) ' beam is (v, x, Column (N)): ' num2str(SizeEleX(MaxXLoc(Counter)))]);
                disp(['CRP Max. value for ' num2str(Counter) ' beam is: ' num2str(ComplexRP(MaxYLoc(Counter), MaxXLoc(Counter)))]);
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
            [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshXSmoothed - v0(Counter)));
            [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshYSmoothed - u0(Counter)));
            MaxXLocRange(Counter, :) = [MaxXLoc(Counter) - DistanceFromCenter : 1 : MaxXLoc(Counter) + DistanceFromCenter];
            MaxYLocRange(Counter, :) = [MaxYLoc(Counter) - DistanceFromCenter : 1 : MaxYLoc(Counter) + DistanceFromCenter];
            BeamMatCPR(:, :, Counter) = ComplexRPSmoothed(MaxYLocRange(Counter, :), MaxXLocRange(Counter, :));
            BeamMaxCRP(Counter) = max(BeamMatCPR(:, :, Counter), [], 'all');
            
            Text = append(PathText, '\ComplexRPMaxSmoothedThetaAndPhiPairNo', num2str(angleCounter), 'L', num2str(Counter), '.mat');
%             CRPSmoothedMAX = ComplexRPSmoothed(MaxYLoc(Counter), MaxXLoc(Counter));
            CRPSmoothedMAX = BeamMaxCRP(Counter);
            save(Text, 'CRPSmoothedMAX');
%             if CRPSmoothedMAX > ComplexRPMaxSmoothed
%                 ComplexRPMaxSmoothed = CRPSmoothedMAX;
%             end
            ComplexRPMaxSmoothed = max(ComplexRPSmoothed, [], 'all');
            
            CorrectionFactorCPR(Counter) = 1 ./ BeamMaxCRP(Counter);
            CorrectionMatrixCRP(:, :, Counter) = CorrectionFactorCPR(Counter) .* BeamMatCPR(:, :, Counter);
            if SCRPKey == 1
                disp('=====================================================================');
                disp(['Smoothed CRP Max. direction for ' num2str(Counter) ' beam is (u, y, Row (M)): ' num2str(UVBaseMeshYSmoothed(MaxYLoc(Counter)))]);
                disp(['Smoothed CRP Max. direction for ' num2str(Counter) ' beam is (v, x, Column (N)): ' num2str(UVBaseMeshXSmoothed(MaxXLoc(Counter)))]);
                disp(['Smoothed CRP Max. value for ' num2str(Counter) ' beam is: ' num2str(CRPSmoothedMAX)]);
            end
        end
        
%         thetaPhiPatCRNormalized(angleCounter, L * 2 + 1 : end, dataSetCounter) = reshape(ComplexRP ./ ComplexRPMax, 1, elementNo);
        thetaPhiPatCRSmoothed(angleCounter, :, dataSetCounter) = reshape(ComplexRPSmoothed, 1, elementNo * SmoothingFactor * SmoothingFactor);
        thetaPhiPatCRSmoothedNormalized(angleCounter, :, dataSetCounter) = reshape(ComplexRPSmoothed ./ ComplexRPMaxSmoothed, 1, elementNo * SmoothingFactor * SmoothingFactor);
        
        if SCRPKey == 1
            % Plot Smoothed 2D and 3D complex value of radiation pattern in xy and uv space
            Fig = figure('Name', 'Sample Pattern in uv space', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            subplot(1, 2, 1);
            UVBaseSurf2DC = surf(UVBaseXSmoothed, UVBaseYSmoothed, ComplexRPSmoothed);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            title('Sample Pattern in UV-Space (2D)', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            set(UVBaseSurf2DC, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            view(2);
            subplot(1, 2, 2);
            UVBaseSurf3D = surf(UVBaseXSmoothed, UVBaseYSmoothed, ComplexRPSmoothed / ComplexRPMaxSmoothed);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            title('Normalized Sample Pattern in UV-Space (3D)', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            set(UVBaseSurf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            % Enhancing the surface characteristics
            axis tight;

            % Spherical plot of the complex radiation pattern
            % Smoothed complex radiation pattern in spherical coordinate
            SphCorAbsCRPSmoothed = zeros(round(SmoothingFactor * M), round(SmoothingFactor * N));
            for kCounter = 1 : 1 : round(SmoothingFactor * M)
                for lCounter = 1 : 1 : round(SmoothingFactor * N)
                    % Calculate u and v values
                    uSphCor = k0 * d * Ymn(kCounter, lCounter);
                    vSphCor = k0 * d * Xmn(kCounter, lCounter);
                    
                    % Calculate complex radiation pattern
                    SphCorAbsCRPSmoothed(kCounter, lCounter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
                end
            end
            Zmn = real(sqrt(1 - power(Xmn, 2) - power(Ymn, 2)));
            
            % Scale Xmn, Ymn, and Zmn
            ScXmn = Xmn .* SphCorAbsCRPSmoothed;
            ScYmn = Ymn .* SphCorAbsCRPSmoothed;
            ScZmn = Zmn .* SphCorAbsCRPSmoothed;
            
            % Plot complex radiation pattern in spherical coordinate
            Fig = figure('Name', 'Complex radiation pattern in spherical space', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            sphericalSurf3D = surf(ScXmn, ScYmn, ScZmn);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('x', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('y', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            zlabel('|Complex Radiation Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            % title('Sampled complex radiation pattern in spherical space', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            set(sphericalSurf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            % Enhancing the surface characteristics
            axis tight;
            axis equal;

            % Plot of smoothed CRP in the xz-plane (phi = 90 or 270 deg, therefore u = 0)
            for Counter = 1 : 1 : SmoothingFactor * elementNo
                % Calculate u and v values
                uSphCor = k0 * d * xzYmn(Counter);
                vSphCor = k0 * d * xzXmn(Counter);
                % Calculate electrical field
                xzSphCorAbsElectricFieldSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
            end
            % Scale Xmn, Ymn, and Zmn
            xzScXmn = xzXmn .* xzSphCorAbsElectricFieldSmoothed;
            xzScYmn = xzYmn .* xzSphCorAbsElectricFieldSmoothed;
            xzScZmn = xzZmn .* xzSphCorAbsElectricFieldSmoothed;
            % Plot of electric pattern in the xz-plane
            Fig = figure('Name', 'CRP in the xz-plane (phi = 90 or 270 deg, therefore u = 0)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            plot(xzScXmn, xzScZmn, 'LineWidth', 1.5, 'Color', 'black');
            grid on;
            % Setting the axes
            xlabel('X', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
            ylabel('Z', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
            axis equal;
            axis([-inf inf -inf inf]);
            
            % Plot of smoothed CRP in the xy-plane (theta = 90 deg)
            for Counter = 1 : 1 : SmoothingFactor * elementNo
                % Calculate u and v values
                uSphCor = k0 * d * xyYmn(Counter);
                vSphCor = k0 * d * xyXmn(Counter);
                % Calculate electrical field
                xySphCorAbsElectricFieldSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
            end
            % Scale Xmn, Ymn, and Zmn
            xyScXmn = xyXmn .* xySphCorAbsElectricFieldSmoothed;
            xyScYmn = xyYmn .* xySphCorAbsElectricFieldSmoothed;
            xyScZmn = xyZmn .* xySphCorAbsElectricFieldSmoothed;
            % Plot of electric pattern in the xy-plane
            Fig = figure('Name', 'CRP in the xy-plane (theta = 90 deg)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            plot(xyScXmn, xyScYmn, 'LineWidth', 1.5, 'Color', 'black');
            grid on;
            % Setting the axes
            xlabel('X', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
            ylabel('Y', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
            axis equal;
            axis([-inf inf -inf inf]);
            
            % Plot of smoothed CRP in the yz-plane (phi = 0 or 180 deg, therefore v = 0)
            for Counter = 1 : 1 : SmoothingFactor * elementNo
                % Calculate u and v values
                uSphCor = k0 * d * yzYmn(Counter);
                vSphCor = k0 * d * yzXmn(Counter);
                % Calculate electrical field
                yzSphCorAbsElectricFieldSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
            end
            % Scale Xmn, Ymn, and Zmn
            yzScXmn = yzXmn .* yzSphCorAbsElectricFieldSmoothed;
            yzScYmn = yzYmn .* yzSphCorAbsElectricFieldSmoothed;
            yzScZmn = yzZmn .* yzSphCorAbsElectricFieldSmoothed;
            % Plot of electric pattern in the yz-plane
            Fig = figure('Name', 'CRP in the yz-plane (phi = 0 or 180 deg, therefore v = 0)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            plot(yzScYmn, yzScZmn, 'LineWidth', 1.5, 'Color', 'black');
            grid on;
            % Setting the axes
            xlabel('Y', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
            ylabel('Z', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
            axis equal;
            axis([-inf inf -inf inf]);
        end
        
        %% =============================================================================================================================================================================================
        %% Inverse Discrete Fourier Transform from complex values and pattern of the antenna
        SincPatternPhase = IDFTComplexRadiationPattern(ComplexRP);
        
        % Confine Function is confining the phase within the [-pi, pi] interval.
        Confine = @(x) pi - mod(pi - x, 2 * pi);
%         Confine = @(x) pi - rem(pi - x, 2 * pi);
        
        % Save phase without noise
        ConfinedSincPatternPhaseWithoutNoise = Confine(SincPatternPhase);
        thetaPhiPhaseWoN(angleCounter, L * 2 + 1 : end, dataSetCounter) = reshape(ConfinedSincPatternPhaseWithoutNoise, 1, elementNo);
        % Noisy phase
        NoisySincPatternPhase = SincPatternPhase + Noise;
%         NoisySincPatternPhase = SincPatternPhase + 0;
        ConfinedSincPatternPhaseNoisy = Confine(NoisySincPatternPhase);
        stochasticPhase = reshape(ConfinedSincPatternPhaseNoisy, 1, elementNo);
        % Save phase without noise
        thetaPhiPhaseWN(angleCounter, L * 2 + 1 : end, dataSetCounter) = reshape(ConfinedSincPatternPhaseNoisy, 1, elementNo);
        
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
            ElemBasePhase = @(x, y) Confine(ConfinedSincPatternPhaseNoisy(1 + round(x / d + (M - 1) / 2), 1 + round(y / d + (N - 1) / 2)));
            
            % Let us caclulate and plot the phase distribution on the MS that realises the given pattern
            ConfinedElemBasePhase2D = zeros(M, N);
            for Row = 1 : 1 : M
                for Column = 1 : 1 : N
                    ConfinedElemBasePhase2D(Row, Column) = ElemBasePhase(SizeEleY(Row), SizeEleX(Column));
                end
            end
            
            Fig = figure('Name', '2D Phase Distributaion', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            subplot(3, 1, 1);
            ElementBaseSurf2DPhaseDis = surf(gridElemBaseX, gridElemBaseY, ConfinedElemBasePhase2D);
            title('Element based');
            xlim([-inf inf]);
            xticks([1 : 1 : N]);
            ylim([-inf inf]);
            yticks([1 : 1 : M]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('Columns', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel('Rows', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            set(ElementBaseSurf2DPhaseDis, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            axis equal;
            view(2);
            subplot(3, 1, 2);
            Surf2D = surf(UVBaseX, UVBaseY, ConfinedElemBasePhase2D);
            title('UV-space');
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            axis equal;
            view(2);
            subplot(3, 1, 3);
            Surf2D = surf(SizeEleX, SizeEleY, ConfinedElemBasePhase2D);
            title('Physical dimension');
            xlim([-inf inf]);
            xticks([SizeEleX]);
            ylim([-inf inf]);
            yticks([SizeEleY]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('Column Dimension (x)', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel('Row Dimension (y)', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            axis equal;
            view(2);
    
%             Fig = figure('Name', '2D Phase Distributaion (Element based)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
%             ElementBaseSurf2DPhaseDis = surf(gridElemBaseX, gridElemBaseY, ConfinedElemBasePhase2D);
%             xlim([-inf inf]);
%             xticks([1 : 1 : N]);
%             ylim([-inf inf]);
%             yticks([1 : 1 : M]);
%             colorbar;
%             colormap(winter(256));
%             % Setting the axes label, type, position, rotation, etc.
%             xlabel('Columns', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             ylabel('Rows', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             % Enhancing the surface characteristics
%             axis tight;
%             set(ElementBaseSurf2DPhaseDis, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
%             axis equal;
%             view(2);
%             
%             Fig = figure('Name', '2D Phase Distributaion (UV-space)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
%             Surf2D = surf(UVBaseX, UVBaseY, ConfinedElemBasePhase2D);
%             xlim([-inf inf]);
%             ylim([-inf inf]);
%             colorbar;
%             colormap(winter(256));
%             % Setting the axes label, type, position, rotation, etc.
%             xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             % Enhancing the surface characteristics
%             axis tight;
%             set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
%             axis equal;
%             view(2);
%             
%             Fig = figure('Name', '2D Phase Distributaion (Physical dimension)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
%             Surf2D = surf(SizeEleX, SizeEleY, ConfinedElemBasePhase2D);
%             xlim([-inf inf]);
%             xticks([SizeEleX]);
%             ylim([-inf inf]);
%             yticks([SizeEleY]);
%             colorbar;
%             colormap(winter(256));
%             % Setting the axes label, type, position, rotation, etc.
%             xlabel('Column Dimension (x)', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             ylabel('Row Dimension (y)', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             % Enhancing the surface characteristics
%             axis tight;
%             set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
%             axis equal;
%             view(2);

        end
        %% =============================================================================================================================================================================================
        %% Calculate Electrical Field
        % Let us calculate and plot 2D far-field radiation pattern (3D electrical field) produced by the MS. Here, we take
        % into account both phase and amplitude variation of the source antenna field on the MS, as for a spherical wave.
        % The Rad(m, n) term at denominator is presenting attenuation due to the propagation

        if NEFieldKey == 1
            AbsElectricField = zeros(M, N);
            for kCounter = 1 : 1 : M
                u = -pi + kValues(kCounter) * 2 * pi / (M - 1);
                % u = 2 * k0 * d * (kValues(kCounter) - (M - 1) / 2) / (M - 1);
                for lCounter = 1 : 1 : N
                    v = -pi + lValues(lCounter) * 2 * pi / (N - 1);
                    % v = 2 * k0 * d * (lValues(kCounter) - (N - 1) / 2) / (N - 1);
                    AbsElectricField(kCounter, lCounter) = ElecFieldFun(ComplexPhaseFactor, [u, v]);
                end
            end
            
            % Element-base Electrical Field
            ElemBaseElectricFieldMax = max(AbsElectricField, [], 'all');
    
            % Save E pattern
            thetaPhiPatE(angleCounter, L * 2 + 1 : L * 2 + elementNo, dataSetCounter) = reshape(AbsElectricField, 1, elementNo);
            thetaPhiPatENormalized(angleCounter, L * 2 + 1 : L * 2 + elementNo, dataSetCounter) = reshape(AbsElectricField ./ ElemBaseElectricFieldMax, 1, elementNo);
            
            % Plot electric field in xy coordinate
            Fig = figure('Name', '2D and 3D Far-Field Radiation Pattern (Element Base)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            subplot(1, 2, 1);
            ElemBase2DElecField = surf(ElemBaseX, ElemBaseY, AbsElectricField);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('Columns', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('Rows', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            zlabel('|Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            set(ElemBase2DElecField, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            view(2);
            
            subplot(1, 2, 2);
            ElemBase3DElecField = surf(ElemBaseX, ElemBaseY, AbsElectricField / ElemBaseElectricFieldMax);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('Columns', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('Rows', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            zlabel('Nor. |Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            set(ElemBase3DElecField, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            
            % Print the max values of Electric Field
            for Counter = 1 : 1 : L
                [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshX - v0(Counter)));
                [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshY - u0(Counter)));
                disp('=====================================================================');
                disp(['Electric Field Max. direction for ' num2str(Counter) ' beam is (u, y, Row (M)): ' num2str(SizeEleY(MaxYLoc(Counter)))]);
                disp(['Electric Field Max. direction for ' num2str(Counter) ' beam is (v, x, Column (N)): ' num2str(SizeEleX(MaxXLoc(Counter)))]);
                disp(['Electric Field Max. value for ' num2str(Counter) ' beam is: ' num2str(ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshY(MaxYLoc(Counter)), MaxLocUVBaseMeshX(MaxXLoc(Counter))]))]);
            end
        end
        
        % Calculate Smoothed Electrical Field
        AbsElectricFieldSmoothed = zeros(round(SmoothingFactor * M), round(SmoothingFactor * N));
        for kCounter = 1 : 1 : round(SmoothingFactor * M)
            u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
            % u = 2 * k0 * d * (kValuesSmoothed(kCounter) - (M - 1) / 2) / (M - 1);
            for lCounter = 1 : 1 : round(SmoothingFactor * N)
                v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
                % v = 2 * k0 * d * (lValuesSmoothed(lCounter) - (N - 1) / 2) / (N - 1);
                AbsElectricFieldSmoothed(kCounter, lCounter) = ElecFieldFun(ComplexPhaseFactor, [u, v]);
            end
        end

        BeamMatEF = zeros(DistanceFromCenter * 2 + 1, DistanceFromCenter * 2 + 1, L);
        BeamMaxEF = zeros(1, L);
        % Print the max values of smoothed complex value of radiation pattern
        AbsElectricFieldSmoothedMax = 0;
        for Counter = 1 : 1 : L
            [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshXSmoothed - v0(Counter)));
            [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshYSmoothed - u0(Counter)));
            BeamMatEF(:, :, Counter) = AbsElectricFieldSmoothed(MaxYLocRange(Counter, :), MaxXLocRange(Counter, :));
            BeamMaxEF(Counter) = max(BeamMatEF(:, :, Counter), [], 'all');
            
            Text = append(PathText, '\EFMaxSmoothedThetaAndPhiPairNo', num2str(angleCounter), 'L', num2str(Counter), '.mat');
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
                disp(['Smoothed Electric Field Max. direction for ' num2str(Counter) ' beam is (u, y, Row (M)): ' num2str(MaxLocUVBaseMeshYSmoothed(MaxYLoc(Counter)))]);
                disp(['Smoothed Electric Field Max. direction for ' num2str(Counter) ' beam is (v, x, Column (N)) :' num2str(MaxLocUVBaseMeshXSmoothed(MaxXLoc(Counter)))]);
                disp(['Smoothed Electric Field Max. value for ' num2str(Counter) ' beam is: ' num2str(EFSmoothedMAX)]);
            end
        end
        
        % Save E pattern smoothed
        thetaPhiPatESmoothed(angleCounter, :, dataSetCounter) = reshape(AbsElectricFieldSmoothed, 1, elementNo * SmoothingFactor * SmoothingFactor);
        thetaPhiPatESmoothedNormalized(angleCounter, :, dataSetCounter) = reshape(AbsElectricFieldSmoothed ./ AbsElectricFieldSmoothedMax, 1, elementNo * SmoothingFactor * SmoothingFactor);
        
        if SEFieldKey == 1
            % Plot electric field in uv coordinate
            Fig = figure('Name', '2D and 3D Far-Field Radiation Pattern (UV Base)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            subplot(1, 2, 1);
            UVBase2DElecField = surf(MaxLocUVBaseXSmoothed, MaxLocUVBaseYSmoothed, AbsElectricFieldSmoothed);
            hold on;
            % Plot visible area in UV-space
            plot3(visibleAreaX, visibleAreaY, repmat(max(AbsElectricFieldSmoothed(:, :), [], 'all'), 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            zlabel('|Electric Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            camlight;
            lighting phong;
            shading interp;
            set(UVBase2DElecField, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            hold off;
            view(2);
            subplot(1, 2, 2);
            UVBase3DElecField = surf(MaxLocUVBaseXSmoothed, MaxLocUVBaseYSmoothed, AbsElectricFieldSmoothed / AbsElectricFieldSmoothedMax, 'DisplayName', 'Nor. |Electric Pattern|');
            hold on;
            Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2, 'DisplayName', 'Visible Area');
            % Plot visible area in UV-space
            for PatternValue = 0.1 : 0.1 : max(AbsElectricFieldSmoothed(:, :), [], 'all') / AbsElectricFieldSmoothedMax
                plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2);
            end
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            zlabel('Nor. |Electric Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            camlight;
            lighting phong;
            shading interp;
            set(UVBase3DElecField, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            legend([UVBase3DElecField Leg2], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
            hold off;
            
%             % Spherical plot of the Electric Field
%             % Smoothed electrical field in spherical coordinate
%             SphCorAbsElectricFieldSmoothed = zeros(round(SmoothingFactor * M), round(SmoothingFactor * N));
%             for kCounter = 1 : 1 : round(SmoothingFactor * M)
%                 for lCounter = 1 : 1 : round(SmoothingFactor * N)
%                     % Calculate u and v values
%                     uSphCor = k0 * d * Ymn(kCounter, lCounter);
%                     vSphCor = k0 * d * Xmn(kCounter, lCounter);
%                     % Calculate electrical field
%                     SphCorAbsElectricFieldSmoothed(kCounter, lCounter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor]);
%                 end
%             end
%             % Scale Xmn, Ymn, and Zmn
%             ScXmn = Xmn .* SphCorAbsElectricFieldSmoothed;
%             ScYmn = Ymn .* SphCorAbsElectricFieldSmoothed;
%             ScZmn = Zmn .* SphCorAbsElectricFieldSmoothed;
%             % Plot electric field in spherical coordinate
%             Fig = figure('Name', 'Electrical field in spherical space', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
%             sphericalSurf3D = surf(ScXmn, ScYmn, ScZmn);
%             xlim([-inf inf]);
%             ylim([-inf inf]);
%             colorbar;
%             colormap(winter(256));
%             % Setting the axes label, type, position, rotation, etc.
%             xlabel('x', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             ylabel('y', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             zlabel('|Electric Field|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             % title('Electric field pattern in spherical space', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%             set(sphericalSurf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
%             % Enhancing the surface characteristics
%             axis tight;
%             axis equal;
%             
%             % Plot of electric pattern in the xz-plane (phi = 90 or 270 deg, therefore u = 0)
%             for Counter = 1 : 1 : SmoothingFactor * elementNo
%                 % Calculate u and v values
%                 uSphCor = k0 * d * xzYmn(Counter);
%                 vSphCor = k0 * d * xzXmn(Counter);
%                 % Calculate electrical field
%                 xzSphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor]);
%             end
%             % Scale Xmn, Ymn, and Zmn
%             xzScXmn = xzXmn .* xzSphCorAbsElectricFieldSmoothed;
%             xzScYmn = xzYmn .* xzSphCorAbsElectricFieldSmoothed;
%             xzScZmn = xzZmn .* xzSphCorAbsElectricFieldSmoothed;
%             % Plot of electric pattern in the xz-plane
%             Fig = figure('Name', 'Electric pattern in the xz-plane (phi = 90 or 270 deg, therefore u = 0)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
%             plot(xzScXmn, xzScZmn, 'LineWidth', 1.5, 'Color', 'black');
%             grid on;
%             % Setting the axes
%             xlabel('X', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
%             ylabel('Z', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
%             axis equal;
%             axis([-inf inf -inf inf]);
%             
%             % Plot of electric pattern in the xy-plane (theta = 90 deg)
%             for Counter = 1 : 1 : SmoothingFactor * elementNo
%                 % Calculate u and v values
%                 uSphCor = k0 * d * xyYmn(Counter);
%                 vSphCor = k0 * d * xyXmn(Counter);
%                 % Calculate electrical field
%                 xySphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor]);
%             end
%             % Scale Xmn, Ymn, and Zmn
%             xyScXmn = xyXmn .* xySphCorAbsElectricFieldSmoothed;
%             xyScYmn = xyYmn .* xySphCorAbsElectricFieldSmoothed;
%             xyScZmn = xyZmn .* xySphCorAbsElectricFieldSmoothed;
%             % Plot of electric pattern in the xy-plane
%             Fig = figure('Name', 'Electric pattern in the xy-plane (theta = 90 deg)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
%             plot(xyScXmn, xyScYmn, 'LineWidth', 1.5, 'Color', 'black');
%             grid on;
%             % Setting the axes
%             xlabel('X', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
%             ylabel('Y', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
%             axis equal;
%             axis([-inf inf -inf inf]);
%             
%             % Plot of electric pattern in the yz-plane (phi = 0 or 180 deg, therefore v = 0)
%             for Counter = 1 : 1 : SmoothingFactor * elementNo
%                 % Calculate u and v values
%                 uSphCor = k0 * d * yzYmn(Counter);
%                 vSphCor = k0 * d * yzXmn(Counter);
%                 % Calculate electrical field
%                 yzSphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor]);
%             end
%             % Scale Xmn, Ymn, and Zmn
%             yzScXmn = yzXmn .* yzSphCorAbsElectricFieldSmoothed;
%             yzScYmn = yzYmn .* yzSphCorAbsElectricFieldSmoothed;
%             yzScZmn = yzZmn .* yzSphCorAbsElectricFieldSmoothed;
%             % Plot of electric pattern in the yz-plane
%             Fig = figure('Name', 'Electric pattern in the yz-plane (phi = 0 or 180 deg, therefore v = 0)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
%             plot(yzScYmn, yzScZmn, 'LineWidth', 1.5, 'Color', 'black');
%             grid on;
%             % Setting the axes
%             xlabel('Y', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
%             ylabel('Z', 'FontWeight', 'bold', 'FontSize' , 18, 'FontName' , 'Times New Roman');
%             axis equal;
%             axis([-inf inf -inf inf]);
        end
        
        
        
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
        
        % Calculate the HPBW or 1 / 10 (is equal to -20 dB) for (u, v)
        syms uvSymbole;
        % uEqu = sinc(sx * (uvSymbole - u0) / 2) - 1 / 2 == 0;
        % uEqu = sinc(sx * (uvSymbole - u0) / 2) - 1 / 10 == 0;
        uEqu = sinc(uvSymbole) - Q == 0;
        uHPBWatOrigin = abs(vpasolve(uEqu, uvSymbole));
        
        % vEqu = sinc(sy * (uvSymbole - v0) / 2) - 1 / 2 == 0;
        % vEqu = sinc(sy * (uvSymbole - v0) / 2) - 1 / 10 == 0;
        vEqu = sinc(uvSymbole) - Q == 0;
        vHPBWatOrigin = abs(vpasolve(vEqu, uvSymbole));
        
        % Apply scalling to the sinc function
        ScaleduHPBWatOrigin = uHPBWatOrigin * 2 / sx;
        ScaledvHPBWatOrigin = vHPBWatOrigin * 2 / sy;

        % Determine bigger raduis to generate sampling points inside the protected areas
        if ScaleduHPBWatOrigin >= ScaledvHPBWatOrigin
            samplingRadius = ScaleduHPBWatOrigin;
        else
            samplingRadius = ScaledvHPBWatOrigin;
        end
        
        % Calculate the preotected area radius
%         ProtectedRadius = sqrt(power(ScaleduHPBWatOrigin, 2) + power(ScaledvHPBWatOrigin, 2));
        
        % Protection area (sqrt((u - u0) ^ 2 + (v - v0) ^ 2) < ProtectedRadius) on the UV space
        % For protection from main beams, inequality (sqrt((u - u0) ^ 2 + (v - v0) ^ 2) - ProtectedRadius) should be positive.
        % PA = @(u, v) sqrt((u - u0) ^ 2 + (v - v0) ^ 2) - ProtectedRadius; (PA = @(u, v) sqrt((u - u0) ^ 2 + (v - v0) ^ 2) - ProtectedRadius;)
        for Counter = 1 : 1 : L
%             PA{Counter} = @(u, v) sqrt((u - u0(Counter)) ^ 2 + (v - v0(Counter)) ^ 2) - ProtectedRadius;
            PA{Counter} = @(u, v) sqrt((u - u0(Counter)) ^ 2 / power(ScaleduHPBWatOrigin, 2) + (v - v0(Counter)) ^ 2 / power(ScaledvHPBWatOrigin, 2)) - 1;
        end
        
        Xunit = zeros(1, length(ProtectedAreaPhi), L);
        Yunit = zeros(1, length(ProtectedAreaPhi), L);
        for Counter = 1 : 1 : L
        %     Xunit(1, :, Counter) = ProtectedRadius(Counter) * cos(ProtectedAreaPhi) + u0(Counter);
        %     Yunit(1, :, Counter) = ProtectedRadius(Counter) * sin(ProtectedAreaPhi) + v0(Counter);
        %     Xunit(1, :, Counter) = ProtectedRadius * sin(ProtectedAreaPhi) + v0(Counter);
        %     Yunit(1, :, Counter) = ProtectedRadius * cos(ProtectedAreaPhi) + u0(Counter);
            
            Xunit(1, :, Counter) = ScaledvHPBWatOrigin * sin(ProtectedAreaPhi) + v0(Counter);
            Yunit(1, :, Counter) = ScaleduHPBWatOrigin * cos(ProtectedAreaPhi) + u0(Counter);
        end
        
        %% Equidistant distribution of sampling points
        % Equidistant diameter [cm]
        OutsidePAEDDiameter = 0.1;
        InsidePAEDDiameter = 0.1;
%         UVBaseMeshXYInsidePA = linspace(-ProtectedRadius, ProtectedRadius, 2);
        UVBaseMeshXYInsidePA = linspace(-samplingRadius, samplingRadius, 2);
        [OutsidePAPotentialUPoints, OutsidePAPotentialVPoints] = EquiDistant2DPoints(MaxLocUVBaseMeshY, MaxLocUVBaseMeshX, OutsidePAEDDiameter);
        [InsidePAPotentialUPoints, InsidePAPotentialVPoints] = EquiDistant2DPoints(UVBaseMeshXYInsidePA, UVBaseMeshXYInsidePA, InsidePAEDDiameter);
        
        % Scattering point number outside and inside of protecting area
        ApproximateScatteringPointNumberOutsidePA = numel(OutsidePAPotentialUPoints);
        ApproximateScatteringPointNumberInsidePA = numel(InsidePAPotentialUPoints);
        
        % Number of scattered point outsice of protecting area
        OutsideProtectedAreaScatteredPointsCounter = 0;
        InsideProtectedAreaScatteredPointsCounter = zeros(1, L);
        % Location of scttered points (u, v, z)
%         OutsideProtectedAreaScatteredPointsProperties = NaN(ApproximateScatteringPointNumberOutsidePA, 3);
        OutsideProtectedAreaScatteredPointsProperties = NaN(ApproximateScatteringPointNumberOutsidePA, 2);
        % Location of scttered points (u, v, Index, Index, z)
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
                OutsideProtectedAreaScatteredPointsProperties(OutsideProtectedAreaScatteredPointsCounter, :) = [OutsidePAPotentialVPoints(PointCounter), OutsidePAPotentialUPoints(PointCounter)];
            end
        end
        
        % Getting samples from edges of the metasurface for outside of protected area
        for Ucounter = 1 : 1 : length(MaxLocUVBaseMeshYSmoothed)
            UandVNeg(Ucounter, :) = [MaxLocUVBaseMeshXSmoothed(1), MaxLocUVBaseMeshYSmoothed(Ucounter), ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshYSmoothed(Ucounter), MaxLocUVBaseMeshXSmoothed(1)])];
            UandVPos(Ucounter, :) = [MaxLocUVBaseMeshXSmoothed(length(MaxLocUVBaseMeshXSmoothed)), MaxLocUVBaseMeshYSmoothed(Ucounter), ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshYSmoothed(Ucounter), MaxLocUVBaseMeshXSmoothed(length(MaxLocUVBaseMeshXSmoothed))])];
        end
        
        for Vcounter = 1 : 1 : length(MaxLocUVBaseMeshXSmoothed) - 2
            VandUNeg(Vcounter, :) = [MaxLocUVBaseMeshXSmoothed(Vcounter + 1), MaxLocUVBaseMeshYSmoothed(1), ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshYSmoothed(1), MaxLocUVBaseMeshXSmoothed(Vcounter + 1)])];
            VandUPos(Vcounter, :) = [MaxLocUVBaseMeshXSmoothed(Vcounter + 1), MaxLocUVBaseMeshYSmoothed(length(MaxLocUVBaseMeshYSmoothed)), ElecFieldFun(ComplexPhaseFactor, [MaxLocUVBaseMeshYSmoothed(length(MaxLocUVBaseMeshYSmoothed)), MaxLocUVBaseMeshXSmoothed(Vcounter + 1)])];
        end
        
        % Sampling points values (E values) regarding fixed sampling points (location) (% Omit the NaN values from outside sampling points)
        % FixedOutsideSamplingPointsValues = OutsideProtectedAreaScatteredPointsProperties((all((~isnan(OutsideProtectedAreaScatteredPointsProperties(:, 1 : 2))), 2)), 1 : 2);
        % FixedOutsideSamplingPointsValues= rmmissing(OutsideProtectedAreaScatteredPointsProperties(:, 1 : 2));
        OutsideSamplingPoints = length([[[OutsideProtectedAreaScatteredPointsProperties((all((~isnan(OutsideProtectedAreaScatteredPointsProperties(:, 1 : 2))), 2)), 1 : 2)]; [UandVNeg(:, 1 : 2)]; [UandVPos(:, 1 : 2)]; [VandUNeg(:, 1 : 2)]; [VandUPos(:, 1 : 2)]]]);
        FixedOutsideSamplingPointsValues(1 : OutsideSamplingPoints, :) = [[[OutsideProtectedAreaScatteredPointsProperties((all((~isnan(OutsideProtectedAreaScatteredPointsProperties(:, 1 : 2))), 2)), 1 : 2)]; [UandVNeg(:, 1 : 2)]; [UandVPos(:, 1 : 2)]; [VandUNeg(:, 1 : 2)]; [VandUPos(:, 1 : 2)]]];
        Text = append(PathText, '\FixedOutsideSamplingPointsValuesThetaAndPhiPairNo', num2str(angleCounter), '.mat');
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
                    InsideProtectedAreaScatteredPointsPropertiesWithNaN(InsideProtectedAreaScatteredPointsCounter(1, Counter), :, Counter) = [InsidePAPotentialVPointsEachBeam(PointCounter, Counter), InsidePAPotentialUPointsEachBeam(PointCounter, Counter), ElecFieldFun(ComplexPhaseFactor, [InsidePAPotentialUPointsEachBeam(PointCounter, Counter), InsidePAPotentialVPointsEachBeam(PointCounter, Counter)])];
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
        Text = append(PathText, '\FixedInsideSamplingPointsValuesThetaAndPhiPairNo', num2str(angleCounter), '.mat');
        save(Text, 'FixedInsideSamplingPointsValues');
        
        if SAreaKey == 1
            % Plot the protected areas and sampling points
            % Let us calculate and plot 2D and 3D far-field radiation patterns produced by the MS.
            Fig = figure('Name', '2D and 3D Far-Field Radiation Pattern (UV base)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
            subplot(1, 2, 1);
            UVBase2DElecFieldProArea = surf(MaxLocUVBaseXSmoothed, MaxLocUVBaseYSmoothed, AbsElectricFieldSmoothed);
            hold on;
            % plot3(Xunit, Yunit, repmat(max(UVCubicE(:, :), [], 'all'), 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
            for Counter = 1 : 1 : L
                % Plot protected area
                plot3(Xunit(1, :, Counter), Yunit(1, :, Counter), repmat(max(AbsElectricFieldSmoothed(:, :), [], 'all'), 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
            end
            % Plot visible area in UV-space
            plot3(visibleAreaX, visibleAreaY, repmat(max(AbsElectricFieldSmoothed(:, :), [], 'all'), 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2);
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            zlabel('|Electric Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            camlight;
            lighting phong;
            shading interp;
            set(UVBase2DElecFieldProArea, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 2), FixedOutsideSamplingPointsValues(:, 1)]), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4]);
            for Counter = 1: 1 : L
                scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1]);
            end
            hold off;
            view(2);
            subplot(1, 2, 2);
            UVBase3DElecFieldProArea = surf(MaxLocUVBaseXSmoothed, MaxLocUVBaseYSmoothed, AbsElectricFieldSmoothed / AbsElectricFieldSmoothedMax, 'DisplayName', 'Nor. |Electric Pattern|');
            hold on;
            Leg2 = plot3(Xunit(1, :, 1), Yunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
            Leg3 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2, 'DisplayName', 'Visible Area');
            % Plot visible area in UV-space
            for PatternValue = 0.1 : 0.1 : max(AbsElectricFieldSmoothed(:, :), [], 'all') / AbsElectricFieldSmoothedMax
                plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2);
            end
            for Counter = 1 : 1 : L
                % Plot protected area
                for PatternValue = 0 : 0.1 : max(AbsElectricFieldSmoothed(:, :), [], 'all') / AbsElectricFieldSmoothedMax
                    plot3(Xunit(1, :, Counter), Yunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                end
            end
            xlim([-inf inf]);
            ylim([-inf inf]);
            colorbar;
            colormap(winter(256));
            % Setting the axes label, type, position, rotation, etc.
            xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            zlabel('Nor. |Electrical Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            camlight;
            lighting phong;
            shading interp;
            set(UVBase3DElecFieldProArea, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
            Leg4 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 2), FixedOutsideSamplingPointsValues(:, 1)]) / AbsElectricFieldSmoothedMax, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
            for Counter = 1 : 1 : L
                Leg5 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / AbsElectricFieldSmoothedMax, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
            end
            legend([UVBase3DElecFieldProArea Leg2 Leg3 Leg4 Leg5], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
            hold off;
        end
        
        %% =============================================================================================================================================================================================
        % Optimization
        
        if OptRunKey == 1
            % Number of iteration
            evaluationNumber = elementNo * 167 + 1;
            
            % Set the initial phase values as an input of an optimization algorithm
            PhaseValuesPerIteration = zeros(M, N, evaluationNumber);
            
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
            optComplexPhaseFactor = 1 * exp(1i * OptPhase);
            
            

            % Print the max values of smoothed complex value of radiation pattern
            optAbsElectricFieldSmoothedMax = 0;
            for Counter = 1 : 1 : L
                [ValX(Counter), MaxXLoc(Counter)] = min(abs(MaxLocUVBaseMeshXSmoothed - v0(Counter)));
                [ValY(Counter), MaxYLoc(Counter)] = min(abs(MaxLocUVBaseMeshYSmoothed - u0(Counter)));
                
                Text = append(PathText, '\OptEFMaxSmoothedThetaAndPhiPairNo', num2str(angleCounter), 'L', num2str(Counter), '.mat');
                OptEFSmoothedMAX = ElecFieldFun(optComplexPhaseFactor, [MaxLocUVBaseMeshYSmoothed(MaxYLoc), MaxLocUVBaseMeshXSmoothed(MaxXLoc)]);
                save(Text, 'OptEFSmoothedMAX');
                if OptEFSmoothedMAX > optAbsElectricFieldSmoothedMax
                    optAbsElectricFieldSmoothedMax = OptEFSmoothedMAX;
                end
                
                disp('=====================================================================');
                disp(['Smoothed Opt. Electric Field Max. direction for ' num2str(Counter) ' beam is (u, y, Row (M)): ' num2str(MaxLocUVBaseMeshYSmoothed(MaxYLoc(Counter)))]);
                disp(['Smoothed Opt. Electric Field Max. direction for ' num2str(Counter) ' beam is (v, x, Column (N)) :' num2str(MaxLocUVBaseMeshXSmoothed(MaxXLoc(Counter)))]);
                disp(['Smoothed Opt. Electric Field Max. value for ' num2str(Counter) ' beam is: ' num2str(OptEFSmoothedMAX)]);
            end
            
            if OptReskey == 1
                %% Plot optimized electric field in uv coordinate
                
                ElemBaseOptPhase = @(x, y) Confine(OptPhase(1 + round(x / d + (M - 1) / 2), 1 + round(y / d + (N - 1) / 2)));
                
                % Let us caclulate and plot the phase distribution on the MS that realises the given pattern
                ConfinedElemBaseOptPhase2D = zeros(M, N);
                for Row = 1 : 1 : M
                    for Column = 1 : 1 : N
                        ConfinedElemBaseOptPhase2D(Row, Column) = ElemBaseOptPhase(SizeEleY(Row), SizeEleX(Column));
                    end
                end
                
                Fig = figure('Name', '2D Optimized Phase Distributaion', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(3, 1, 1);
                ElementBaseSurf2DPhaseDis = surf(gridElemBaseX, gridElemBaseY, ConfinedElemBaseOptPhase2D);
                title('Element based');
                xlim([-inf inf]);
                xticks([1 : 1 : N]);
                ylim([-inf inf]);
                yticks([1 : 1 : M]);
                colorbar;
                colormap(winter(256));
                % Setting the axes label, type, position, rotation, etc.
                xlabel('Columns', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
                ylabel('Rows', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
                % Enhancing the surface characteristics
                axis tight;
                set(ElementBaseSurf2DPhaseDis, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                axis equal;
                view(2);
                subplot(3, 1, 2);
                Surf2D = surf(UVBaseX, UVBaseY, ConfinedElemBaseOptPhase2D);
                title('UV-space');
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(winter(256));
                % Setting the axes label, type, position, rotation, etc.
                xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
                ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
                % Enhancing the surface characteristics
                axis tight;
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                axis equal;
                view(2);
                subplot(3, 1, 3);
                Surf2D = surf(SizeEleX, SizeEleY, ConfinedElemBaseOptPhase2D);
                title('Physical dimension');
                xlim([-inf inf]);
                xticks([SizeEleX]);
                ylim([-inf inf]);
                yticks([SizeEleY]);
                colorbar;
                colormap(winter(256));
                % Setting the axes label, type, position, rotation, etc.
                xlabel('Column Dimension (x)', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
                ylabel('Row Dimension (y)', 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
                % Enhancing the surface characteristics
                axis tight;
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                axis equal;
                view(2);
                
                % Calculate Smoothed Electrical Field
                optAbsElectricFieldSmoothed = zeros(round(SmoothingFactor * M), round(SmoothingFactor * N));
                for kCounter = 1 : 1 : round(SmoothingFactor * M)
                    u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
                    %  u = 2 * k0 * d * (kValuesSmoothed(kCounter) - (M - 1) / 2) / (M - 1);
                    for lCounter = 1 : 1 : round(SmoothingFactor * N)
                        v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
                        %  v = 2 * k0 * d * (lValuesSmoothed(lCounter) - (N - 1) / 2) / (N - 1);
                        optAbsElectricFieldSmoothed(kCounter, lCounter) = ElecFieldFun(optComplexPhaseFactor, [u, v]);
                    end
                end
                
                % Element-base Electrical Field
%                 optAbsElectricFieldSmoothedMax = max(optAbsElectricFieldSmoothed, [], 'all');
                
                Fig = figure('Name', '2D and 3D Optimized Far-Field Radiation Pattern (UV Base)', 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                UVBase2DElecField = surf(MaxLocUVBaseXSmoothed, MaxLocUVBaseYSmoothed, optAbsElectricFieldSmoothed);
                hold on;
                % Plot visible area in UV-space
                plot3(visibleAreaX, visibleAreaY, repmat(max(optAbsElectricFieldSmoothed(InVisAreaX, InVisAreaY), [], 'all'), 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2);
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(winter(256));
                % Setting the axes label, type, position, rotation, etc.
                xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel('|Electric Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                set(UVBase2DElecField, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                hold off;
                view(2);
                subplot(1, 2, 2);
                UVBase3DElecField = surf(MaxLocUVBaseXSmoothed, MaxLocUVBaseYSmoothed, optAbsElectricFieldSmoothed / max(optAbsElectricFieldSmoothed(InVisAreaX, InVisAreaY), [], 'all'), 'DisplayName', 'Nor. |Electric Pattern|');
                hold on;
                Leg2 = plot3(visibleAreaX(1, :, 1), visibleAreaY(1, :, 1), repmat(0, 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2, 'DisplayName', 'Visible Area');
                % Plot visible area in UV-space
                for PatternValue = 0.1 : 0.1 : max(optAbsElectricFieldSmoothed(:, :), [], 'all') / max(optAbsElectricFieldSmoothed(InVisAreaX, InVisAreaY), [], 'all')
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'cyan', 'LineWidth', 2);
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(winter(256));
                % Setting the axes label, type, position, rotation, etc.
                xlabel('v-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel('u-axis', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel('Nor. |Electric Pattern|', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                set(UVBase3DElecField, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([UVBase3DElecField Leg2], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                
            end
            
            
            
            % Save data
            Text = append(PathText, '\NFEThetaAndPhiPairNo', num2str(angleCounter), '.mat');
            save(Text, 'NFE');
            Text = append(PathText, '\PhaseValuesPerIterationThetaAndPhiPairNo', num2str(angleCounter), '.mat');
            save(Text, 'PhaseValuesPerIteration');
            Text = append(PathText, '\CostValuePerIterThetaAndPhiPairNo', num2str(angleCounter), '.mat');
            save(Text, 'CostValuePerIter');
            thetaPhiPhaseOpt(angleCounter, L * 2 + 1 : end, dataSetCounter) = reshape(OptPhase, 1, elementNo);
            %  thetaPhiPatEOpt(angleCounter, L * 2 + 1 : end, dataSetCounter) = reshape(PatternE(:, :, NFE), 1, elementNo);
            %  thetaPhiPatEOptNormalized(angleCounter, L * 2 + 1 : end, dataSetCounter) = reshape(PatternENormalized(:, :, NFE), 1, elementNo);
            thetaPhiPatESmoothedOpt(angleCounter, :, dataSetCounter) = reshape(PatternESmooth(:, :, NFE), 1, elementNo * SmoothingFactor * SmoothingFactor);
            thetaPhiPatESmoothedOptNormalized(angleCounter, :, dataSetCounter) = reshape(PatternESmoothedNormalized(:, :, NFE), 1, elementNo * SmoothingFactor * SmoothingFactor);

            PlotVisibleAndProtectedCRP

        end
    end
end

save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPhaseWoN.mat', 'thetaPhiPhaseWoN');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPhaseWN.mat', 'thetaPhiPhaseWN');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPhaseOpt.mat', 'thetaPhiPhaseOpt');

% save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatCR.mat', 'thetaPhiPatCR');
% save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatCRNormalized.mat', 'thetaPhiPatCRNormalized');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatCRSmoothed.mat', 'thetaPhiPatCRSmoothed');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatCRSmoothedNormalized.mat', 'thetaPhiPatCRSmoothedNormalized');


% save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatE.mat', 'thetaPhiPatE');
% save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatENormalized.mat', 'thetaPhiPatENormalized');
% save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatEOpt.mat', 'thetaPhiPatEOpt');
% save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatEOptNormalized.mat', 'thetaPhiPatEOptNormalized');

save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatESmoothed.mat', 'thetaPhiPatESmoothed');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatESmoothedNormalized.mat', 'thetaPhiPatESmoothedNormalized');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatESmoothedOpt.mat', 'thetaPhiPatESmoothedOpt');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\thetaPhiPatESmoothedOptNormalized.mat', 'thetaPhiPatESmoothedOptNormalized');

save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\Coeff1.mat', 'Coeff1');
save('C:\Users\AI-BEAM\Documents\MATLAB\Optimization Data\Coeff2.mat', 'Coeff2');
