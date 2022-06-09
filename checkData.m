clear all;
close all;
clc;

%% Set axes color for figures
global left_color right_color;
left_color = [0 0 0];
right_color = [0 0 0];

global Rf;
Rf = 50.000000;
global M N d k0 SmoothingFactor;
global ComplexPhaseFactor optComplexPhaseFactor;

%% Load data
% Load M and N
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\M.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\N.mat');

global sx sy;
sx = M / 1;
sy = N / 1;

% Confine Function is confining the phase within the [-pi, pi] interval.
Confine = @(x) pi - mod(pi - x, 2 * pi);

% Load L, d, and Lambda
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\L.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\d.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\Lambda.mat');

load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPhaseWoN.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPhaseWN.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPhaseOpt.mat');

% Load project pattern
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPatCRSmoothed.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPatCRSmoothedNormalized.mat');

% Load electric field (before and after optimization)
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPatESmoothed.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPatESmoothedNormalized.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPatESmoothedOpt.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\Two Beam ([50, 0], [50, 180])\thetaPhiPatESmoothedOptNormalized.mat');

% Number of data sets
DataSet = 100;
% Number of points per data set
PointsNum = 1;
% Number of antenna elements
AntennaElements = N * M;

% Smoothing Factor
SmoothingFactor = 5;

%% Select which Theta and Phi to plot
selectThetaPhi = 1;
% Select data set
selectDataSet = 1;

ComplexPhaseFactor = 1 * exp(1i * reshape(thetaPhiPhaseWN(selectThetaPhi, 2 * L + 1 : end, selectDataSet), N, M));
optComplexPhaseFactor = 1 * exp(1i * reshape(thetaPhiPhaseOpt(selectThetaPhi, 2 * L + 1 : end, selectDataSet), N, M));

% ComplexPhaseFactor = transpose(1 * exp(1i * reshape(thetaPhiPhaseWN(selectThetaPhi, 2 * L + 1 : end, selectDataSet), M, N)));
% optComplexPhaseFactor = transpose(1 * exp(1i * reshape(thetaPhiPhaseOpt(selectThetaPhi, 2 * L + 1 : end, selectDataSet), M, N)));

display('[Theta Phi] in rad. =');
for Counter = 1 : 1 : L
    Theta(Counter) = thetaPhiPhaseWoN(selectThetaPhi, (2 * Counter - 1), selectDataSet)
    Phi(Counter) = thetaPhiPhaseWoN(selectThetaPhi, 2 * Counter, selectDataSet)
end

display('[Theta; Phi] in degree =');
[(Theta * 180 / pi); (Phi * 180 / pi)]

k0 = 2 * pi / Lambda;
syms U0 V0;
global u0 v0;
for Counter = 1 : 1 : L
    Equations = [k0 * d * sin(Theta(Counter)) .* cos(Phi(Counter)) - U0 == 0, k0 * d * sin(Theta(Counter)) .* sin(Phi(Counter)) - V0 == 0];
    [u0(Counter), v0(Counter)] = solve(Equations, [U0 V0]);
end

display('[u0, v0] =');
double([u0, v0])

loadedCRSmoothed = reshape(thetaPhiPatCRSmoothed(selectThetaPhi, :, selectDataSet), N * SmoothingFactor, M * SmoothingFactor);
loadedCRSmoothedNormalized = reshape(thetaPhiPatCRSmoothedNormalized(selectThetaPhi, :, selectDataSet),  N * SmoothingFactor, M * SmoothingFactor);
loadedAbsElectricFieldSmoothed = reshape(thetaPhiPatESmoothed(selectThetaPhi, :, selectDataSet),  N * SmoothingFactor, M * SmoothingFactor);
loadedAbsElectricFieldSmoothedNormalized = reshape(thetaPhiPatESmoothedNormalized(selectThetaPhi, :, selectDataSet),  N * SmoothingFactor, M * SmoothingFactor);
loadedAbsElectricFieldSmoothedOpt = reshape(thetaPhiPatESmoothedOpt(selectThetaPhi, :, selectDataSet),  N * SmoothingFactor, M * SmoothingFactor);
loadedAbsElectricFieldSmoothedOptNormalized = reshape(thetaPhiPatESmoothedOptNormalized(selectThetaPhi, :, selectDataSet),  N * SmoothingFactor, M * SmoothingFactor);

% Complex value of radiation pattern's samples at different values of u and v
% Mesh regarding Massive MIMO Antenna rows and columns indices (Element Base)
% Column
ElemBaseMeshX = linspace(-floor(M / 2), floor(M / 2), M);
% Row
ElemBaseMeshY = linspace(-floor(N / 2), floor(N / 2), N);
[ElemBaseX, ElemBaseY] = meshgrid(ElemBaseMeshX, ElemBaseMeshY);
% k and l are number of elements in each Row and Column, respectively.
kValues = linspace(0, M - 1, M);
lValues = linspace(0, N - 1, N);
% Convert coordinate to UV space
% Column
UVBaseMeshX = 2 * pi * (kValues - (M - 1) / 2) / (M - 1);
% Row
UVBaseMeshY = 2 * pi * (lValues - (N - 1) / 2) / (N - 1);
[UVBaseX, UVBaseY] = meshgrid(UVBaseMeshX, UVBaseMeshY);

% Calculated smoothed CRP
kValuesSmoothed = linspace(0, M - 1, round(SmoothingFactor * M));
lValuesSmoothed = linspace(0, N - 1, round(SmoothingFactor * N));

% Column
CRPUVBaseMeshXSmoothed = 2 * k0 * d * (kValuesSmoothed - (M - 1) / 2) / (M - 1);
% Row
CRPUVBaseMeshYSmoothed = 2 * k0 * d * (lValuesSmoothed - (N - 1) / 2) / (N - 1);
% Generate the mesh grid
[gridUVBaseXSmoothed, gridUVBaseYSmoothed] = meshgrid(CRPUVBaseMeshXSmoothed, CRPUVBaseMeshYSmoothed);

% For plotting E-field in uv-space and finding exact loacation of beams
MaxLocUVBaseMeshXSmoothed = 2 * pi * (kValuesSmoothed - (M - 1) / 2) / (M - 1);
MaxLocUVBaseMeshYSmoothed = 2 * pi * (lValuesSmoothed - (N - 1) / 2) / (N - 1);
% Generate mesh grid for plotting E-field in uv-space and finding exact loacation of beams
[gridMaxLocUVBaseXSmoothed, gridMaxLocUVBaseYSmoothed] = meshgrid(MaxLocUVBaseMeshXSmoothed, MaxLocUVBaseMeshYSmoothed);

% Calculate distance from center of the visible area (indeed, MS)
centerDistance = sqrt(power(gridMaxLocUVBaseXSmoothed, 2) + power(gridMaxLocUVBaseYSmoothed, 2));

% Define visible area
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

% Max. of the selected opt. EF
optAbsElectricFieldSmoothedMax = max(loadedAbsElectricFieldSmoothedOpt(centerDistance < visibleAreaRadius));

%% Combine Data
% Phase format conversion
for DataSetCounter = 1 : 1 : DataSet
    for PointCounter = 1 : 1 : PointsNum
        loadedPhaseWoN = reshape(thetaPhiPhaseWoN(PointCounter, 2 * L + 1 : end, DataSetCounter), M, N);
        loadedPhaseWN = reshape(thetaPhiPhaseWN(PointCounter, 2 * L + 1 : end, DataSetCounter), M, N);
        loadedPhaseOpt = reshape(thetaPhiPhaseOpt(PointCounter, 2 * L + 1 : end, DataSetCounter), M, N);

        thetaPhiPhaseWoN(PointCounter, 2 * L + 1 : end, DataSetCounter) = reshape(transpose(loadedPhaseWoN), 1, AntennaElements);
        thetaPhiPhaseWN(PointCounter, 2 * L + 1 : end, DataSetCounter) = reshape(transpose(loadedPhaseWN), 1, AntennaElements);
        thetaPhiPhaseOpt(PointCounter, 2 * L + 1 : end, DataSetCounter) = reshape(transpose(loadedPhaseOpt), 1, AntennaElements);
    end
end

PhaseWithoutNoise = [];
for DataSetCounter = 1 : 1 : DataSet
    PhaseWithoutNoise = [PhaseWithoutNoise; thetaPhiPhaseWoN(:, :, DataSetCounter)];
end
save('C:\Users\k.kaboutari\Documents\MATLAB\Optimization Data\PhaseWithoutNoise.mat', 'PhaseWithoutNoise');

PhaseWithNoise = [];
for DataSetCounter = 1 : 1 : DataSet
    PhaseWithNoise = [PhaseWithNoise; thetaPhiPhaseWN(:, :, DataSetCounter)];
end
save('C:\Users\k.kaboutari\Documents\MATLAB\Optimization Data\PhaseWithNoise.mat', 'PhaseWithNoise');

PhaseOpt = [];
for DataSetCounter = 1 : 1 : DataSet
    PhaseOpt = [PhaseOpt; thetaPhiPhaseOpt(:, :, DataSetCounter)];
end
save('C:\Users\k.kaboutari\Documents\MATLAB\Optimization Data\PhaseOpt.mat', 'PhaseOpt');

%% Plot results related to the project pattern
% Plot 2D and 3d complex value of radiation pattern in xy and uv space
ComplexRPMaxSmoothed = max(transpose(loadedCRSmoothed), [], 'all');
Plot2Dand3D('Sample Pattern in uv-space', gridUVBaseXSmoothed, 'u-axis', gridUVBaseYSmoothed, 'v-axis', loadedCRSmoothed, '|C(u,v)|', ComplexRPMaxSmoothed, [1, 1], [0, 0], 0);
% Smoothed complex radiation pattern in spherical coordinate
SphericalPlot3D('Complex radiation pattern in spherical space', '|C(x,y)|', 'CRP');
% Plot of smoothed CRP in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Project pattern on the xy-plane (Theta = 90 deg)', 'xy-plane', 'CRP');
% Plot of smoothed CRP in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Project pattern on the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'CRP');
% Plot of smoothed CRP in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Project pattern on the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'CRP');

%% Plot results related to the electrical field
AbsElectricFieldSmoothedMax = max(loadedAbsElectricFieldSmoothed, [], 'all');
% Plot electric field in uv coordinate
Plot2Dand3D('2D and 3D Far-Field Radiation Pattern (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', loadedAbsElectricFieldSmoothed, '|E(u,v)|', AbsElectricFieldSmoothedMax, [1, 1], [0, 0], 0);
% Smoothed electrical field in spherical coordinate
SphericalPlot3D('Electrical field in spherical space', 'Normalized |E(x,y)|', 'EP');
% Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Electric pattern on the xy-plane (Theta = 90 deg)', 'xy-plane', 'EP');
% Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Electric pattern on the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'EP');
% Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Electric pattern on the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'EP');

%% Plot results related to the optimized electrical field
% Plot 2D and 3D Optimized Far-Field Radiation Pattern in uv-Base
Plot2Dand3D('2D and 3D Optimized Far-Field Radiation Pattern (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', loadedAbsElectricFieldSmoothedOpt, '|E(u,v)|', optAbsElectricFieldSmoothedMax, [1, 1], [0, 0], 0);
% Smoothed electrical field in spherical coordinate
SphericalPlot3D('Optimized electrical field in spherical space', 'Normalized |E_{Opt.}(x,y)|', 'OptEP');
% Plot of optimized electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Optimized electric pattern on the xy-plane (Theta = 90 deg)', 'xy-plane', 'OptEP');
% Plot of optimized electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Optimized electric pattern on the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'OptEP');
% Plot of optimized electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Optimized electric pattern on the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'OptEP');
% Plot of smoothed CRP, EF, and OptEF in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Project, EF, and Optimized EF patterns on the xy-plane (Theta = 90 deg)', 'xy-plane', 'All');
% Plot of smoothed CRP, EF, and OptEF in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Project, EF, and Optimized EF patterns on the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'All');
% Plot of smoothed CRP, EF, and OptEF in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Project, EF, and Optimized EF patterns on the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'All');

%% Plot phases without noise
for figureCounter = 1 : 1 : PointsNum
    Text = append('Phase without noise for Theta and Phi no. ', num2str(figureCounter));
    figure('Name', Text);
    hold on;
    grid on;
    for dataCounter = 1 : 1 : DataSet
        plot(1 : 1 : AntennaElements, reshape(transpose(reshape(thetaPhiPhaseWoN(figureCounter, 2 * L + 1 : end, dataCounter), N, M)), 1, N * M));
        ylim([-pi pi]);
    end
end
hold off;

%% Plot phases with noise
for figureCounter = 1 : 1 : PointsNum
    Text = append('Phase with noise for Theta and Phi no. ', num2str(figureCounter));
    figure('Name', Text);
    hold on;
    grid on;
    for dataCounter = 1 : 1 : DataSet
        plot(1 : 1 : AntennaElements, reshape(transpose(reshape(thetaPhiPhaseWN(figureCounter, 2 * L + 1 : end, dataCounter), N, M)), 1, N * M));
        ylim([-pi pi]);
    end
end
hold off;

%% Plot optimized phases
for figureCounter = 1 : 1 : PointsNum
    Text = append('Optimized phase for Theta and Phi no. ', num2str(figureCounter));
    figure('Name', Text);
    hold on;
    grid on;
    for dataCounter = 1 : 1 : DataSet
        plot(1 : 1 : AntennaElements, reshape(transpose(reshape(thetaPhiPhaseOpt(figureCounter, 2 * L + 1 : end, dataCounter), N, M)), 1, N * M));
        ylim([-pi pi]);
    end
end
hold off;
