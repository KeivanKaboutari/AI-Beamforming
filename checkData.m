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
global ComplexPhaseFactor optComplexPhaseFactor modifiedOptComplexPhaseFactor;

% Amplitude and phase of the beams, with the same dimension L (number of beams)
amp = [1, 1];
phase = [0, pi];

%% Load data
% Load M and N
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\M.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\N.mat');

% Load Theta and Phi of the source
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\thetaSource.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\phiSource.mat');

global sx sy;
sx = M / 1;
sy = N / 1;

% Confine Function is confining the phase within the [-pi, pi] interval.
Confine = @(x) pi - mod(pi - x, 2 * pi);

% Load L, d, and Lambda
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\L.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\d.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Lambda.mat');

load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPhaseWoN.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPhaseWN.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPhaseOpt.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\modifiedPhaseOpt.mat');

% Load project pattern
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPatCRSmoothed.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPatCRSmoothedNormalized.mat');

% Load electric field (before and after optimization)
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPatESmoothed.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPatESmoothedNormalized.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPatESmoothedOpt.mat');
load('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\thetaPhiPatESmoothedOptNormalized.mat');

% Number of data sets
DataSet = 1;
% Number of points per data set (Number of angle samples)
PointsNum = 1;
% Number of antenna elements
AntennaElements = N * M;

% Smoothing Factor
SmoothingFactor = 10;

% Convert the angles to distance regarding period of the surface
Lx = Rf * tan(pi * phiSource / 180) / d;
Ly = Rf * tan(pi * thetaSource / 180) / d;

%% Select which Theta and Phi to plot
selectThetaPhi = 1;
% Select data set
selectDataSet = 1;

ComplexPhaseFactor = 1 * exp(1i * reshape(thetaPhiPhaseWN(selectThetaPhi, 2 * L + 1 : end, selectDataSet), N, M));
optComplexPhaseFactor = 1 * exp(1i * reshape(thetaPhiPhaseOpt(selectThetaPhi, 2 * L + 1 : end, selectDataSet), N, M));
modifiedOptComplexPhaseFactor = 1 * exp(1i * reshape(modifiedPhaseOpt(selectThetaPhi, :, selectDataSet), N, M));

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


% Generate grid to plot in real dimension on the MS
SizeEleX = linspace(-(M - 1) * d / 2, (M - 1) * d / 2, M);
SizeEleY = linspace(-(N - 1) * d / 2, (N - 1) * d / 2, N);

% Complex value of radiation pattern's samples at different values of u and v
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

% k and l are number of elements in each Row and Column, respectively.
kValues = linspace(0, M - 1, M);
lValues = linspace(0, N - 1, N);
% Convert coordinate to UV space
% Column
UVBaseMeshX = 2 * pi * (kValues - (M - 1) / 2) / (M - 1);
% Row
UVBaseMeshY = 2 * pi * (lValues - (N - 1) / 2) / (N - 1);
[UVBaseX, UVBaseY] = meshgrid(UVBaseMeshX, UVBaseMeshY);

% Convert coordinate to UV space
% Column
UVBaseMeshXk0 = 2 * k0 * d * (kValues - (M - 1) / 2) / (M - 1);
% Row
UVBaseMeshYk0 = 2 * k0 * d * (lValues - (N - 1) / 2) / (N - 1);
% Generate mesh grid
[gridUVBaseX, gridUVBaseY] = meshgrid(UVBaseMeshXk0, UVBaseMeshYk0);

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
        loadedPhaseWoN = reshape(thetaPhiPhaseWoN(PointCounter, 2 * L + 1 : end, DataSetCounter), N, M);
        loadedPhaseWN = reshape(thetaPhiPhaseWN(PointCounter, 2 * L + 1 : end, DataSetCounter), N, M);
        loadedPhaseOpt = reshape(thetaPhiPhaseOpt(PointCounter, 2 * L + 1 : end, DataSetCounter), N, M);
        loadedModifiedPhaseOpt = reshape(modifiedPhaseOpt(PointCounter, :, DataSetCounter), N, M);
        
        thetaPhiPhaseWoN(PointCounter, 2 * L + 1 : end, DataSetCounter) = reshape(transpose(loadedPhaseWoN), 1, AntennaElements);
        thetaPhiPhaseWN(PointCounter, 2 * L + 1 : end, DataSetCounter) = reshape(transpose(loadedPhaseWN), 1, AntennaElements);
        thetaPhiPhaseOpt(PointCounter, 2 * L + 1 : end, DataSetCounter) = reshape(transpose(loadedPhaseOpt), 1, AntennaElements);
        modifiedPhaseOpt(PointCounter, :, DataSetCounter) = reshape(transpose(loadedModifiedPhaseOpt), 1, AntennaElements);
    end
end

PhaseWithoutNoise = [];
for DataSetCounter = 1 : 1 : DataSet
    PhaseWithoutNoise = [PhaseWithoutNoise; thetaPhiPhaseWoN(:, :, DataSetCounter)];
end
save('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\Modified\PhaseWithoutNoise.mat', 'PhaseWithoutNoise');

PhaseWithNoise = [];
for DataSetCounter = 1 : 1 : DataSet
    PhaseWithNoise = [PhaseWithNoise; thetaPhiPhaseWN(:, :, DataSetCounter)];
end
save('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\Modified\PhaseWithNoise.mat', 'PhaseWithNoise');

PhaseOpt = [];
ModifiedPhaseOpt = [];
for DataSetCounter = 1 : 1 : DataSet
    PhaseOpt = [PhaseOpt; thetaPhiPhaseOpt(:, :, DataSetCounter)];
    ModifiedPhaseOpt = [ModifiedPhaseOpt; modifiedPhaseOpt(:, :, DataSetCounter)];
end
save('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\Modified\PhaseOpt.mat', 'PhaseOpt');
save('C:\Users\k.kaboutari\Desktop\13by13 MS by Abdel\Data\twoBeam13by13ElementPhi0Theta-30Phi0Theta30\Final Results\Modified\ModifiedPhaseOpt.mat', 'ModifiedPhaseOpt');

%% Plot results related to the project pattern
% Plot 2D and 3d complex value of radiation pattern in xy and uv space
ComplexRPMaxSmoothed = max(transpose(loadedCRSmoothed), [], 'all');
Plot2Dand3D('Sample Pattern in uv-space', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', loadedCRSmoothed, '|C(u,v)|', ComplexRPMaxSmoothed, [1, 1], [0, 0], 0);
% Smoothed complex radiation pattern in spherical coordinate
SphericalPlot3D('Complex radiation pattern in spherical space', '|C(x,y)|', 'CRP', Lx, Ly, amp, phase);
% Plot of smoothed CRP in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Project pattern on the xy-plane (Theta = 90 deg)', 'xy-plane', 'CRP', Lx, Ly, amp, phase);
% Plot of smoothed CRP in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Project pattern on the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'CRP', Lx, Ly, amp, phase);
% Plot of smoothed CRP in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Project pattern on the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'CRP', Lx, Ly, amp, phase);

%% Plot results related to the electrical field
ConfinedSincPatternPhaseNoisy = atan2(imag(ComplexPhaseFactor), real(ComplexPhaseFactor));

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

AbsElectricFieldSmoothedMax = max(loadedAbsElectricFieldSmoothed, [], 'all');
% Plot electric field in uv coordinate
Plot2Dand3D('2D and 3D Far-Field Radiation Pattern (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', loadedAbsElectricFieldSmoothed, '|E(u,v)|', AbsElectricFieldSmoothedMax, [1, 1], [0, 0], 0);
% Smoothed electrical field in spherical coordinate
SphericalPlot3D('Electrical field in spherical space', 'Normalized |E(x,y)|', 'EP', Lx, Ly, amp, phase);
% Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Electric pattern on the xy-plane (Theta = 90 deg)', 'xy-plane', 'EP', Lx, Ly, amp, phase);
% Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Electric pattern on the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'EP', Lx, Ly, amp, phase);
% Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Electric pattern on the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'EP', Lx, Ly, amp, phase);

%% Plot results related to the optimized electrical field
ConfinedSincPatternPhaseNoisy = atan2(imag(optComplexPhaseFactor), real(optComplexPhaseFactor));

ElemBasePhase = @(y, x) Confine(ConfinedSincPatternPhaseNoisy(1 + round(y / d + (N - 1) / 2), 1 + round(x / d + (M - 1) / 2)));

% Let us caclulate and plot the phase distribution on the MS that realises the given pattern
ConfinedElemBasePhase2D = zeros(N, M);
for Row = 1 : 1 : N
    for Column = 1 : 1 : M
        ConfinedElemBasePhase2D(Row, Column) = ElemBasePhase(SizeEleY(Row), SizeEleX(Column));
    end
end

% Plot phase distribution on MS
PlotPhase('2D Optimized Phase Distributaion', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', ConfinedElemBasePhase2D, 'Element based', ...
                                    gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
                                    SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
                                                                                                                              [1, 1, 1]);

% Plot 2D and 3D Optimized Far-Field Radiation Pattern in uv-Base
Plot2Dand3D('2D and 3D Optimized Far-Field Radiation Pattern (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', loadedAbsElectricFieldSmoothedOpt, '|E(u,v)|', optAbsElectricFieldSmoothedMax, [1, 1], [0, 0], 0);
% Smoothed electrical field in spherical coordinate
SphericalPlot3D('Optimized electrical field in spherical space', 'Normalized |E_{Opt.}(x,y)|', 'OptEP', Lx, Ly, amp, phase);
% Plot of optimized electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Optimized electric pattern on the xy-plane (Theta = 90 deg)', 'xy-plane', 'OptEP', Lx, Ly, amp, phase);
% Plot of optimized electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Optimized electric pattern on the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'OptEP', Lx, Ly, amp, phase);
% Plot of optimized electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Optimized electric pattern on the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'OptEP', Lx, Ly, amp, phase);
% Plot of smoothed CRP, EF, and OptEF in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Project, EF, and Optimized EF patterns on the xy-plane (Theta = 90 deg)', 'xy-plane', 'All', Lx, Ly, amp, phase);
% Plot of smoothed CRP, EF, and OptEF in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Project, EF, and Optimized EF patterns on the yz-plane (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'All', Lx, Ly, amp, phase);
% Plot of smoothed CRP, EF, and OptEF in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Project, EF, and Optimized EF patterns on the xz-plane (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'All', Lx, Ly, amp, phase);

%% Plot results related to the optimized electrical field related to the modified phases
ConfinedSincPatternPhaseNoisy = atan2(imag(modifiedOptComplexPhaseFactor), real(modifiedOptComplexPhaseFactor));

ElemBasePhase = @(y, x) Confine(ConfinedSincPatternPhaseNoisy(1 + round(y / d + (N - 1) / 2), 1 + round(x / d + (M - 1) / 2)));

% Let us caclulate and plot the phase distribution on the MS that realises the given pattern
ConfinedElemBasePhase2D = zeros(N, M);
for Row = 1 : 1 : N
    for Column = 1 : 1 : M
        ConfinedElemBasePhase2D(Row, Column) = ElemBasePhase(SizeEleY(Row), SizeEleX(Column));
    end
end

% Plot phase distribution on MS
PlotPhase('2D Modified and Optimized Phase Distributaion', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', ConfinedElemBasePhase2D, 'Element based', ...
                                    gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
                                    SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
                                                                                                                              [1, 1, 1]);

%% Plot phases without Noise
for figureCounter = 1 : 1 : PointsNum
    Text = append('Phase Without Noise for Theta and Phi no. ', num2str(figureCounter));
    figure('Name', Text);
    hold on;
    grid on;
    for dataCounter = 1 : 1 : DataSet
        plot(1 : 1 : AntennaElements, reshape(transpose(reshape(thetaPhiPhaseWoN(figureCounter, 2 * L + 1 : end, dataCounter), M, N)), 1, AntennaElements));
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
        plot(1 : 1 : AntennaElements, reshape(transpose(reshape(thetaPhiPhaseWN(figureCounter, 2 * L + 1 : end, dataCounter), M, N)), 1, AntennaElements));
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
        plot(1 : 1 : AntennaElements, reshape(transpose(reshape(thetaPhiPhaseOpt(figureCounter, 2 * L + 1 : end, dataCounter), M, N)), 1, AntennaElements));
        ylim([-pi pi]);
    end
end
hold off;

%% Plot modified optimized phases
for figureCounter = 1 : 1 : PointsNum
    Text = append('Optimized and modified phase for Theta and Phi no. ', num2str(figureCounter));
    figure('Name', Text);
    hold on;
    grid on;
    for dataCounter = 1 : 1 : DataSet
        plot(1 : 1 : AntennaElements, reshape(transpose(reshape(modifiedPhaseOpt(figureCounter, :, dataCounter), M, N)), 1, AntennaElements));
        ylim([-pi pi]);
    end
end
hold off;
