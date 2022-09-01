clear all;
close all;
clc;

%% Load selected data
BeamNumber = 1;

%% Set axes color for figures
global left_color right_color;
left_color = [0 0 0];
right_color = [0 0 0];

global Rf;
Rf = 50.000000;
global M N d k0 SmoothingFactor;
global ComplexPhaseFactor;

%% Load data
% Load M and N
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\M.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\N.mat');

global sx sy;
sx = M / 1;
sy = N / 1;

% Confine Function is confining the phase within the [-pi, pi] interval.
Confine = @(x) pi - mod(pi - x, 2 * pi);

% Load L, d, and Lambda
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\L.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\d.mat');
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\Lambda.mat');

% Calculate k0
k0 = 2 * pi / Lambda;

% Number of antenna elements
AntennaElements = N * M;

% Smoothing Factor
SmoothingFactor = 5;

% Smoothed mesh regarding Massive MIMO Antenna
kValuesSmoothed = linspace(0, M - 1, round(SmoothingFactor * M));
lValuesSmoothed = linspace(0, N - 1, round(SmoothingFactor * N));

% Generate mesh grid for plotting E-field in uv-space and finding exact loacation of beams
% For plotting E-field in uv-space and finding exact loacation of beams
MaxLocUVBaseMeshXSmoothed = 2 * pi * (kValuesSmoothed - (M - 1) / 2) / (M - 1);
MaxLocUVBaseMeshYSmoothed = 2 * pi * (lValuesSmoothed - (N - 1) / 2) / (N - 1);
[gridMaxLocUVBaseXSmoothed, gridMaxLocUVBaseYSmoothed] = meshgrid(MaxLocUVBaseMeshXSmoothed, MaxLocUVBaseMeshYSmoothed);

% Generate grid to plot in real dimension on the MS
SizeEleX = linspace(-(M - 1) * d / 2, (M - 1) * d / 2, M);
SizeEleY = linspace(-(N - 1) * d / 2, (N - 1) * d / 2, N);

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

%%
% Load different phases to check the directions (original data generated by algorithm)
% load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\thetaPhiPhaseOpt.mat');
% ComplexPhaseFactor = 1 * exp(1i * reshape(thetaPhiPhaseOpt(BeamNumber, 2 * L + 1 : end, 1), N, M));
% 
% % Calculate and plot
% ConfinedSincPatternPhaseNoisy = atan(real(ComplexPhaseFactor) ./ imag(ComplexPhaseFactor));
% 
% ElemBasePhase = @(y, x) Confine(ConfinedSincPatternPhaseNoisy(1 + round(y / d + (N - 1) / 2), 1 + round(x / d + (M - 1) / 2)));
% 
% % Let us caclulate and plot the phase distribution on the MS that realises the given pattern
% ConfinedElemBasePhase2D = zeros(N, M);
% for Row = 1 : 1 : N
%     for Column = 1 : 1 : M
%         ConfinedElemBasePhase2D(Row, Column) = ElemBasePhase(SizeEleY(Row), SizeEleX(Column));
%     end
% end
% 
% % Plot phase distribution on MS
% PlotPhase('2D Phase Distributaion by original data generated by algorithm', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', ConfinedElemBasePhase2D, 'Element based', ...
%                                     gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
%                                     SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
%                                                                                                                               [1, 1, 1]);
% 
% % Calculate Smoothed Electrical Field
% AbsElectricFieldSmoothedOriginalData = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
% for kCounter = 1 : 1 : round(SmoothingFactor * M)
%     u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
%     % u = 2 * k0 * d * (kValuesSmoothed(kCounter) - (M - 1) / 2) / (M - 1);
%     for lCounter = 1 : 1 : round(SmoothingFactor * N)
%         v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
%         % v = 2 * k0 * d * (lValuesSmoothed(lCounter) - (N - 1) / 2) / (N - 1);
%         AbsElectricFieldSmoothedOriginalData(lCounter, kCounter) = ElecFieldFun(ComplexPhaseFactor, [u, v]);
%     end
% end
% 
% % Plot electric field in uv coordinate
% Plot2Dand3D('2D and 3D Far-Field Radiation Pattern by original data generated by algorithm (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', AbsElectricFieldSmoothedOriginalData, '|E(u,v)|', max(AbsElectricFieldSmoothedOriginalData, [], 'all'), [1, 1], [0, 0], 0);
% 
% % Smoothed electrical field in spherical coordinate
% SphericalPlot3D('Electrical field in spherical space by original data generated by algorithm', 'Normalized |E(x,y)|', 'EP');
% 
% % Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
% PolarPlot2D('Electric pattern in the xy-plane by original data generated by algorithm (Theta = 90 deg)', 'xy-plane', 'EP');
% 
% % Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
% PolarPlot2D('Electric pattern in the yz-plane by original data generated by algorithm (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'EP');
% 
% % Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
% PolarPlot2D('Electric pattern in the xz-plane by original data generated by algorithm (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'EP');

%%
% Load different phases to check the directions (shaped data by checkData script)
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\PhaseData\PhaseOpt.mat');
ComplexPhaseFactor = 1 * exp(1i * reshape(transpose(reshape(PhaseOpt(BeamNumber, 2 * L + 1 : end), M, N)), N, M));

% Calculate and plot
ConfinedSincPatternPhaseNoisy = atan(real(ComplexPhaseFactor) ./ imag(ComplexPhaseFactor));

ElemBasePhase = @(y, x) Confine(ConfinedSincPatternPhaseNoisy(1 + round(y / d + (N - 1) / 2), 1 + round(x / d + (M - 1) / 2)));

% Let us caclulate and plot the phase distribution on the MS that realises the given pattern
ConfinedElemBasePhase2D = zeros(N, M);
for Row = 1 : 1 : N
    for Column = 1 : 1 : M
        ConfinedElemBasePhase2D(Row, Column) = ElemBasePhase(SizeEleY(Row), SizeEleX(Column));
    end
end

% Plot phase distribution on MS
PlotPhase('2D Phase Distributaion using shaped data by checkData script', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', ConfinedElemBasePhase2D, 'Element based', ...
                                    gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
                                    SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
                                                                                                                              [1, 1, 1]);

% Calculate Smoothed Electrical Field
AbsElectricFieldSmoothedShapedData = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
for kCounter = 1 : 1 : round(SmoothingFactor * M)
    u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
    % u = 2 * k0 * d * (kValuesSmoothed(kCounter) - (M - 1) / 2) / (M - 1);
    for lCounter = 1 : 1 : round(SmoothingFactor * N)
        v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
        % v = 2 * k0 * d * (lValuesSmoothed(lCounter) - (N - 1) / 2) / (N - 1);
        AbsElectricFieldSmoothedShapedData(lCounter, kCounter) = ElecFieldFun(ComplexPhaseFactor, [u, v]);
    end
end

% Plot electric field in uv coordinate
Plot2Dand3D('2D and 3D Far-Field Radiation Pattern using shaped data by checkData script (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', AbsElectricFieldSmoothedShapedData, '|E(u,v)|', max(AbsElectricFieldSmoothedShapedData, [], 'all'), [1, 1], [0, 0], 0);

% Smoothed electrical field in spherical coordinate
SphericalPlot3D('Electrical field in spherical space using shaped data by checkData script', 'Normalized |E(x,y)|', 'EP');

% Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Electric pattern in the xy-plane using shaped data by checkData script (Theta = 90 deg)', 'xy-plane', 'EP');

% Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Electric pattern in the yz-plane using shaped data by checkData script (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'EP');

% Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Electric pattern in the xz-plane using shaped data by checkData script (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'EP');

%%
% Load different modified phases to check the directions (shaped data by checkData script)
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\PhaseData\ModifiedPhaseOpt.mat');
ComplexPhaseFactor = 1 * exp(1i * reshape(transpose(reshape(ModifiedPhaseOpt(BeamNumber, :), M, N)), N, M));

% Calculate and plot
ConfinedSincPatternPhaseNoisy = atan(real(ComplexPhaseFactor) ./ imag(ComplexPhaseFactor));

ElemBasePhase = @(y, x) Confine(ConfinedSincPatternPhaseNoisy(1 + round(y / d + (N - 1) / 2), 1 + round(x / d + (M - 1) / 2)));

% Let us caclulate and plot the phase distribution on the MS that realises the given pattern
ConfinedElemBasePhase2D = zeros(N, M);
for Row = 1 : 1 : N
    for Column = 1 : 1 : M
        ConfinedElemBasePhase2D(Row, Column) = ElemBasePhase(SizeEleY(Row), SizeEleX(Column));
    end
end

% Plot phase distribution on MS
PlotPhase('2D Modified Phase Distributaion using shaped data by checkData script', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', ConfinedElemBasePhase2D, 'Element based', ...
                                    gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
                                    SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
                                                                                                                              [1, 1, 1]);

% Calculate Smoothed Electrical Field
AbsElectricFieldSmoothedShapedData = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
for kCounter = 1 : 1 : round(SmoothingFactor * M)
    u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
    % u = 2 * k0 * d * (kValuesSmoothed(kCounter) - (M - 1) / 2) / (M - 1);
    for lCounter = 1 : 1 : round(SmoothingFactor * N)
        v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
        % v = 2 * k0 * d * (lValuesSmoothed(lCounter) - (N - 1) / 2) / (N - 1);
        AbsElectricFieldSmoothedShapedData(lCounter, kCounter) = ElecFieldFun(ComplexPhaseFactor, [u, v]);
    end
end

% Plot electric field in uv coordinate
Plot2Dand3D('2D and 3D Far-Field Radiation Pattern using modified shaped data by checkData script (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', AbsElectricFieldSmoothedShapedData, '|E(u,v)|', max(AbsElectricFieldSmoothedShapedData, [], 'all'), [1, 1], [0, 0], 0);

% Smoothed electrical field in spherical coordinate
SphericalPlot3D('Electrical field in spherical space using modified shaped data by checkData script', 'Normalized |E(x,y)|', 'EP');

% Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Electric pattern in the xy-plane using modified shaped data by checkData script (Theta = 90 deg)', 'xy-plane', 'EP');

% Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Electric pattern in the yz-plane using modified shaped data by checkData script (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'EP');

% Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Electric pattern in the xz-plane using modified shaped data by checkData script (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'EP');

%%
% Load neural networks data
load('C:\Users\k.kaboutari\Desktop\Intelligent Beamforming Metasurfaces for Future Telecommunications (MATLAB Codes)\Data\One beam in (-300,0) for experimental studies\PhaseData\NN_42Dir_1BEAM_val_WON.mat');
ComplexPhaseFactor = 1 * (reshape(transpose(reshape(NN_data_val(BeamNumber, :), M, N)), N, M));

% Calculate and plot
ConfinedSincPatternPhaseNoisy = atan(real(ComplexPhaseFactor) ./ imag(ComplexPhaseFactor));

ElemBasePhase = @(y, x) Confine(ConfinedSincPatternPhaseNoisy(1 + round(y / d + (N - 1) / 2), 1 + round(x / d + (M - 1) / 2)));

% Let us caclulate and plot the phase distribution on the MS that realises the given pattern
ConfinedElemBasePhase2D = zeros(N, M);
for Row = 1 : 1 : N
    for Column = 1 : 1 : M
        ConfinedElemBasePhase2D(Row, Column) = ElemBasePhase(SizeEleY(Row), SizeEleX(Column));
    end
end

% Plot phase distribution on MS
PlotPhase('2D Phase Distributaion by Neural Networks', gridElemBaseX, 'Columns (m)', gridElemBaseY, 'Row (n)', ConfinedElemBasePhase2D, 'Element based', ...
                                    gridUVBaseX, 'u-axis', gridUVBaseY, 'v-axis', 'UV-space', ...
                                    SizeEleX, 'Column Dimension (x)', SizeEleY, 'Row Dimension (y)', 'Physical dimension', ...
                                                                                                                              [1, 1, 1]);

% Calculate Smoothed Electrical Field
AbsElectricFieldSmoothedNeuralNetworks = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
for kCounter = 1 : 1 : round(SmoothingFactor * M)
    u = -pi + kValuesSmoothed(kCounter) * 2 * pi / (M - 1);
    % u = 2 * k0 * d * (kValuesSmoothed(kCounter) - (M - 1) / 2) / (M - 1);
    for lCounter = 1 : 1 : round(SmoothingFactor * N)
        v = -pi + lValuesSmoothed(lCounter) * 2 * pi / (N - 1);
        % v = 2 * k0 * d * (lValuesSmoothed(lCounter) - (N - 1) / 2) / (N - 1);
        AbsElectricFieldSmoothedNeuralNetworks(lCounter, kCounter) = ElecFieldFun(ComplexPhaseFactor, [u, v]);
    end
end

% Plot electric field in uv coordinate
Plot2Dand3D('2D and 3D Far-Field Radiation Pattern by Neural Networks (UV Base)', gridMaxLocUVBaseXSmoothed, 'u-axis', gridMaxLocUVBaseYSmoothed, 'v-axis', AbsElectricFieldSmoothedNeuralNetworks, '|E(u,v)|', max(AbsElectricFieldSmoothedNeuralNetworks, [], 'all'), [1, 1], [0, 0], 0);

% Smoothed electrical field in spherical coordinate
SphericalPlot3D('Electrical field in spherical space by Neural Networks', 'Normalized |E(x,y)|', 'EP');

% Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
PolarPlot2D('Electric pattern in the xy-plane by Neural Networks (Theta = 90 deg)', 'xy-plane', 'EP');

% Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
PolarPlot2D('Electric pattern in the yz-plane by Neural Networks (Phi = 90 or 270 deg, therefore u = 0)', 'yz-plane', 'EP');

% Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
PolarPlot2D('Electric pattern in the xz-plane by Neural Networks (Phi = 0 or 180 deg, therefore v = 0)', 'xz-plane', 'EP');
