function SphericalPlot3D(FigureName, Zlabel, PatternType, Lx, Ly)
    global M N d k0 SmoothingFactor left_color right_color;

    % Element number
    EleNum = N * M;
    
    % Generate mesh grid for plotting in spherical coordinate using Phi and Theta values
    meshThetaSmoothed = linspace(0 , pi / 2, SmoothingFactor * N);
    meshPhiSmoothed = linspace(0 , 2 * pi, SmoothingFactor * M);
    
    % Grids for plotting in polar and spherical coordinates (2D and 3D)
    Xmn = zeros(SmoothingFactor * N, SmoothingFactor * M);
    Ymn = zeros(SmoothingFactor * N, SmoothingFactor * M);
    for thetaCounter = 1 : 1 : round(SmoothingFactor * N)
        for phiCounter = 1 : 1 : round(SmoothingFactor * M)
            % Generating mesh in spherical coordinate
            % Number of rows
            Xmn(thetaCounter, phiCounter) = sin(meshThetaSmoothed(thetaCounter)) * cos(meshPhiSmoothed(phiCounter));
            % Number of columns
            Ymn(thetaCounter, phiCounter) = sin(meshThetaSmoothed(thetaCounter)) * sin(meshPhiSmoothed(phiCounter));
        end
    end
    Zmn = real(sqrt(1 - power(Xmn, 2) - power(Ymn, 2)));
    
    switch PatternType
        case 'CRP'
            % Spherical plot of the complex radiation pattern
            % Smoothed complex radiation pattern in spherical coordinate
            SphCorAbsCRPSmoothed = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
            for kCounter = 1 : 1 : round(SmoothingFactor * M)
                for lCounter = 1 : 1 : round(SmoothingFactor * N)
                    % Calculate u and v values
                    uSphCor = k0 * d * Xmn(lCounter, kCounter);
                    vSphCor = k0 * d * Ymn(lCounter, kCounter);
                    
                    % Calculate complex radiation pattern
                    SphCorAbsCRPSmoothed(lCounter, kCounter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
                end
            end
            % Scale Xmn, Ymn, and Zmn
            Xaxis = Xmn .* SphCorAbsCRPSmoothed;
            Yaxis = Ymn .* SphCorAbsCRPSmoothed;
            Zvalue = Zmn .* SphCorAbsCRPSmoothed;
        case 'EP'
            global ComplexPhaseFactor;

            % Spherical plot of the Electric Field
            % Smoothed electrical field in spherical coordinate
            SphCorAbsElectricFieldSmoothed = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
            for kCounter = 1 : 1 : round(SmoothingFactor * M)
                for lCounter = 1 : 1 : round(SmoothingFactor * N)
                    % Calculate u and v values
                    uSphCor = k0 * d * Xmn(lCounter, kCounter);
                    vSphCor = k0 * d * Ymn(lCounter, kCounter);
                    % Calculate electrical field
                    SphCorAbsElectricFieldSmoothed(lCounter, kCounter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                end
            end
            % Scale Xmn, Ymn, and Zmn
            Xaxis = Xmn .* SphCorAbsElectricFieldSmoothed;
            Yaxis = Ymn .* SphCorAbsElectricFieldSmoothed;
            Zvalue = Zmn .* SphCorAbsElectricFieldSmoothed;
        case 'OptEP'
            global optComplexPhaseFactor;
            
            % Calculate Smoothed and Optimized Electrical Field
            optAbsElectricFieldSmoothed = zeros(round(SmoothingFactor * N), round(SmoothingFactor * M));
            for kCounter = 1 : 1 : round(SmoothingFactor * M)
                for lCounter = 1 : 1 : round(SmoothingFactor * N)
                    % Calculate u and v values
                    uSphCor = k0 * d * Xmn(lCounter, kCounter);
                    vSphCor = k0 * d * Ymn(lCounter, kCounter);
                    % Calculate electrical field
                    optAbsElectricFieldSmoothed(lCounter, kCounter) = ElecFieldFun(optComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                end
            end

            % Scale Xmn, Ymn, and Zmn
            Xaxis = Xmn .* optAbsElectricFieldSmoothed;
            Yaxis = Ymn .* optAbsElectricFieldSmoothed;
            Zvalue = Zmn .* optAbsElectricFieldSmoothed;
    end
    
    % Change color of the surf regarding direction of pattern
    Radius = sqrt(power(Xaxis, 2) + power(Yaxis, 2) + power(Zvalue, 2));
    
    % Normalize figure factor
    NormalizingFactor = max(sqrt(power(Xaxis, 2) + power(Yaxis, 2) + power(Zvalue, 2)), [], 'all');

    % Plot 3D in spherical coordinate
    Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
    % Setting the axes label, type, position, rotation, etc.
    Text = append('', Zlabel);
%     Surf3D = surf(Xaxis, Yaxis, Zvalue, Radius, 'DisplayName', Text);
    Surf3D = surf(Xaxis / NormalizingFactor, Yaxis / NormalizingFactor, Zvalue / NormalizingFactor, Radius / NormalizingFactor, 'DisplayName', Text);
    colorbar;
    colormap(cool(256));
    hold on;
    % Plot X, Y, and Z axes, respectively
    Arrow3D([0 0 0], [max(Xaxis, [], 'all') + 0.1 * max(Xaxis, [], 'all') 0 0] / NormalizingFactor, '3q', 1.2);
    text((max(Xaxis, [], 'all') + 0.2 * max(Xaxis, [], 'all')) / NormalizingFactor, 0, 0, 'x', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman');
    Arrow3D([0 0 0], [0 max(Yaxis, [], 'all') + 0.1 * max(Yaxis, [], 'all') 0] / NormalizingFactor, '3q', 1.2);
    text(0, (max(Yaxis, [], 'all') + 0.2 * max(Yaxis, [], 'all')) / NormalizingFactor, 0, 'y', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman');
    Arrow3D([0 0 0], [0 0 max(Zvalue, [], 'all') + 0.1 * max(Zvalue, [], 'all')] / NormalizingFactor, '3q', 1.2);
    text(0, 0, (max(Zvalue, [], 'all') + 0.2 * max(Zvalue, [], 'all')) / NormalizingFactor, 'z', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman');
    
%     AxisX = linspace(0, max(Xaxis, [], 'all'), 360);
%     plot3(AxisX, zeros(1, 360), zeros(1, 360), 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '-');
%     AxisY = linspace(0, max(Yaxis, [], 'all'), 360);
%     plot3(zeros(1, 360), AxisY, zeros(1, 360), 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '-');
%     AxisZ = linspace(0, max(Zvalue, [], 'all'), 360);
%     plot3(zeros(1, 360), zeros(1, 360), AxisZ, 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '-');
    
%     % Set the location of X,Y, and Z axes
%     % X crossover with Y axis
%     hAxis.XRuler.FirstCrossoverValue  = 0;
%     % X crossover with Z axis
%     hAxis.XRuler.SecondCrossoverValue  = 0;
%     % Y crossover with X axis
%     hAxis.YRuler.FirstCrossoverValue  = 0;
%     % Y crossover with Z axis
%     hAxis.YRuler.SecondCrossoverValue  = 0;
%     % Z crossover with X axis
%     hAxis.ZRuler.FirstCrossoverValue  = 0;
%     % Z crossover with Y axis
%     hAxis.ZRuler.SecondCrossoverValue = 0;
    
    box off;
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', [])
    set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
    xlim([-inf inf]);
    ylim([-inf inf]);
    zlim([-inf inf]);
    
    % Enhancing the surface characteristics
    axis tight;
    axis equal;
%     camlight;
    lighting phong;
    shading interp;
%     title(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%     xlabel('x', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%     ylabel('y', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
%     zlabel('z', 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
    set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
    legend(Surf3D, 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'best');
    grid off;
    hold off;
    view(150, 10);
end