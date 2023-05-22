function Plot2Dand3D(FigureName, Xaxis, Xlabel, Yaxis, Ylabel, Zvalue, Zlabel, MaxZ, SelectPlot, VisibleProtectedAreas, SamplingArea, Lx, Ly)
% SelectPlot = [?, ?] means plot 2D and/or 3D surfs , 1 means plot and 0 means do not plot
% VisibleProtectedAreas = [?, ?] means plot Visible area and/or Protected area, 1 means plot and 0 means do not plot
% SamplingArea = 1 means plot sampling points
    global k0 d L;
    global ComplexPhaseFactor;
    global EFXunit EFYunit;
    global visibleAreaX visibleAreaY;
    global visibleAreaTheta;
    global FixedOutsideSamplingPointsValues InsideProtectedAreaScatteredPointsProperties;
    global ProtectedAreaPhi;
    global left_color right_color;
    
    if SamplingArea == 0
        if VisibleProtectedAreas == [1, 1]
            if SelectPlot == [1, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                end
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(MaxZ, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg1 Leg2], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
                subplot(1, 2, 2);
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                for PatternValue = 0.1 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2);
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg1 Leg2 ], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            elseif SelectPlot == [1, 0]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(1, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg1 Leg2], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
            elseif SelectPlot == [0, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(1, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(1, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(1, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(1, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                end
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                % Plot visible area in UV-space
                for PatternValue = 0.1 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2);
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg1 Leg2], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            end
        elseif VisibleProtectedAreas == [1, 0]
            if SelectPlot == [1, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(1, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg2], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
                subplot(1, 2, 2);
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                % Plot visible area in UV-space
                for PatternValue = 0.1 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2);
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg2], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            elseif SelectPlot == [1, 0]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(1, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg2], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
            elseif SelectPlot == [0, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                % Plot visible area in UV-space
                for PatternValue = 0.1 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2);
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg2 ], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            end
        elseif VisibleProtectedAreas == [0, 1]
            if SelectPlot == [1, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg1], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
                subplot(1, 2, 2);
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg1], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            elseif SelectPlot == [1, 0]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg1], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
            elseif SelectPlot == [0, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg1], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            end
        elseif VisibleProtectedAreas == [0, 0]
            if SelectPlot == [1, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D ], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
                subplot(1, 2, 2);
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D ], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            elseif SelectPlot == [1, 0]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D ], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
            elseif SelectPlot == [0, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            end
        end
    elseif SamplingArea == 1
        if VisibleProtectedAreas == [1, 1]
            if SelectPlot == [1, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                end
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(MaxZ, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4], 'DisplayName', 'Outside PA Samples');
                for Counter = 1: 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg1 Leg2 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
                subplot(1, 2, 2);
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                for PatternValue = 0.1 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2);
                end
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly) / MaxZ, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
                for Counter = 1 : 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / MaxZ, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg1 Leg2 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            elseif SelectPlot == [1, 0]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                end
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(1, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4], 'DisplayName', 'Outside PA Samples');
                for Counter = 1: 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg1 Leg2 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
            elseif SelectPlot == [0, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                for PatternValue = 0.1 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2);
                end
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly) / MaxZ, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
                for Counter = 1 : 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / MaxZ, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg1 Leg2 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            end
        elseif VisibleProtectedAreas == [1, 0]
            if SelectPlot == [1, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(1, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4], 'DisplayName', 'Outside PA Samples');
                for Counter = 1: 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg2 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
                subplot(1, 2, 2);
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                % Plot visible area in UV-space
                for PatternValue = 0.1 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2);
                end
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly) / MaxZ, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
                for Counter = 1 : 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / MaxZ, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg2 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            elseif SelectPlot == [1, 0]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(1, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4], 'DisplayName', 'Outside PA Samples');
                for Counter = 1: 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg2 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
            elseif SelectPlot == [0, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot visible area in UV-space
                Leg2 = plot3(visibleAreaX, visibleAreaY, repmat(0, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2, 'DisplayName', 'Visible Area (VA)');
                % Plot visible area in UV-space
                for PatternValue = 0.1 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                    plot3(visibleAreaX, visibleAreaY, repmat(PatternValue, 1, length(visibleAreaTheta)), 'green', 'LineWidth', 2);
                end
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly) / MaxZ, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
                for Counter = 1 : 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / MaxZ, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg2 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            end
        elseif VisibleProtectedAreas == [0, 1]
            if SelectPlot == [1, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                end
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4], 'DisplayName', 'Outside PA Samples');
                for Counter = 1: 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg1 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
                subplot(1, 2, 2);
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly) / MaxZ, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
                for Counter = 1 : 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / MaxZ, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg1 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            elseif SelectPlot == [1, 0]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 2 : 1 : L
                        plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(MaxZ, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2);
                    end
                end
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4], 'DisplayName', 'Outside PA Samples');
                for Counter = 1: 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg1 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
            elseif SelectPlot == [0, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                % Plot protected area in UV-space
                if Zlabel == "|C(u,v)|"
                    Leg1 = plot3(EFXunit(1, :, 1) / pi * k0 * d, EFYunit(1, :, 1) / pi * k0 * d, repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter) / pi * k0 * d, EFYunit(1, :, Counter) / pi * k0 * d, repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                else
                    Leg1 = plot3(EFXunit(1, :, 1), EFYunit(1, :, 1), repmat(0, 1, length(ProtectedAreaPhi)), 'red', 'LineWidth', 2, 'DisplayName', 'Protected Area (PA)');
                    for Counter = 1 : 1 : L
                        for PatternValue = 0 : 0.5 : max(Zvalue(:, :), [], 'all') / MaxZ
                            plot3(EFXunit(1, :, Counter), EFYunit(1, :, Counter), repmat(PatternValue, 1, length(ProtectedAreaPhi)), 'r', 'LineWidth', 2);
                        end
                    end
                end
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly) / MaxZ, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
                for Counter = 1 : 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / MaxZ, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg1 Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            end
        elseif VisibleProtectedAreas == [0, 0]
            if SelectPlot == [1, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                subplot(1, 2, 1);
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4], 'DisplayName', 'Outside PA Samples');
                for Counter = 1: 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
                subplot(1, 2, 2);
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly) / MaxZ, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
                for Counter = 1 : 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / MaxZ, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 9, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            elseif SelectPlot == [1, 0]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Surf2D = surf(Xaxis, Yaxis, Zvalue, 'DisplayName', Zlabel);
                hold on;
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2)], Lx, Ly), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.3 0.2 0.4], 'DisplayName', 'Outside PA Samples');
                for Counter = 1: 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter), '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colorbar;
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                title(Zlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf2D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf2D Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
                view(2);
            elseif SelectPlot == [0, 1]
                % Plot 2D and 3D
                Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
                Text = append('Normalized ', Zlabel);
                Surf3D = surf(Xaxis, Yaxis, Zvalue / MaxZ, 'DisplayName', Text);
                hold on;
                Leg3 = scatter3(FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), ElecFieldFun(ComplexPhaseFactor, [FixedOutsideSamplingPointsValues(:, 1), FixedOutsideSamplingPointsValues(:, 2), Lx, Ly]) / MaxZ, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1 0.1 0.1], 'DisplayName', 'Outside PA Samples');
                for Counter = 1 : 1 : L
                    Leg4 = scatter3(InsideProtectedAreaScatteredPointsProperties(:, 1, Counter), InsideProtectedAreaScatteredPointsProperties(:, 2, Counter), InsideProtectedAreaScatteredPointsProperties(:, 3, Counter) / MaxZ, '*', 'MarkerEdgeColor', [1 0.5 1], 'MarkerFaceColor', [1 0.5 1], 'DisplayName', 'Inside PA Samples');
                end
                xlim([-inf inf]);
                ylim([-inf inf]);
                colormap(cool(256));
                % Enhancing the surface characteristics
                axis tight;
                camlight;
                lighting phong;
                shading interp;
                % Setting the axes label, type, position, rotation, etc.
                xlabel(Xlabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                ylabel(Ylabel, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                zlabel(Text, 'fontweight', 'bold', 'FontSize' , 20 , 'FontName' , 'Times New Roman');
                set(Surf3D, 'edgecolor', [0 0 0.4], 'meshstyle', 'both', 'linewidth', 0.15);
                legend([Surf3D Leg3 Leg4], 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'northeast');
                hold off;
            end
        end
    end
end