function PlotPhase(FigureName, Xaxis1, Xlabel1, Yaxis1, Ylabel1, Zvalue, Title1, ...
                               Xaxis2, Xlabel2, Yaxis2, Ylabel2, Title2, ...
                               SizeEleX, Xlabel3, SizeEleY, Ylabel3, Title3, ...
                                                                             SelectPlot)
% SelectPlot = [?, ?, ?] means plot which phases, (Element-base, uv-space, Physical dimension), 1 means plot and 0 means do not plot
    global d;
    global left_color right_color;
    [N, M] = size(Xaxis1);
    
    % MS dimension in uv-space
    UaxisMax = max(Xaxis2, [], 'all');
    UaxisMin = min(Xaxis2, [], 'all');
    Uaxis = linspace(UaxisMin, UaxisMax, M);
    VaxisMax = max(Yaxis2, [], 'all');
    VaxisMin = min(Yaxis2, [], 'all');
    Vaxis = linspace(VaxisMin, VaxisMax, N);
    
    resizeScale = 500;
    Phase = imresize(Zvalue, resizeScale, "method", "nearest", "Antialiasing", false);
    [Nsize, Msize] = size(Phase);
    % Dimension of the MS
    MaxPhase = max(Phase, [], 'all');
    MinPhase = min(Phase, [], 'all');
    
    % Define boundary lines
    LineX = linspace(0, Msize, N * M);
    LineYUp = zeros(1, N * M);
    LineYDown = ones(1, N * M) * Nsize;
    LineY = linspace(0, Nsize, N * M);
    LineXLeft = zeros(1, N * M);
    LineXRight = ones(1, N * M) * Msize;
    
    % Definge figure window
    Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
    
    if SelectPlot == [1, 0, 0]
        Fig2DElem = imshow(Phase);
        set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
        title(Title1);
        colorbar;
        colormap(cool(256));
        clim([MinPhase MaxPhase]);
        % Setting the axes label, type, position, rotation, etc.
        xlabel(Xlabel1, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
        ylabel(Ylabel1, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
        % Enhancing the surface characteristics
        axis tight;
        axis equal;
        hold on;
        % Plot outer boundary lines
        plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
        plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
        plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
        plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
        % Plot inner boundary lines
        for Counter = 1 : 1 : N - 1
            plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
        end
        for Counter = 1 : 1 : M - 1
            plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
        end
        hold off;
        xlim([-inf inf]);
        xticks([resizeScale / 2 : resizeScale : Msize]);
        for Counter = 1 : 1 : M
            LabelX(1, Counter) = {num2str(Counter)};
        end
        xticklabels(LabelX);
        ylim([-inf inf]);
        yticks([resizeScale / 2 : resizeScale : Nsize]);
        for Counter = 1 : 1 : N
            LabelY(1, Counter) = {num2str(Counter)};
        end
        yticklabels(LabelY);
    elseif SelectPlot == [1, 1, 0]
        subplot(2, 1, 1);
        Fig2DElem = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title1);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel1, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel1, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(Counter)};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(Counter)};
            end
            yticklabels(LabelY);
        subplot(2, 1, 2);
            Fig2DuvBase = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title2);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel2, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel2, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(Uaxis(Counter))};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(Vaxis(Counter))};
            end
            yticklabels(LabelY);
    elseif SelectPlot == [1, 1, 1]
        subplot(3, 1, 1);
            Fig2DElem = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title1);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel1, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel1, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(Counter)};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(Counter)};
            end
            yticklabels(LabelY);
        subplot(3, 1, 2);
            Fig2DuvBase = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title2);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel2, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel2, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(Uaxis(Counter))};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(Vaxis(Counter))};
            end
            yticklabels(LabelY);
        subplot(3, 1, 3);
            Fig2DPhyDim = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title3);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel3, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel3, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(SizeEleX(Counter))};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(SizeEleY(Counter))};
            end
            yticklabels(LabelY);
    elseif SelectPlot == [0, 1, 0]
            Fig2DuvBase = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title2);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel2, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel2, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(Uaxis(Counter))};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(Vaxis(Counter))};
            end
            yticklabels(LabelY);
    elseif SelectPlot == [0, 1, 1]
        subplot(2, 1, 1);
            Fig2DuvBase = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title2);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel2, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel2, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(Uaxis(Counter))};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(Vaxis(Counter))};
            end
            yticklabels(LabelY);
        subplot(2, 1, 2);
            Fig2DPhyDim = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title3);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel3, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel3, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(SizeEleX(Counter))};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(SizeEleY(Counter))};
            end
            yticklabels(LabelY);
    elseif SelectPlot == [0, 0, 1]
        Fig2DPhyDim = imshow(Phase);
        set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
        title(Title3);
        colorbar;
        colormap(cool(256));
        clim([MinPhase MaxPhase]);
        % Setting the axes label, type, position, rotation, etc.
        xlabel(Xlabel3, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
        ylabel(Ylabel3, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
        % Enhancing the surface characteristics
        axis tight;
        axis equal;
        hold on;
        % Plot outer boundary lines
        plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
        plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
        plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
        plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
        % Plot inner boundary lines
        for Counter = 1 : 1 : N - 1
            plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
        end
        for Counter = 1 : 1 : M - 1
            plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
        end
        hold off;
        xlim([-inf inf]);
        xticks([resizeScale / 2 : resizeScale : Msize]);
        for Counter = 1 : 1 : M
            LabelX(1, Counter) = {num2str(SizeEleX(Counter))};
        end
        xticklabels(LabelX);
        ylim([-inf inf]);
        yticks([resizeScale / 2 : resizeScale : Nsize]);
        for Counter = 1 : 1 : N
            LabelY(1, Counter) = {num2str(SizeEleY(Counter))};
        end
        yticklabels(LabelY);
    elseif SelectPlot == [1, 0, 1]
        subplot(2, 1, 1);
            Fig2DElem = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title1);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel1, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel1, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(Counter)};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(Counter)};
            end
            yticklabels(LabelY);
        subplot(2, 1, 2);
            Fig2DPhyDim = imshow(Phase);
            set(gca, 'Visible', 'On', 'XAxisLocation', 'top', 'YAxisLocation', 'left');
            title(Title3);
            colorbar;
            colormap(cool(256));
            clim([MinPhase MaxPhase]);
            % Setting the axes label, type, position, rotation, etc.
            xlabel(Xlabel3, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            ylabel(Ylabel3, 'fontweight', 'bold', 'FontSize' , 12 , 'FontName' , 'Times New Roman');
            % Enhancing the surface characteristics
            axis tight;
            axis equal;
            hold on;
            % Plot outer boundary lines
            plot(LineX, LineYUp, 'Color', 'black', 'LineWidth', 2);
            plot(LineX, LineYDown, 'Color', 'black', 'LineWidth', 2);
            plot(LineXLeft, LineY, 'Color', 'black', 'LineWidth', 2);
            plot(LineXRight, LineY, 'Color', 'black', 'LineWidth', 2);
            % Plot inner boundary lines
            for Counter = 1 : 1 : N - 1
                plot(LineX, LineYUp + Counter * resizeScale, 'Color', 'black', 'LineWidth', 1.5);
            end
            for Counter = 1 : 1 : M - 1
                plot(LineXLeft + Counter * resizeScale + 0, LineY, 'Color', 'black', 'LineWidth', 1.5);
            end
            hold off;
            xlim([-inf inf]);
            xticks([resizeScale / 2 : resizeScale : Msize]);
            for Counter = 1 : 1 : M
                LabelX(1, Counter) = {num2str(SizeEleX(Counter))};
            end
            xticklabels(LabelX);
            ylim([-inf inf]);
            yticks([resizeScale / 2 : resizeScale : Nsize]);
            for Counter = 1 : 1 : N
                LabelY(1, Counter) = {num2str(SizeEleY(Counter))};
            end
            yticklabels(LabelY);
    end
    
end