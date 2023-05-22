function PolarPlot2D(FigureName, Plane, PatternType, Lx, Ly)
    % Plane determines to plot which plane
    % PatternType determines to plot CRP (project pattern) or EF or OptEF of all of them
    global M N d k0 SmoothingFactor left_color right_color;

    % MS number of elements
    elementNo = N * M;
    
    % Generate circles radius and Theta
    PolarRadius = 0.2 : 0.2 : 1;
    Step = 3600;
    PolarTheta = linspace(0, 2 * pi, Step);
    
    % Generate inner circles of the polar coordinate
    for raduisCounter = 1 : 1 : numel(PolarRadius) - 1
        InnerPolarX(:, raduisCounter) = PolarRadius(raduisCounter) * cos(PolarTheta);
        InnerPolarY(:, raduisCounter) = PolarRadius(raduisCounter) * sin(PolarTheta);
    end
    
    % Generate outer circles of the polar coordinate
    OuterPolarXBoarder = PolarRadius(raduisCounter + 1) * cos(PolarTheta);
    OuterPolarYBoarder = PolarRadius(raduisCounter + 1) * sin(PolarTheta);
    % Generates straight lines    
    % 0 degree
    Theta0 = zeros(1, Step) * pi / 180;
    % 30 degree
    Theta30 = 30 * ones(1, Step) * pi / 180;
    % 60 degree
    Theta60 = 60 * ones(1, Step) * pi / 180;
    % 90 degree
    Theta90 = 90 * ones(1, Step) * pi / 180;
    % 120 degree
    Theta120 = 120 * ones(1, Step) * pi / 180;
    % 150 degree
    Theta150 = 150 * ones(1, Step) * pi / 180;
    
    % Generate inner circles of the polar coordinate
    % 0 degree
    PolarX0 = linspace(-1, 1, Step) .* cos(Theta0);
    PolarY0 = linspace(-1, 1, Step) .* sin(Theta0);
    % 30 degree
    PolarX30 = linspace(-1, 1, Step) .* cos(Theta30);
    PolarY30 = linspace(-1, 1, Step) .* sin(Theta30);
    % 60 degree
    PolarX60 = linspace(-1, 1, Step) .* cos(Theta60);
    PolarY60 = linspace(-1, 1, Step) .* sin(Theta60);
    % 90 degree
    PolarX90 = linspace(-1, 1, Step) .* cos(Theta90);
    PolarY90 = linspace(-1, 1, Step) .* sin(Theta90);
    % 120 degree
    PolarX120 = linspace(-1, 1, Step) .* cos(Theta120);
    PolarY120 = linspace(-1, 1, Step) .* sin(Theta120);
    % 150 degree
    PolarX150 = linspace(-1, 1, Step) .* cos(Theta150);
    PolarY150 = linspace(-1, 1, Step) .* sin(Theta150);
    
    switch Plane
            case {'xy-plane', 'yx-plane'}
                titleText = 'on xy-plane.';
                % xy-plane (Theta = 90 deg)
                ThetaValueXY = 90 * pi / 180;
                
                % Generate vector of Phi for plotting electric pattern in the planes
                SmoothedMeshPhiPlane = linspace(0, 2 * pi, SmoothingFactor * elementNo);
                
                xyXmn = zeros(1, SmoothingFactor * elementNo);
                xyYmn = zeros(1, SmoothingFactor * elementNo);
                for Counter = 1 : 1 : SmoothingFactor * elementNo
                    % Number of rows (XY)
                    xyXmn(Counter) = sin(ThetaValueXY) * cos(SmoothedMeshPhiPlane(Counter));
                    % Number of columns (XY)
                    xyYmn(Counter) = sin(ThetaValueXY) * sin(SmoothedMeshPhiPlane(Counter));
                end
                xyZmn = real(sqrt(1 - power(xyXmn, 2) - power(xyYmn, 2)));
                
                switch PatternType
                    case 'CRP'
                        % Setting the legend
                        Text = "Normalized project pattern";
                        
                        % Plot of smoothed pattern in the xy-plane (Theta = 90 deg) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xyXmn(Counter);
                            vSphCor = k0 * d * xyYmn(Counter);
                            % Calculate electrical field
                            xySphCorAbsCRPSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
                        end
                        % Scale Xmn, Ymn, and Zmn
                        Xaxis = xyXmn .* xySphCorAbsCRPSmoothed;
                        Yaxis = xyYmn .* xySphCorAbsCRPSmoothed;
    %                     xyScZmn = xyZmn .* xySphCorAbsCRPSmoothed;
            
                        % OR Plot of smoothed pattern in the xy-plane (Theta = 90 deg) calculated in uv-space and converted to polar system
    %                     for Counter = 1 : 1 : SmoothingFactor * elementNo
    %                         Xaxis(Counter) = abs(ComplexRadiationPattern(k0 * d * cos(SmoothedMeshPhiPlane(Counter)), k0 * d * sin(SmoothedMeshPhiPlane(Counter)))) * cos(SmoothedMeshPhiPlane(Counter));
    %                         Yaxis(Counter) = abs(ComplexRadiationPattern(k0 * d * cos(SmoothedMeshPhiPlane(Counter)), k0 * d * sin(SmoothedMeshPhiPlane(Counter)))) * sin(SmoothedMeshPhiPlane(Counter));
    %                     end
                    case 'EP'
                        % Setting the legend
                        Text = "Normalized electric field pattern before optimization";
                        
                        global ComplexPhaseFactor;
                        
                        % Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xyXmn(Counter);
                            vSphCor = k0 * d * xyYmn(Counter);
                            % Calculate electrical field
                            xySphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        Xaxis = xyXmn .* xySphCorAbsElectricFieldSmoothed;
                        Yaxis = xyYmn .* xySphCorAbsElectricFieldSmoothed;
    %                     xyScZmn = xyZmn .* xySphCorAbsElectricFieldSmoothed;
                    case 'OptEP'
                        % Setting the legend
                        Text = "Normalized electric field pattern after optimization";
                        
                        global optComplexPhaseFactor;
    
                        % Plot of optimized electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xyXmn(Counter);
                            vSphCor = k0 * d * xyYmn(Counter);
                            % Calculate electrical field
                            xySphCorAbsElectricFieldSmoothedOpt(Counter) = ElecFieldFun(optComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        Xaxis = xyXmn .* xySphCorAbsElectricFieldSmoothedOpt;
                        Yaxis = xyYmn .* xySphCorAbsElectricFieldSmoothedOpt;
    %                     xyScZmn = xyZmn .* xySphCorAbsElectricFieldSmoothedOpt;

                    case 'All'
                        % Plot of smoothed pattern in the xy-plane (Theta = 90 deg) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xyXmn(Counter);
                            vSphCor = k0 * d * xyYmn(Counter);
                            % Calculate electrical field
                            xySphCorAbsCRPSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisCRP = xyXmn .* xySphCorAbsCRPSmoothed;
                        YaxisCRP = xyYmn .* xySphCorAbsCRPSmoothed;
                        
                        global ComplexPhaseFactor;
                        % Plot of electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xyXmn(Counter);
                            vSphCor = k0 * d * xyYmn(Counter);
                            % Calculate electrical field
                            xySphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisEF = xyXmn .* xySphCorAbsElectricFieldSmoothed;
                        YaxisEF = xyYmn .* xySphCorAbsElectricFieldSmoothed;
                        
                        global optComplexPhaseFactor;
                        % Plot of optimized electric pattern in the xy-plane (Theta = 90 deg) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xyXmn(Counter);
                            vSphCor = k0 * d * xyYmn(Counter);
                            % Calculate electrical field
                            xySphCorAbsElectricFieldSmoothedOpt(Counter) = ElecFieldFun(optComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisOptEF = xyXmn .* xySphCorAbsElectricFieldSmoothedOpt;
                        YaxisOptEF = xyYmn .* xySphCorAbsElectricFieldSmoothedOpt;
                end
                
            case {'yz-plane', 'zy-plane'}
                titleText = 'on yz-plane.';
                % yz-plane (Phi = 90 or 270 deg, therefore u = 0)
                phiValueYZ = 90 * pi / 180;
                
                % Generate vector of Theta and Phi for plotting electric pattern in the planes
                SmoothedMeshThetaPlane = linspace(-pi / 2, pi / 2, SmoothingFactor * elementNo);
                
                yzXmn = zeros(1, SmoothingFactor * elementNo);
                yzYmn = zeros(1, SmoothingFactor * elementNo);
                for Counter = 1 : 1 : SmoothingFactor * elementNo
                    % Number of rows (YZ)
                    yzXmn(Counter) = sin(SmoothedMeshThetaPlane(Counter)) * cos(phiValueYZ);
                    % Number of columns (YZ)
                    yzYmn(Counter) = sin(SmoothedMeshThetaPlane(Counter)) * sin(phiValueYZ);
                end
                yzZmn = real(sqrt(1 - power(yzXmn, 2) - power(yzYmn, 2)));
                
                switch PatternType
                    case 'CRP'
                        % Setting the legend
                        Text = "Normalized project pattern";
                        
                        % Plot of smoothed pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * yzXmn(Counter);
                            vSphCor = k0 * d * yzYmn(Counter);
                            % Calculate electrical field
                            yzSphCorAbsCRPSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
                        end
                        % Scale Xmn, Ymn, and Zmn
    %                     yzScXmn = yzXmn .* yzSphCorAbsCRPSmoothed;
                        Xaxis = yzYmn .* yzSphCorAbsCRPSmoothed;
                        Yaxis = yzZmn .* yzSphCorAbsCRPSmoothed;
                        
                        % Plot of smoothed pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in uv-space and converted to polar system
    %                     for Counter = 1 : 1 : SmoothingFactor * elementNo
    %                         Xaxis(Counter) = abs(ComplexRadiationPattern(0, k0 * d * sin(SmoothedMeshThetaPlane(Counter)))) * sin(SmoothedMeshThetaPlane(Counter));
    %                         Yaxis(Counter) = abs(ComplexRadiationPattern(0, k0 * d * sin(SmoothedMeshThetaPlane(Counter)))) * cos(SmoothedMeshThetaPlane(Counter));
    %                     end
                    case 'EP'
                        % Setting the legend
                        Text = "Normalized electric field pattern before optimization";
                        
                        global ComplexPhaseFactor;
                        
                        % Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * yzXmn(Counter);
                            vSphCor = k0 * d * yzYmn(Counter);
                            % Calculate electrical field
                            yzSphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
    %                     yzScXmn = yzXmn .* yzSphCorAbsElectricFieldSmoothed;
                        Xaxis = yzYmn .* yzSphCorAbsElectricFieldSmoothed;
                        Yaxis = yzZmn .* yzSphCorAbsElectricFieldSmoothed;
                    case 'OptEP'
                        % Setting the legend
                        Text = "Normalized electric field pattern after optimization";

                        global optComplexPhaseFactor;
                        
                        % Plot of optimized electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * yzXmn(Counter);
                            vSphCor = k0 * d * yzYmn(Counter);
                            % Calculate electrical field
                            yzSphCorAbsElectricFieldSmoothedOpt(Counter) = ElecFieldFun(optComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
    %                     yzScXmn = yzXmn .* yzSphCorAbsElectricFieldSmoothedOpt;
                        Xaxis = yzYmn .* yzSphCorAbsElectricFieldSmoothedOpt;
                        Yaxis = yzZmn .* yzSphCorAbsElectricFieldSmoothedOpt;
                    case 'All'
                        % Plot of smoothed pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * yzXmn(Counter);
                            vSphCor = k0 * d * yzYmn(Counter);
                            % Calculate electrical field
                            yzSphCorAbsCRPSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisCRP = yzYmn .* yzSphCorAbsCRPSmoothed;
                        YaxisCRP = yzZmn .* yzSphCorAbsCRPSmoothed;
        
                        global ComplexPhaseFactor;
                        % Plot of electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * yzXmn(Counter);
                            vSphCor = k0 * d * yzYmn(Counter);
                            % Calculate electrical field
                            yzSphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisEF = yzYmn .* yzSphCorAbsElectricFieldSmoothed;
                        YaxisEF = yzZmn .* yzSphCorAbsElectricFieldSmoothed;
        
                        global optComplexPhaseFactor;
                        % Plot of optimized electric pattern in the yz-plane (Phi = 90 or 270 deg, therefore u = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * yzXmn(Counter);
                            vSphCor = k0 * d * yzYmn(Counter);
                            % Calculate electrical field
                            yzSphCorAbsElectricFieldSmoothedOpt(Counter) = ElecFieldFun(optComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisOptEF = yzYmn .* yzSphCorAbsElectricFieldSmoothedOpt;
                        YaxisOptEF = yzZmn .* yzSphCorAbsElectricFieldSmoothedOpt;
                end
                
            case {'xz-plane', 'zx-plane'}
                titleText = 'on zx-plane.';
                % xz-plane (Phi = 0 or 180 deg, therefore v = 0)
                phiValueXZ = 0 * pi / 180;
                
                % Generate vector of Theta for plotting electric pattern in the planes
                SmoothedMeshThetaPlane = linspace(-pi / 2, pi / 2, SmoothingFactor * elementNo);
                
                xzXmn = zeros(1, SmoothingFactor * elementNo);
                xzYmn = zeros(1, SmoothingFactor * elementNo);
                for Counter = 1 : 1 : SmoothingFactor * elementNo
                    % Number of rows (XZ)
                    xzXmn(Counter) = sin(SmoothedMeshThetaPlane(Counter)) * cos(phiValueXZ);
                    % Number of columns (XZ)
                    xzYmn(Counter) = sin(SmoothedMeshThetaPlane(Counter)) * sin(phiValueXZ);
                end
                xzZmn = real(sqrt(1 - power(xzXmn, 2) - power(xzYmn, 2)));
                
                switch PatternType
                    case 'CRP'
                        % Setting the legend
                        Text = "Normalized project pattern";
                        
                        % Plot of smoothed pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xzXmn(Counter);
                            vSphCor = k0 * d * xzYmn(Counter);
                            % Calculate electrical field
                            xzSphCorAbsCRPSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
                        end
                        % Scale Xmn, Ymn, and Zmn
                        Xaxis = xzXmn .* xzSphCorAbsCRPSmoothed;
    %                     xzScYmn = xzYmn .* xzSphCorAbsCRPSmoothed;
                        Yaxis = xzZmn .* xzSphCorAbsCRPSmoothed;
                        
            %             % Plot of smoothed pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in uv-space and converted to polar system
    %                     for Counter = 1 : 1 : SmoothingFactor * elementNo
    %                         xzScXmn(Counter) = abs(ComplexRadiationPattern(k0 * d * sin(SmoothedMeshThetaPlane(Counter)), 0)) * sin(SmoothedMeshThetaPlane(Counter));
    %                         xzScZmn(Counter) = abs(ComplexRadiationPattern(k0 * d * sin(SmoothedMeshThetaPlane(Counter)), 0)) * cos(SmoothedMeshThetaPlane(Counter));
    %                     end
                    case 'EP'
                        % Setting the legend
                        Text = "Normalized electric field pattern before optimization";
                        
                        global ComplexPhaseFactor;
    
                        % Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xzXmn(Counter);
                            vSphCor = k0 * d * xzYmn(Counter);
                            % Calculate electrical field
                            xzSphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        Xaxis = xzXmn .* xzSphCorAbsElectricFieldSmoothed;
    %                     xzScYmn = xzYmn .* xzSphCorAbsElectricFieldSmoothed;
                        Yaxis = xzZmn .* xzSphCorAbsElectricFieldSmoothed;
                    case 'OptEP'
                        % Setting the legend
                        Text = "Normalized electric field pattern after optimization";
                        
                        global optComplexPhaseFactor;
                        
                        % Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xzXmn(Counter);
                            vSphCor = k0 * d * xzYmn(Counter);
                            % Calculate electrical field
                            xzSphCorAbsElectricFieldSmoothedOpt(Counter) = ElecFieldFun(optComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        Xaxis = xzXmn .* xzSphCorAbsElectricFieldSmoothedOpt;
    %                     xzScYmn = xzYmn .* xzSphCorAbsElectricFieldSmoothedOpt;
                        Yaxis = xzZmn .* xzSphCorAbsElectricFieldSmoothedOpt;
                    case 'All'
                        % Plot of smoothed pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xzXmn(Counter);
                            vSphCor = k0 * d * xzYmn(Counter);
                            % Calculate electrical field
                            xzSphCorAbsCRPSmoothed(Counter) = abs(ComplexRadiationPattern(uSphCor, vSphCor));
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisCRP = xzXmn .* xzSphCorAbsCRPSmoothed;
                        YaxisCRP = xzZmn .* xzSphCorAbsCRPSmoothed;
                        
                        global ComplexPhaseFactor;
                        % Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xzXmn(Counter);
                            vSphCor = k0 * d * xzYmn(Counter);
                            % Calculate electrical field
                            xzSphCorAbsElectricFieldSmoothed(Counter) = ElecFieldFun(ComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisEF = xzXmn .* xzSphCorAbsElectricFieldSmoothed;
                        YaxisEF = xzZmn .* xzSphCorAbsElectricFieldSmoothed;
                        
                        global optComplexPhaseFactor;
                        % Plot of electric pattern in the xz-plane (Phi = 0 or 180 deg, therefore v = 0) calculated in polar system
                        for Counter = 1 : 1 : SmoothingFactor * elementNo
                            % Calculate u and v values
                            uSphCor = k0 * d * xzXmn(Counter);
                            vSphCor = k0 * d * xzYmn(Counter);
                            % Calculate electrical field
                            xzSphCorAbsElectricFieldSmoothedOpt(Counter) = ElecFieldFun(optComplexPhaseFactor, [uSphCor, vSphCor], Lx, Ly);
                        end
                        % Scale Xmn, Ymn, and Zmn
                        XaxisOptEF = xzXmn .* xzSphCorAbsElectricFieldSmoothedOpt;
                        YaxisOptEF = xzZmn .* xzSphCorAbsElectricFieldSmoothedOpt;
                end
    end

    % Plot
    if PatternType == "All"
        % Normalize figure factor
        NormalizingFactorCRP = max(sqrt(power(XaxisCRP, 2) + power(YaxisCRP, 2)));
        % Normalize figure factor
        NormalizingFactorEF = max(sqrt(power(XaxisEF, 2) + power(YaxisEF, 2)));
        % Normalize figure factor
        NormalizingFactorOptEF = max(sqrt(power(XaxisOptEF, 2) + power(YaxisOptEF, 2)));
        
        % Normalize figure factor
        [NormalizingFactor, MaxLocation] = max([NormalizingFactorCRP, NormalizingFactorEF, NormalizingFactorOptEF]);
        
        % Plot 3D in spherical coordinate
        Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
        switch Plane
            case {'xy-plane', 'yx-plane'}
                plot(OuterPolarXBoarder * NormalizingFactor, OuterPolarYBoarder * NormalizingFactor, 'Color', 'black', 'LineWidth', 1.5);
                hold on;
                plot(InnerPolarX * NormalizingFactor, InnerPolarY * NormalizingFactor, 'Color', '#77AC30');
                % Plot lines
                plot(PolarX0 * NormalizingFactor, PolarY0 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX30 * NormalizingFactor, PolarY30 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX60 * NormalizingFactor, PolarY60 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX90 * NormalizingFactor, PolarY90 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX120 * NormalizingFactor, PolarY120 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX150 * NormalizingFactor, PolarY150 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                % Plot center point
                scatter(0, 0, 'filled');
                xlim([-inf inf]);
                ylim([-inf inf]);
                % Enhancing the surface characteristics
                axis tight;
                axis equal;
                axis([-inf inf -inf inf]);
                set(gca, 'XColor', 'none', 'YColor', 'none');
                text(NormalizingFactor *  cos(0 * pi / 180) + 0.05 * NormalizingFactor, NormalizingFactor *  sin(0 * pi / 180) + 0.01 * NormalizingFactor, '\phi = 0^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(180 * pi / 180) - 0.15 * NormalizingFactor, NormalizingFactor *  sin(180 * pi / 180) + 0.01 * NormalizingFactor, '180^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(30 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(30 * pi / 180) + 0.02 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(210 * pi / 180) - 0.15 * NormalizingFactor, NormalizingFactor *  sin(210 * pi / 180) - 0.01 * NormalizingFactor, '210^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(60 * pi / 180), NormalizingFactor *  sin(60 * pi / 180) + 0.045 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(240 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(240 * pi / 180) - 0.06 * NormalizingFactor, '240^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(90 * pi / 180) - 0.03 * NormalizingFactor, NormalizingFactor *  sin(90 * pi / 180) + 0.05 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(270 * pi / 180) - 0.06 * NormalizingFactor, NormalizingFactor *  sin(270 * pi / 180) - 0.05 * NormalizingFactor, '270^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(120 * pi / 180) - 0.07 * NormalizingFactor, NormalizingFactor *  sin(120 * pi / 180) + 0.06 * NormalizingFactor, '120^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(300 * pi / 180) - 0.05 * NormalizingFactor, NormalizingFactor *  sin(300 * pi / 180) - 0.06 * NormalizingFactor, '300^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(150 * pi / 180) - 0.13 * NormalizingFactor, NormalizingFactor *  sin(150 * pi / 180) + 0.025 * NormalizingFactor, '150^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(330 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(330 * pi / 180), '330^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.1 * NormalizingFactor *  cos(0 * pi / 180), 0.1 * NormalizingFactor *  sin(0 * pi / 180), '0.2', 'Color', '#77AC30', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.3 * NormalizingFactor *  cos(0 * pi / 180), 0.3 * NormalizingFactor *  sin(0 * pi / 180), '0.4', 'Color', '#77AC30', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.5 * NormalizingFactor *  cos(0 * pi / 180), 0.5 * NormalizingFactor *  sin(0 * pi / 180), '0.6', 'Color', '#77AC30', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.7 * NormalizingFactor *  cos(0 * pi / 180), 0.7 * NormalizingFactor *  sin(0 * pi / 180), '0.8', 'Color', '#77AC30', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.9 * NormalizingFactor *  cos(0 * pi / 180), 0.9 * NormalizingFactor *  sin(0 * pi / 180), '1.0', 'Color', '#77AC30', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                
                disp(append('Project pattern with ', num2str(NormalizingFactorCRP), ' Max. in the main direction/s ', titleText));
                disp(append('EF pattern with ', num2str(NormalizingFactorEF), ' Max. in the main direction/s ', titleText));
                disp(append('Optimized EF beam with ', num2str(NormalizingFactorOptEF), ' Max. in the main direction/s ', titleText));
                
                % Setting the legend
                TextCRP = append("Normalized project pattern");
                TextEF = append("Normalized electric field pattern before optimization");
                TextOptEF = append("Normalized electric field pattern after optimization");
                
                % Adjustment of the magnitudes with new scale
                Plot2DCRP = plot(XaxisCRP / NormalizingFactorCRP * NormalizingFactor, YaxisCRP / NormalizingFactorCRP * NormalizingFactor, 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--', 'DisplayName', TextCRP);
                Plot2DEF = plot(XaxisEF / NormalizingFactorEF * NormalizingFactor, YaxisEF / NormalizingFactorEF * NormalizingFactor, 'LineWidth', 2, 'Color', 'magenta', 'LineStyle', ':', 'DisplayName', TextEF);
                Plot2DOptEF = plot(XaxisOptEF / NormalizingFactorOptEF * NormalizingFactor, YaxisOptEF / NormalizingFactorOptEF * NormalizingFactor, 'LineWidth', 2, 'Color', 'blue', 'DisplayName', TextOptEF);
                text(1.05 * NormalizingFactor *  cos(45 * pi / 180), 1.05 * NormalizingFactor *  sin(45 * pi / 180), '\theta = 90^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                legend([Plot2DCRP, Plot2DEF, Plot2DOptEF], 'FontWeight', 'bold', 'FontSize' , 16, 'FontName' , 'Times New Roman', 'Location', 'best');
                hold off;
            case {'yz-plane', 'zy-plane'}
                plot(OuterPolarXBoarder(1 : Step / 2) * NormalizingFactor, OuterPolarYBoarder(1 : Step / 2) * NormalizingFactor, 'Color', 'black', 'LineWidth', 1.5);
                hold on;
                plot(InnerPolarX(1 : Step / 2, :) * NormalizingFactor, InnerPolarY(1 : Step / 2, :) * NormalizingFactor, 'Color', '#77AC30');
                % Plot lines
                plot(PolarX0 * NormalizingFactor, PolarY0 * NormalizingFactor, 'Color', 'black', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX30(Step / 2 : end) * NormalizingFactor, PolarY30(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX60(Step / 2 : end) * NormalizingFactor, PolarY60(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX90(Step / 2 : end) * NormalizingFactor, PolarY90(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX120(Step / 2 : end) * NormalizingFactor, PolarY120(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX150(Step / 2 : end) * NormalizingFactor, PolarY150(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                % Plot center point
                scatter(0, 0, 'filled');
                xlim([-inf inf]);
                ylim([-inf inf]);
                % Enhancing the surface characteristics
                axis tight;
                axis equal;
                axis([-inf inf -inf inf]);
                set(gca, 'XColor', 'none', 'YColor', 'none');
                text(NormalizingFactor *  cos(0 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(0 * pi / 180) + 0.01 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(180 * pi / 180) - 0.11 * NormalizingFactor, NormalizingFactor *  sin(180 * pi / 180) + 0.01 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(30 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(30 * pi / 180) + 0.02 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(210 * pi / 180) - 0.15 * NormalizingFactor, NormalizingFactor *  sin(210 * pi / 180) - 0.01 * NormalizingFactor, '210^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(60 * pi / 180), NormalizingFactor *  sin(60 * pi / 180) + 0.045 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(240 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(240 * pi / 180) - 0.06 * NormalizingFactor, '240^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(90 * pi / 180) - 0.065 * NormalizingFactor, NormalizingFactor *  sin(90 * pi / 180) + 0.05 * NormalizingFactor, '\theta = 0^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(270 * pi / 180) - 0.06 * NormalizingFactor, NormalizingFactor *  sin(270 * pi / 180) - 0.05 * NormalizingFactor, '270^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(120 * pi / 180) - 0.055 * NormalizingFactor, NormalizingFactor *  sin(120 * pi / 180) + 0.055 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(300 * pi / 180) - 0.05 * NormalizingFactor, NormalizingFactor *  sin(300 * pi / 180) - 0.06 * NormalizingFactor, '300^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(150 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(150 * pi / 180) + 0.03 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(330 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(330 * pi / 180), '330^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.1 * NormalizingFactor *  cos(0 * pi / 180), 0.1 * NormalizingFactor *  sin(0 * pi / 180), '0.2', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.3 * NormalizingFactor *  cos(0 * pi / 180), 0.3 * NormalizingFactor *  sin(0 * pi / 180), '0.4', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.5 * NormalizingFactor *  cos(0 * pi / 180), 0.5 * NormalizingFactor *  sin(0 * pi / 180), '0.6', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.7 * NormalizingFactor *  cos(0 * pi / 180), 0.7 * NormalizingFactor *  sin(0 * pi / 180), '0.8', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.9 * NormalizingFactor *  cos(0 * pi / 180), 0.6 * NormalizingFactor *  sin(0 * pi / 180), '1.0', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.05 * NormalizingFactor *  cos(45 * pi / 180), 1.05 * NormalizingFactor *  sin(45 * pi / 180), '\phi = 90^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.2 * NormalizingFactor *  cos(140 * pi / 180), 1.2 * NormalizingFactor *  sin(140 * pi / 180), '\phi = 270^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                
                disp(append('Project pattern with ', num2str(NormalizingFactorCRP), ' Max. in the main direction/s ', titleText));
                disp(append('EF pattern with ', num2str(NormalizingFactorEF), ' Max. in the main direction/s ', titleText));
                disp(append('Optimized EF beam with ', num2str(NormalizingFactorOptEF), ' Max. in the main direction/s ', titleText));
                
                % Setting the legend
                TextCRP = append("Normalized project pattern");
                TextEF = append("Normalized electric field pattern before optimization");
                TextOptEF = append("Normalized electric field pattern after optimization");
                
                % Adjustment of the magnitudes with new scale
                Plot2DCRP = plot(XaxisCRP / NormalizingFactorCRP * NormalizingFactor, YaxisCRP / NormalizingFactorCRP * NormalizingFactor, 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--', 'DisplayName', TextCRP);
                Plot2DEF = plot(XaxisEF / NormalizingFactorEF * NormalizingFactor, YaxisEF / NormalizingFactorEF * NormalizingFactor, 'LineWidth', 2, 'Color', 'magenta', 'LineStyle', ':', 'DisplayName', TextEF);
                Plot2DOptEF = plot(XaxisOptEF / NormalizingFactorOptEF * NormalizingFactor, YaxisOptEF / NormalizingFactorOptEF * NormalizingFactor, 'LineWidth', 2, 'Color', 'blue', 'DisplayName', TextOptEF);
                legend([Plot2DCRP, Plot2DEF, Plot2DOptEF], 'FontWeight', 'bold', 'FontSize' , 16, 'FontName' , 'Times New Roman', 'Location', 'south');
                hold off;
            case {'xz-plane', 'zx-plane'}
                plot(OuterPolarXBoarder(1 : Step / 2) * NormalizingFactor, OuterPolarYBoarder(1 : Step / 2) * NormalizingFactor, 'Color', 'black', 'LineWidth', 1.5);
                hold on;
                plot(InnerPolarX(1 : Step / 2, :) * NormalizingFactor, InnerPolarY(1 : Step / 2, :) * NormalizingFactor, 'Color', '#77AC30');
                % Plot lines
                plot(PolarX0 * NormalizingFactor, PolarY0 * NormalizingFactor, 'Color', 'black', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX30(Step / 2 : end) * NormalizingFactor, PolarY30(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX60(Step / 2 : end) * NormalizingFactor, PolarY60(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX90(Step / 2 : end) * NormalizingFactor, PolarY90(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX120(Step / 2 : end) * NormalizingFactor, PolarY120(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX150(Step / 2 : end) * NormalizingFactor, PolarY150(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                % Plot center point
                scatter(0, 0, 'filled');
                xlim([-inf inf]);
                ylim([-inf inf]);
                % Enhancing the surface characteristics
                axis tight;
                axis equal;
                axis([-inf inf -inf inf]);
                set(gca, 'XColor', 'none', 'YColor', 'none');
                text(NormalizingFactor *  cos(0 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(0 * pi / 180) + 0.01 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(180 * pi / 180) - 0.11 * NormalizingFactor, NormalizingFactor *  sin(180 * pi / 180) + 0.01 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(30 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(30 * pi / 180) + 0.02 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(210 * pi / 180) - 0.15 * NormalizingFactor, NormalizingFactor *  sin(210 * pi / 180) - 0.01 * NormalizingFactor, '210^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(60 * pi / 180), NormalizingFactor *  sin(60 * pi / 180) + 0.045 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(240 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(240 * pi / 180) - 0.06 * NormalizingFactor, '240^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(90 * pi / 180) - 0.065 * NormalizingFactor, NormalizingFactor *  sin(90 * pi / 180) + 0.05 * NormalizingFactor, '\theta = 0^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(270 * pi / 180) - 0.06 * NormalizingFactor, NormalizingFactor *  sin(270 * pi / 180) - 0.05 * NormalizingFactor, '270^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(120 * pi / 180) - 0.055 * NormalizingFactor, NormalizingFactor *  sin(120 * pi / 180) + 0.055 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(300 * pi / 180) - 0.05 * NormalizingFactor, NormalizingFactor *  sin(300 * pi / 180) - 0.06 * NormalizingFactor, '300^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(150 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(150 * pi / 180) + 0.03 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(330 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(330 * pi / 180), '330^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.1 * NormalizingFactor *  cos(0 * pi / 180), 0.1 * NormalizingFactor *  sin(0 * pi / 180), '0.2', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.3 * NormalizingFactor *  cos(0 * pi / 180), 0.3 * NormalizingFactor *  sin(0 * pi / 180), '0.4', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.5 * NormalizingFactor *  cos(0 * pi / 180), 0.5 * NormalizingFactor *  sin(0 * pi / 180), '0.6', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.7 * NormalizingFactor *  cos(0 * pi / 180), 0.7 * NormalizingFactor *  sin(0 * pi / 180), '0.8', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.9 * NormalizingFactor *  cos(0 * pi / 180), 0.9 * NormalizingFactor *  sin(0 * pi / 180), '1.0', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.05 * NormalizingFactor *  cos(45 * pi / 180), 1.05 * NormalizingFactor *  sin(45 * pi / 180), '\phi = 0^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.2 * NormalizingFactor *  cos(140 * pi / 180), 1.2 * NormalizingFactor *  sin(140 * pi / 180), '\phi = 180^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                
                disp(append('Project pattern with ', num2str(NormalizingFactorCRP), ' Max. in the main direction/s ', titleText));
                disp(append('EF pattern with ', num2str(NormalizingFactorEF), ' Max. in the main direction/s ', titleText));
                disp(append('Optimized EF beam with ', num2str(NormalizingFactorOptEF), ' Max. in the main direction/s ', titleText));
                
                % Setting the legend
                TextCRP = append("Normalized project pattern");
                TextEF = append("Normalized electric field pattern before optimization");
                TextOptEF = append("Normalized electric field pattern after optimization");
                
                % Adjustment of the magnitudes with new scale
                Plot2DCRP = plot(XaxisCRP / NormalizingFactorCRP * NormalizingFactor, YaxisCRP / NormalizingFactorCRP * NormalizingFactor, 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--', 'DisplayName', TextCRP);
                Plot2DEF = plot(XaxisEF / NormalizingFactorEF * NormalizingFactor, YaxisEF / NormalizingFactorEF * NormalizingFactor, 'LineWidth', 2, 'Color', 'magenta', 'LineStyle', ':', 'DisplayName', TextEF);
                Plot2DOptEF = plot(XaxisOptEF / NormalizingFactorOptEF * NormalizingFactor, YaxisOptEF / NormalizingFactorOptEF * NormalizingFactor, 'LineWidth', 2, 'Color', 'blue', 'DisplayName', TextOptEF);
                legend([Plot2DCRP, Plot2DEF, Plot2DOptEF], 'FontWeight', 'bold', 'FontSize' , 16, 'FontName' , 'Times New Roman', 'Location', 'south');
                hold off;
        end
        
    else
        % Normalize figure factor
        NormalizingFactor = max(sqrt(power(Xaxis, 2) + power(Yaxis, 2)));
        
        % Plot 3D in spherical coordinate
        Fig = figure('Name', FigureName, 'DefaultAxesColorOrder' , [left_color; right_color], 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold', 'Color', 'white', 'WindowState', 'maximized');
        switch Plane
            case {'xy-plane', 'yx-plane'}
                plot(OuterPolarXBoarder * NormalizingFactor, OuterPolarYBoarder * NormalizingFactor, 'Color', 'black', 'LineWidth', 1.5);
                hold on;
                plot(InnerPolarX * NormalizingFactor, InnerPolarY * NormalizingFactor, 'Color', 'blue');
                % Plot lines
                plot(PolarX0 * NormalizingFactor, PolarY0 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX30 * NormalizingFactor, PolarY30 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX60 * NormalizingFactor, PolarY60 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX90 * NormalizingFactor, PolarY90 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX120 * NormalizingFactor, PolarY120 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX150 * NormalizingFactor, PolarY150 * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                % Plot center point
                scatter(0, 0, 'filled');
                Plot2D = plot(Xaxis, Yaxis, 'LineWidth', 2, 'Color', 'red', 'DisplayName', Text);
                xlim([-inf inf]);
                ylim([-inf inf]);
                % Enhancing the surface characteristics
                axis tight;
                axis equal;
                axis([-inf inf -inf inf]);
                set(gca, 'XColor', 'none', 'YColor', 'none');
                text(NormalizingFactor *  cos(0 * pi / 180) + 0.05 * NormalizingFactor, NormalizingFactor *  sin(0 * pi / 180) + 0.01 * NormalizingFactor, '\phi = 0^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(180 * pi / 180) - 0.15 * NormalizingFactor, NormalizingFactor *  sin(180 * pi / 180) + 0.01 * NormalizingFactor, '180^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(30 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(30 * pi / 180) + 0.02 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(210 * pi / 180) - 0.15 * NormalizingFactor, NormalizingFactor *  sin(210 * pi / 180) - 0.01 * NormalizingFactor, '210^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(60 * pi / 180), NormalizingFactor *  sin(60 * pi / 180) + 0.045 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(240 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(240 * pi / 180) - 0.06 * NormalizingFactor, '240^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(90 * pi / 180) - 0.03 * NormalizingFactor, NormalizingFactor *  sin(90 * pi / 180) + 0.05 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(270 * pi / 180) - 0.06 * NormalizingFactor, NormalizingFactor *  sin(270 * pi / 180) - 0.05 * NormalizingFactor, '270^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(120 * pi / 180) - 0.07 * NormalizingFactor, NormalizingFactor *  sin(120 * pi / 180) + 0.06 * NormalizingFactor, '120^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(300 * pi / 180) - 0.05 * NormalizingFactor, NormalizingFactor *  sin(300 * pi / 180) - 0.06 * NormalizingFactor, '300^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(150 * pi / 180) - 0.13 * NormalizingFactor, NormalizingFactor *  sin(150 * pi / 180) + 0.025 * NormalizingFactor, '150^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(330 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(330 * pi / 180), '330^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.1 * NormalizingFactor *  cos(0 * pi / 180), 0.1 * NormalizingFactor *  sin(0 * pi / 180), '0.2', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.3 * NormalizingFactor *  cos(0 * pi / 180), 0.3 * NormalizingFactor *  sin(0 * pi / 180), '0.4', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.5 * NormalizingFactor *  cos(0 * pi / 180), 0.5 * NormalizingFactor *  sin(0 * pi / 180), '0.6', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.7 * NormalizingFactor *  cos(0 * pi / 180), 0.7 * NormalizingFactor *  sin(0 * pi / 180), '0.8', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.9 * NormalizingFactor *  cos(0 * pi / 180), 0.9 * NormalizingFactor *  sin(0 * pi / 180), '1.0', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.05 * NormalizingFactor *  cos(45 * pi / 180), 1.05 * NormalizingFactor *  sin(45 * pi / 180), '\theta = 90^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                legend(Plot2D, 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'best');
                hold off;
            case {'yz-plane', 'zy-plane'}
                plot(OuterPolarXBoarder(1 : Step / 2) * NormalizingFactor, OuterPolarYBoarder(1 : Step / 2) * NormalizingFactor, 'Color', 'black', 'LineWidth', 1.5);
                hold on;
                plot(InnerPolarX(1 : Step / 2, :) * NormalizingFactor, InnerPolarY(1 : Step / 2, :) * NormalizingFactor, 'Color', 'blue');
                % Plot lines
                plot(PolarX0 * NormalizingFactor, PolarY0 * NormalizingFactor, 'Color', 'black', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX30(Step / 2 : end) * NormalizingFactor, PolarY30(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX60(Step / 2 : end) * NormalizingFactor, PolarY60(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX90(Step / 2 : end) * NormalizingFactor, PolarY90(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX120(Step / 2 : end) * NormalizingFactor, PolarY120(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX150(Step / 2 : end) * NormalizingFactor, PolarY150(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                % Plot center point
                scatter(0, 0, 'filled');
                Plot2D = plot(Xaxis, Yaxis, 'LineWidth', 2, 'Color', 'red', 'DisplayName', Text);
                xlim([-inf inf]);
                ylim([-inf inf]);
                % Enhancing the surface characteristics
                axis tight;
                axis equal;
                axis([-inf inf -inf inf]);
                set(gca, 'XColor', 'none', 'YColor', 'none');
                text(NormalizingFactor *  cos(0 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(0 * pi / 180) + 0.01 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(180 * pi / 180) - 0.11 * NormalizingFactor, NormalizingFactor *  sin(180 * pi / 180) + 0.01 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(30 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(30 * pi / 180) + 0.02 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(210 * pi / 180) - 0.15 * NormalizingFactor, NormalizingFactor *  sin(210 * pi / 180) - 0.01 * NormalizingFactor, '210^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(60 * pi / 180), NormalizingFactor *  sin(60 * pi / 180) + 0.045 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(240 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(240 * pi / 180) - 0.06 * NormalizingFactor, '240^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(90 * pi / 180) - 0.065 * NormalizingFactor, NormalizingFactor *  sin(90 * pi / 180) + 0.05 * NormalizingFactor, '\theta = 0^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(270 * pi / 180) - 0.06 * NormalizingFactor, NormalizingFactor *  sin(270 * pi / 180) - 0.05 * NormalizingFactor, '270^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(120 * pi / 180) - 0.055 * NormalizingFactor, NormalizingFactor *  sin(120 * pi / 180) + 0.055 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(300 * pi / 180) - 0.05 * NormalizingFactor, NormalizingFactor *  sin(300 * pi / 180) - 0.06 * NormalizingFactor, '300^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(150 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(150 * pi / 180) + 0.03 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(330 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(330 * pi / 180), '330^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.1 * NormalizingFactor *  cos(0 * pi / 180), 0.1 * NormalizingFactor *  sin(0 * pi / 180), '0.2', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.3 * NormalizingFactor *  cos(0 * pi / 180), 0.3 * NormalizingFactor *  sin(0 * pi / 180), '0.4', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.5 * NormalizingFactor *  cos(0 * pi / 180), 0.5 * NormalizingFactor *  sin(0 * pi / 180), '0.6', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.7 * NormalizingFactor *  cos(0 * pi / 180), 0.7 * NormalizingFactor *  sin(0 * pi / 180), '0.8', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.9 * NormalizingFactor *  cos(0 * pi / 180), 0.6 * NormalizingFactor *  sin(0 * pi / 180), '1.0', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.05 * NormalizingFactor *  cos(45 * pi / 180), 1.05 * NormalizingFactor *  sin(45 * pi / 180), '\phi = 90^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.2 * NormalizingFactor *  cos(140 * pi / 180), 1.2 * NormalizingFactor *  sin(140 * pi / 180), '\phi = 270^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                legend(Plot2D, 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'south');
            case {'xz-plane', 'zx-plane'}
                plot(OuterPolarXBoarder(1 : Step / 2) * NormalizingFactor, OuterPolarYBoarder(1 : Step / 2) * NormalizingFactor, 'Color', 'black', 'LineWidth', 1.5);
                hold on;
                plot(InnerPolarX(1 : Step / 2, :) * NormalizingFactor, InnerPolarY(1 : Step / 2, :) * NormalizingFactor, 'Color', 'blue');
                % Plot lines
                plot(PolarX0 * NormalizingFactor, PolarY0 * NormalizingFactor, 'Color', 'black', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX30(Step / 2 : end) * NormalizingFactor, PolarY30(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX60(Step / 2 : end) * NormalizingFactor, PolarY60(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX90(Step / 2 : end) * NormalizingFactor, PolarY90(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX120(Step / 2 : end) * NormalizingFactor, PolarY120(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                plot(PolarX150(Step / 2 : end) * NormalizingFactor, PolarY150(Step / 2 : end) * NormalizingFactor, 'Color', '#4DBEEE', 'LineWidth', 1, 'LineStyle', '-');
                % Plot center point
                scatter(0, 0, 'filled');
                Plot2D = plot(Xaxis, Yaxis, 'LineWidth', 2, 'Color', 'red', 'DisplayName', Text);
                xlim([-inf inf]);
                ylim([-inf inf]);
                % Enhancing the surface characteristics
                axis tight;
                axis equal;
                axis([-inf inf -inf inf]);
                set(gca, 'XColor', 'none', 'YColor', 'none');
                text(NormalizingFactor *  cos(0 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(0 * pi / 180) + 0.01 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(180 * pi / 180) - 0.11 * NormalizingFactor, NormalizingFactor *  sin(180 * pi / 180) + 0.01 * NormalizingFactor, '90^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(30 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(30 * pi / 180) + 0.02 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(210 * pi / 180) - 0.15 * NormalizingFactor, NormalizingFactor *  sin(210 * pi / 180) - 0.01 * NormalizingFactor, '210^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(60 * pi / 180), NormalizingFactor *  sin(60 * pi / 180) + 0.045 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(240 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(240 * pi / 180) - 0.06 * NormalizingFactor, '240^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(90 * pi / 180) - 0.065 * NormalizingFactor, NormalizingFactor *  sin(90 * pi / 180) + 0.05 * NormalizingFactor, '\theta = 0^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(270 * pi / 180) - 0.06 * NormalizingFactor, NormalizingFactor *  sin(270 * pi / 180) - 0.05 * NormalizingFactor, '270^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(120 * pi / 180) - 0.055 * NormalizingFactor, NormalizingFactor *  sin(120 * pi / 180) + 0.055 * NormalizingFactor, '30^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(300 * pi / 180) - 0.05 * NormalizingFactor, NormalizingFactor *  sin(300 * pi / 180) - 0.06 * NormalizingFactor, '300^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(NormalizingFactor *  cos(150 * pi / 180) - 0.08 * NormalizingFactor, NormalizingFactor *  sin(150 * pi / 180) + 0.03 * NormalizingFactor, '60^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
%                 text(NormalizingFactor *  cos(330 * pi / 180) + 0.02 * NormalizingFactor, NormalizingFactor *  sin(330 * pi / 180), '330^\circ', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.1 * NormalizingFactor *  cos(0 * pi / 180), 0.1 * NormalizingFactor *  sin(0 * pi / 180), '0.2', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.3 * NormalizingFactor *  cos(0 * pi / 180), 0.3 * NormalizingFactor *  sin(0 * pi / 180), '0.4', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.5 * NormalizingFactor *  cos(0 * pi / 180), 0.5 * NormalizingFactor *  sin(0 * pi / 180), '0.6', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.7 * NormalizingFactor *  cos(0 * pi / 180), 0.7 * NormalizingFactor *  sin(0 * pi / 180), '0.8', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(0.9 * NormalizingFactor *  cos(0 * pi / 180), 0.9 * NormalizingFactor *  sin(0 * pi / 180), '1.0', 'Color', 'blue', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.05 * NormalizingFactor *  cos(45 * pi / 180), 1.05 * NormalizingFactor *  sin(45 * pi / 180), '\phi = 0^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                text(1.2 * NormalizingFactor *  cos(140 * pi / 180), 1.2 * NormalizingFactor *  sin(140 * pi / 180), '\phi = 180^\circ', 'Color', 'magenta', 'FontSize', 18, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontWeight', 'bold');
                legend(Plot2D, 'FontWeight', 'bold', 'FontSize' , 20, 'FontName' , 'Times New Roman', 'Location', 'south');
        end
    end
end
