% function [Gaussian3D] = modifiedGaussianFilter(U, V)
function [modifiedGaussian3D] = modifiedGaussianFilter(U, V, BeamNum)
    global L u0 v0 Lambda d;

    % Initialize the outpout of the modified Gaussian function
    modifiedGaussian3D = 0;
    
    % Standard deviation of the distribution
    k0 = 2 * pi / Lambda;
    STD = k0 * d;
    
    %% 2D customized Gaussian filter
%     Gaussian2D = @(x) 1 / sqrt(2 * pi * power(STD, 2)) .* exp(-power(x - Mean, 2) / (2 * power(STD, 2)));
    
    %% 3D modified Gaussian Filter
    if BeamNum == 0
        for Counter = 1 : 1 : L
            % Gaussian3D = @(x, y) 1 / (2 * pi * power(STD, 2)) .* exp(-(power(y - u0, 2) + power(x - v0, 2)) ./ (2 * power(STD, 2)));
            modifiedGaussian3D = modifiedGaussian3D + exp(-(power(U - u0(Counter), 2) + power(V - v0(Counter), 2)) ./ (2 * power(STD, 2)));
        end
        % Calculate the mean of modified Gaussian filter for points outside of protected area
        modifiedGaussian3D = modifiedGaussian3D / L;
    else
        % Gaussian3D = @(x, y) 1 / (2 * pi * power(STD, 2)) .* exp(-(power(y - u0, 2) + power(x - v0, 2)) ./ (2 * power(STD, 2)));
        modifiedGaussian3D = exp(-(power(U - u0(BeamNum), 2) + power(V - v0(BeamNum), 2)) ./ (2 * power(STD, 2)));
    end
end