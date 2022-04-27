function [UPointsVec, VPointsVec] = EquiDistant2DPoints(UArray, VArray, Diameter)
    %% This function produces NP points which distributed equi-distantly on the 2D UV Coordinates
    % NP is the number of points
    % ULimit and VLimit are range of U and V axes in the UV coordinate, respectively.
    
    UMax = max(UArray);
    UMin = min(UArray);
    UStep = sqrt(3) * Diameter / 2;
    UGrid = UMin : UStep : UMax;
    ULength = numel(UGrid);
    
    VMax = max(VArray);
    VMin = min(VArray);
    VStep = Diameter;
    VGrid = VMin : VStep : VMax;
    VLength = numel(VGrid);
    
    DistributedPoints = NaN(ULength * VLength, 2);
    IndexModifier = 1;
    for UCounter = 1 : 1 : ULength
        StartingU = UMin + (UCounter - 1) * (sqrt(3) * Diameter / 2);
        % Odd UCounter
        if (mod(UCounter, 2) ~= 0)
            StartingV = VMin - Diameter;
            for VCounter = 1 : 1 : VLength
                StartingV = StartingV + Diameter;
                if StartingV > VMax
                    continue;
                end
                DistributedPoints(IndexModifier, :) = [StartingU, StartingV];
                IndexModifier = IndexModifier + 1;
            end
        % Even UCounter
        else
            StartingV = VMin + Diameter / 2;
            if StartingV > VMax
                continue;
            end
            for VCounter = 1 : 1 : VLength
                if StartingV > VMax
                    continue;
                end
                DistributedPoints(IndexModifier, :) = [StartingU, StartingV];
                IndexModifier = IndexModifier + 1;
                StartingV = StartingV + Diameter;
            end
        end
    end
    
    DistributedPointsUModification = (UMax - max(DistributedPoints(:, 1))) / 2;
    DistributedPointsVModification = (VMax - max(DistributedPoints(:, 2))) / 2; 
    
    % Relocate the points at the center of the 2D plane
    UPointsVec = DistributedPointsUModification + transpose(DistributedPoints(~isnan(DistributedPoints(:, 1)), 1));
    VPointsVec = DistributedPointsVModification + transpose(DistributedPoints(~isnan(DistributedPoints(:, 2)), 2));
end