% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates ratio of currents that flow through two resistors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% Rx, Ry (Double) = Resistor Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% currentRatioX, currentRatioY (Double) = Current ratios for both resistors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [currentRatioX, currentRatioY] = fCurrentRatio(Rx, Ry)
[currentRatioX, currentRatioY] = deal(...
    (Rx^-1)/((Rx^-1)+(Ry^-1)), (Ry^-1)/((Rx^-1)+(Ry^-1)));

end