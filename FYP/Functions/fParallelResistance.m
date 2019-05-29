% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates equivalent resistance value of two resistors in parallel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% Rx, Ry (Double) = Resistance values with Rx and Ry in parallel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% Req (Double) = Equivalent resistance of circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Req = fParallelResistance(Rx, Ry)
% Req = (Rx*Ry)/(Rx + Ry);
if isinf(Rx)
    if isinf(Ry)
        Req = Inf;
    else
        Req = (Ry)/(1 + Ry/Rx);
    end
else
    Req = (Rx)/(1 + Rx/Ry);
end

end