% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a unique number from an integer tuple
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% a,b (Double)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% f (Double) = Tuple's Unique Numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = fCantorTupleFunction(v)
%Cantor Tuple Function
if length(v)==1
    f = v;
    return
end
f = fPairFunction(fCantorTupleFunction(v(1:end-1)), v(end));
end
