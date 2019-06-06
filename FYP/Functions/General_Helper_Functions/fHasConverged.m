% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts a simulated circuit from the time domain to frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% x (MxN Double) = Matrix of values
% x_prev (MxN Double) = Previous state of matrix x
% tol (Double) = tolerance values i.e. measure of similarity between x and 
% x_prev (kinda vaguely defined tbh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% isConverged (Boolean) = True if the variable x has converged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isConverged = fHasConverged(x, x_prev, tol)

if ~exist('tol', 'var')
    tol = 1e-3;
end
tol = repmat(tol, size(x));

z = arrayfun(@hasC, x, x_prev, tol);
isConverged = all(z(:));

end

function isC = hasC(x, x_prev, tol)
if (x==0)
    x = x+1;
    x_prev = x_prev+1;
end
convergence = (x-x_prev)/(tol*x);
isC = (abs(convergence)<tol);

end