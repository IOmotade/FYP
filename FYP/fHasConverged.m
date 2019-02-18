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