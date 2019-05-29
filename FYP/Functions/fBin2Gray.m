% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts Binary number to its Gray code equivalent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% b (1xN Char) = Character array for input binary number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% g (1xN Char) = Character array for output gray code equivalent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = fBin2Gray(b)
[r, ~] = size(b);

g = char(zeros(size(b)));
for idx = 1:r
    g(idx, :)  = internalBin2Gray(b(idx, :));
end

end

function g = internalBin2Gray(b)
g(1) = b(1);
for i = 2 : length(b);
    x = xor(str2num(b(i-1)), str2num(b(i)));
    g(i) = num2str(x);
end
end