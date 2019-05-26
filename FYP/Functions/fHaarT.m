% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gives the Haar Transform of a signal x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% x (1xN Double) = time series values of signal x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% hft (1xN Double) = haar transform series values of signal x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hft = fHaarT(x)

% if(mod(length(x), 2)==1)
%     x = [x, 0];
% end

hft = zeros(size(x));

t = 0:2*pi/(length(x)-1):2*pi;
for idx = 1:length(x)
    h = square((idx-1)*t);
    if sum(size(h) == size(x)) == 0
        h = transpose(h);
    end
    hft(idx) = sum(x./h);
end
hft = hft/length(hft);

end