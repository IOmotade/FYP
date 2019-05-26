% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gives the Inverse Haar Transform of hft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% hft (1xN Double) = haar transform of a signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% h (1xN Double) = inverse haar transform/time series values of signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = fInvHaarT(hft)

% if(mod(length(hft), 2)==1)
%     hft = [hft, 0];
% end

h = zeros(size(hft));

t = 0:2*pi/(length(hft)-1):2*pi;
for idx = 1:length(hft)
    h = h + hft(idx)*square((idx-1)*t)/2;
%     plot(square((idx-1)*t));
%     title("Current Figure")
end
% h = h/2;
end