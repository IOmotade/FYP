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
function hft = fHaarT(x, method)

% Default method is 'clean OFDM' method
if ~exist('method', 'var')
    method = '2';
end

if method=='1'
    hft = version1(x);
elseif method=='2'
    hft = version2(x);
end

end

%% Version 1: Using Square Waves as composite frequencies
function hft = version1(x)

if(mod(length(x), 2)==1)
    x = [x, 0];
end
x = x(:);

hft = zeros(size(x));

t = 0:2*pi/(length(x)-1):2*pi;
for idx = 1:length(x)
    h = square((idx-1)*t);
    if sum(size(h) == size(x)) == 0
        h = transpose(h);
    end
    hft(idx) = sum(x(:)./h);
%     H(idx, :) = h;
end
hft = hft/length(hft);

end

%% Version 2: Directly using Transformation Matrix
function hft = version2(x)
N = length(x);
if (N<2 || (log2(N)-floor(log2(N)))~=0)
    error('The input argument should be of form 2^k');
end

H = fOFDMTransformMtrx(N);

intermediate_hft = (H*x(:))/N;

hft = [intermediate_hft(1); zeros(N-1, 1)];

for idx = 2:length(intermediate_hft)
%     idx-1
    hft_idx = (2^((idx-1)-1) + 1);% - 1
    hft(hft_idx) = intermediate_hft(idx);
end

end