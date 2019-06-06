function fAlert()
load('gong.mat', 'y');
N = length(y);

% w = window(@blackmanharris, N);
w = window(@gausswin, N, 5);
w = w(ceil(end/2):end);
w = [w; zeros(N-length(w), 1)];

sound(y.*w)
end