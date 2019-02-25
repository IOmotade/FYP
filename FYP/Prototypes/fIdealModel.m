
%% simple model
N = 3; 
R = 10e3; R = repmat(R, [N, N]);
vs = 5; vs = repmat(vs, [N, 1]);
G = 1./R;
I = transpose(G)*vs;
V_ = fUnits(vs, 'V')
R_ = fUnits(R, 'Ohm')
I_ = fUnits(transpose(G)*vs, 'A')