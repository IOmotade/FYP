function i = fIdealModel(N, v, MemR)
%% simple model
G = 1./MemR;
I = transpose(G)*v;
end