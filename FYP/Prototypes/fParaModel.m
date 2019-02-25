function [Io, Vo] = fParaModel(R, Rr, Rc, Vs, Io, Vo)
%%


end
%%
function scheme_v1(R, Rr, Rc, Vs)
%define operation matrices
J = fShiftingMatrix(N);
U = fUTriangularMatrix(N);
%define signals
%voltages
vs = 5; vs = repmat(vs, [N, 1]);
Vs = [vs zeros(N, N-1)];
vc = zeros(N, 1);
Vc = ones(N, 1)*transpose(vc);
Vo = ones(N, 1)*transpose(vs); %Vo = [Vo(:, 1), zeros(N, N-1)]
Vi = Vo*J' + Vs;
X = Vc;
%currents
Is = zeros(N);
Io = zeros(N);
Ii = J*Io;
iter = 0
total_iter = 0;
iter_limit = 2000;
% %% Give answer:
% Io(N,:) = [1.03249e-3 1.03271e-3 2.7e-3];

% %% calculations for second model v2
iter = 0;
hasConverged = false;
while ((~hasConverged) && (iter<iter_limit))||(iter<5)
    Io_prev = Io;
    Vo_prev = Vo;
    Vo = Vi - Rr.*(Is - Io + Ii);
    % Io = G.*(Vi - X + (R.*Ii));
    Io = G.*(Vi - X) + Ii;
    iter = iter + 1;
    %recalculated inputs
    X = Vc + U*(Rc.*Io);        % with x method 2 (Vc has all non-zero elements) better!
    % X = (J'*X) + Vc + Rc.*Io; % with x method 2 (Vc has zero elements)
    Ii = J*Io;
    Vi = Vo*J' + Vs;
    Is = Io - Ii + (Is*J);
    %check for convergence
    hasConverged = fHasConverged([Io Vo], [Io_prev Vo_prev], 1e-3);
    [Io Vo]
end
end