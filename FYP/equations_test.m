%% simple model
N = 3;
R = 10e3; R = repmat(R, [N, N]);
vs = 5; vs = repmat(vs, [N, 1]);
G = 1./R;
I = transpose(G)*vs;
V_ = fUnits(vs, 'V')
R_ = fUnits(R, 'Ohm')
I_ = fUnits(transpose(G)*vs, 'A')

%% second model
N = 3;
%define array
R = 10e3; R = repmat(R, [N, N]); R(3, 3) = 10e3;
% Rc = R; Rr = R;
Rc = 1e3*ones(size(R)); Rr = 1e3*ones(size(R));
%define operation matrices
J = fShiftingMatrix(N);
U = fUTriangularMatrix(N);
%define signals
%voltages
vs = 5; vs = repmat(vs, [N, 1]);
Vs = [vs zeros(N, N-1)];
vc = zeros(N, 1);
Vc = ones(N, 1)*transpose(vc);
Vo = ones(N, 1)*transpose(vs);
%currents
Is = zeros(N);
ivs = [0; 0; 0];
Ivs = [ivs zeros(N, N-1)];
Io = zeros(size(Vo));

%perform calculations
G = 1./R;
Vo = (Vo*J' + Vs) - Rr.*(Is - Io + (Ivs + (Io*J' + Ivs)));
Io = (G.*...
    ((Vo*J' + Vs) - (Vc+U*(Rc.*Io))) + Is) ...
    /(eye(N) - J);
ivs = sum(Io, 2);
Ivs = [ivs zeros(N, N-1)];
%display results
Vo_ = fUnits(Vo, 'V')
Io_ = fUnits(Io, 'A')

%% second model version 2
rst
N = 3;
%define array
R = 10e3; R = repmat(R, [N, N]); %R(3, 3) = 50;% R(3, 2) = 5e3; R(3, 3) = 30e3;
% R = 10e3*abs(randn(size(R)))
G = 1./R;
% Rc = R; Rr = R;
Rc = 1e3*ones(size(R)); Rr = 1e3*ones(size(R));
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
    Val = - Rr.*(Is - Io + Ii);
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
total_iter = total_iter + iter
%display results
total_iter_ = fUnits(total_iter, 'iterations')
Vo_ = fUnits(Vo, 'V')
Io_prev_ = fUnits(Io_prev, 'A')
Io_ = fUnits(Io, 'A')



%% second model version 3
rst
N = 3;
%define array
R = 10e3; R = repmat(R, [N, N]);% R(3, 3) = 1.2e3;% R(3, 2) = 5e3; R(3, 3) = 30e3;
% R = 10e3*abs(randn(size(R)))
G = 1./R;
% Rc = R; Rr = R;
Rc = 1e3*ones(size(R)); Rr = 1e3*ones(size(R));
%define operation matrices
In = eye(N);
J = fShiftingMatrix(N);
U = fUTriangularMatrix(N);
A = In - J;
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
iter_limit = 5000;
% %% Give answer:
% Io(N,:) = [1.03249e-3 1.03271e-3 2.7e-3];

% %% calculations for second model v2
iter = 0;
hasConverged = false;
while ((~hasConverged) && (iter<iter_limit))||(iter<5)
    Io_prev = Io;
    Vo_prev = Vo;
    Val = - Rr.*(Is - Io + Ii);
    Vo = Vi - Rr.*(Is - Io + Ii);
    % Io = G.*(Vi - X + (R.*Ii));
    Io = G.*(Vi - X) + Ii;
    iter = iter + 1;
    %recalculated inputs
    X = Vc + U*(Rc.*Io);        % with x method 2 (Vc has all non-zero elements) better!
    % X = (J'*X) + Vc + Rc.*Io; % with x method 2 (Vc has zero elements)
    Ii = J*Io;
    Vi = Vo*J' + Vs;
    % Is = Io - Ii + (Is*J);            %1st approx method for Is
    % Is = (Io - Ii)*(eye(N) - J)^(-1);   %2nd approx method for Is (faster!, better!)
    % Is = (eye(N) - J)*(Io)*(eye(N) - J)^(-1);   %3rd approx method for Is
    Is = A*(Io)*(A^(-1));   %4th approx method for Is (more compact)
    %check for convergence
    hasConverged = fHasConverged([Io Vo], [Io_prev Vo_prev], 1e-1);
    [Io Vo]
end
total_iter = total_iter + iter
%display results
total_iter_ = fUnits(total_iter, 'iterations')
Vo_ = fUnits(Vo, 'V')
Io_prev_ = fUnits(Io_prev, 'A')
Io_ = fUnits(Io, 'A')

%% second model version 4
rst
N = 3;
%define array
R = 10e3; R = repmat(R, [N, N]); %R(3, 3) = 2.5e3;% R(3, 2) = 5e3; R(3, 3) = 30e3;
% R = 10e3*abs(randn(size(R)))
G = 1./R;
% Rc = R; Rr = R;
Rc = 1e3*ones(size(R)); Rr = 1e3*ones(size(R));
%define operation matrices
In = eye(N);
J = fShiftingMatrix(N);
U = fUTriangularMatrix(N);
A = In - J;
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
%program limits
iter = 0
total_iter = 0;
iter_limit = 5000;
% %% Give answer:
% Io(N,:) = [1.03249e-3 1.03271e-3 2.7e-3];

% %% calculations for second model v2
iter = 0;
hasConverged = false;
while ((~hasConverged) && (iter<iter_limit))||(iter<5)
    Io_prev = Io;
    Vo_prev = Vo;
    % Io = G.*(Vi - X) + Ii;
    Vo = (Vs - Rr.*(A*Io*((A^-1) - In))) * ((In - J')^(-1));
    Io = (A^(-1))*(G.*((Vo*J') + Vs - Vc - (U*(Rc.*Io))));
    iter = iter + 1;
    %check for convergence
    hasConverged = fHasConverged([Io Vo], [Io_prev Vo_prev], 1e-1);
    [Io Vo];
end
total_iter = total_iter + iter
%display results
total_iter_ = fUnits(total_iter, 'iterations')
Vo_ = fUnits(Vo, 'V')
Io_prev_ = fUnits(Io_prev, 'A')
Io_ = fUnits(Io, 'A')

%% second model version 5 with Vs slowly increasing
tic
rst
N = 3;
%define array
R = 10e3; R = repmat(R, [N, N]); R(3, 3) = 3e3;% R(3, 2) = 5e3; R(3, 3) = 30e3;
% R = 10e3*abs(randn(size(R)))
G = 1./R;
% Rc = R; Rr = R;
Rc = 1e-9*ones(size(R)); Rr = 1e-9*ones(size(R));
Rc_f = 1e3*ones(size(R)); Rr_f = 1e3*ones(size(R));
%define operation matrices
In = eye(N);
J = fShiftingMatrix(N);
U = fUTriangularMatrix(N);
A = In - J;
%define signals
%voltages
vs = 5; vs = repmat(vs, [N, 1]);
% Vs = [zeros(size(vs)) zeros(N, N-1)];
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
%program limits
iter = 0
total_iter = 0;
iter_limit = 1000;
% %% Give answer:
% Io(N,:) = [1.03249e-3 1.03271e-3 2.7e-3];

% %% calculations for second model v2
while ~(fHasConverged(Rc, Rc_f, 1e-6) && fHasConverged(Rr, Rr_f, 1e-6))
    iter = 0;
    hasConverged = false;
    while ((~hasConverged) && (iter<iter_limit))||(iter<5)
        Io_prev = Io;
        Vo_prev = Vo;
        % Io = G.*(Vi - X) + Ii;
        Vo = (Vs - Rr.*(A*Io*((A^-1) - In))) * ((In - J')^(-1));
        Io = (A^(-1))*(G.*((Vo*J') + Vs - Vc - (U*(Rc.*Io))));
        iter = iter + 1;
        %check for convergence
        hasConverged = fHasConverged([Io Vo], [Io_prev Vo_prev], 1e-6);
        [Io Vo];
    end
    Rr = Rr + 1e-3*Rr_f;
    Rc = Rc + 1e-3*Rc_f;
    if(~hasConverged || any(any(isnan(Io))) || any(any(Vo>5)))
        hasConverged;
        iter_limit = iter_limit + 1000;
    end
    total_iter = total_iter + iter;
end
%display results
total_iter_ = fUnits(total_iter, 'iterations')
Vo_ = fUnits(Vo, 'V')
Io_prev_ = fUnits(Io_prev, 'A')
Io_ = fUnits(Io, 'A')
toc

%% second model version 6 with Rr & Rc slowly increasing
tic
rst
N = 3;
%define array
R = 100e3; R = repmat(R, [N, N]); %R(3, 3) = 10e3;% R(3, 2) = 1e3; R(3, 3) = 30e3;
% R = 10e3*abs(randn(size(R)))
G = 1./R;
Rc = 1e-9*ones(size(R)); Rr = 1e-9*ones(size(R));
Rc_f = 10e3*ones(size(R)); Rr_f = 10e3*ones(size(R));
%define operation matrices
In = eye(N);
J = fShiftingMatrix(N);
U = fUTriangularMatrix(N);
A = In - J;
%define signals
%voltages
vs = 5; vs = repmat(vs, [N, 1]);
Vs = [zeros(size(vs)) zeros(N, N-1)];
Vs_f = [vs zeros(N, N-1)];
vc = zeros(N, 1);
Vc = ones(N, 1)*transpose(vc);
Vo = ones(N, 1)*transpose(vs); %Vo = [Vo(:, 1), zeros(N, N-1)]
Vi = Vo*J' + Vs;
X = Vc;
%currents
Is = zeros(N);
Io = zeros(N);
Ii = J*Io;
%program limits
iter = 0
total_iter = 0;
iter_limit = 1000;
stop_now=false;
% %% Give answer:
% Io(N,:) = [1.03249e-3 1.03271e-3 2.7e-3];

% %% calculations for second model v2
while (~(fHasConverged(Rc, Rc_f, 1e-6) && fHasConverged(Rr, Rr_f, 1e-6)) &&...
        ~stop_now)
    iter = 0;
    hasConverged = false;
    while ((~hasConverged) && (iter<iter_limit))||(iter<5)
        Io_prev = Io;
        Vo_prev = Vo;
        Val = - Rr.*(Is - Io + Ii);
        Vo = Vi - Rr.*(Is - Io + Ii);
        % Io = G.*(Vi - X + (R.*Ii));
        Io = G.*(Vi - X) + Ii;
        iter = iter + 1;
        %recalculated inputs
        X = Vc + U*(Rc.*Io);        % with x method 2 (Vc has all non-zero elements) better!
        % X = (J'*X) + Vc + Rc.*Io; % with x method 2 (Vc has zero elements)
        Ii = J*Io;
        Vi = Vo*J' + Vs;
        % Is = Io - Ii + (Is*J);            %1st approx method for Is
        % Is = (Io - Ii)*(eye(N) - J)^(-1);   %2nd approx method for Is (faster!, better!)
        % Is = (eye(N) - J)*(Io)*(eye(N) - J)^(-1);   %3rd approx method for Is
        Is = A*(Io)*(A^(-1));   %4th approx method for Is (more compact)
        %check for convergence
        hasConverged = fHasConverged([Io Vo], [Io_prev Vo_prev], 1e-1);
        [Io Vo];
    end
    Rr = Rr + 1e-5*Rr_f;
    Rc = Rc + 1e-5*Rc_f;
    if ~fHasConverged(Vs, Vs_f)
        Vs = Vs + 1e-3*Vs_f;
    end
    if(~hasConverged)
        hasConverged
        iter_limit = iter_limit + 1000;
        if(any(any(isnan(Io))) || any(any(Vo>5)))
            stop_now=true;
        end
    end
    total_iter = total_iter + iter;
end
%display results
total_iter_ = fUnits(total_iter, 'iterations')
Vo_ = fUnits(Vo, 'V')
Io_prev_ = fUnits(Io_prev, 'A')
Io_ = fUnits(Io, 'A')
toc