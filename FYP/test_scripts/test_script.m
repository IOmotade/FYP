%% Simple Circuit
% Testing simple Line Resistance Estimation 
% rst
addpath(genpath('../FYP/Functions/'))
[fig_az, fig_el] = deal(145, 30);

%% Setup Circuit System
N = 8;
p = 0.5; %1e-9+(0:0.1:1);
Rp_estimate = zeros(size(p));
for p_idx = 1:length(p)
    MemR  = abs(10e3*ones(N) + 0e3*randn(N));
    LRowR = p(p_idx)*MemR; %abs(1e3*ones(N)  + 0e3*randn(N));
    LColR = p(p_idx)*MemR; %abs(1e3*ones(N)  + 0e3*randn(N));
    
    s_orig = [5, zeros(1, N-1)]';       %Source Magnitude
    c  = zeros(2*N, 1);                 %Connection Matrix
    c(1:N-1) = 1;
%     c(N+1:end) = 1;
    R = zeros(N);
    LRe = zeros(N);
    for idx = N%1:N
        % idx = N;
        s = fShiftingMatrix(length(s_orig))^(idx-1) * s_orig;
        
        %% Simulation
        Circuit = fMacSpiceSim(N, s, c, MemR, LRowR, LColR);
        
        %% Plot Results
        close all
        
        %     figure;
        %     surf(Circuit.VO.value);
        %     title('VO Output');
        %     view([fig_az fig_el]);
        %     set(gca,'xdir','reverse')
        %
        %     figure;
        %     surf(Circuit.IO.value);
        %     title('IO Output');
        %     view([fig_az fig_el]);
        %     set(gca,'xdir','reverse')
        %
        %     figure;
        %     plot(1:N, Circuit.IO.value(N, :));
        %     title('IO Output');
        %     % view([fig_az fig_el]);
        %
        %     figure;
        %     plot(1:N-1, diff(Circuit.IO.value(N, :)));
        %     title('IO Output Gradient');
        
        %     close all
        io = Circuit.IO.value(N, :);
        
        %% Simple estimation of Line Resistances
        %     R = zeros(N);
        %     LRe = zeros(N);
        for i = 1:N
            for j = 1:N
                if i==idx
                    R(i, j) = s(i)/io(j);
                    LRe(i, j) = (R(i, j) - MemR(i, j))/(2*(N+1)-(i+j));
                end
            end
        end
    end
    Rp_estimate(p_idx) = mean(LRe(:));
end

Rp_error = Rp_estimate - p.*MemR(1);
Rp_err_percent = (abs(Rp_error)./Rp_estimate)*100;
% plot (p, Rp_error)