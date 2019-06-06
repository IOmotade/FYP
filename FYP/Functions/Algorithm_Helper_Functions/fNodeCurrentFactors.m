% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates current factors perceived from memristor's bottom node to
% Vdd/gnd roughly. (Assumes column terminating nodes are connected to gnd 
% or open)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% N (Integer) = Order of memristor array i.e. NxN array
% i, j (Integer) = Row, Column of interest
% node (Integer) = current node being calculated (i.e. node offset from
% column index of interest
% MemR (NxN Double) = Memristor array values
% LRowR (NxN Double) = Line Row Resistances array values
% LColR (NxN Double) = Line Column Resistances array values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% nCF (Nx2) = Upward & Downward current factors through current path from
% memristor bottom plate to Vdd/Vss roughly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nCF = fNodeCurrentFactors(nCF, N, i, j, node, MemR, LRowR, LColR)
%% Instantiate Variables
[downIdx, upIdx] = deal(1, 2);
if node==0
    nCF = NaN*zeros(N, 2);
end

% if node+i==2
%     x = 0;
% end

%% Recursion Termination Condition
if ((node+i==1) || (node+i==N))
    if (node~=0)
        if (node+i==1)
        [Rdown, Rup] = deal(Inf, MemR(node+i, j));
        elseif (node+i==N)
            [Rdown, Rup] = deal(LColR(node+i, j), MemR(node+i, j));
        end
        [nCF(node+i, downIdx), nCF(node+i, upIdx)] = fCurrentRatio(Rdown, Rup);
    else
        if i==1
            [nCF(node+i, downIdx), nCF(node+i, upIdx)] = deal(1, 0);
        elseif i==N
            [Rdown, Rup] = deal(LColR(node+i, j), ...
                fEquivalentResistance("Up", N, i, j, MemR, LRowR, LColR));
            [nCF(node+i, downIdx), nCF(node+i, upIdx)] = fCurrentRatio(Rdown, Rup);
        end
    end
    
    if node~=0
        return
    end
end

%% Recursive Section
% Going towards SMU
if (node>=0) && (node+1+i<=N)
    nCF = fNodeCurrentFactors(nCF, N, i, j, node+1, MemR, LRowR, LColR);
end
% Going away from SMU
if (node<=0) && (node-1+i>0)
    nCF = fNodeCurrentFactors(nCF, N, i, j, node-1, MemR, LRowR, LColR);
end

%% %Calculate Rdown and Rup for general nodes
if ~(((node+i==1) || (node+i==N)))
    % If traversing upwards in circuit
    if node<0
        [Rdown, Rup] = deal(fEquivalentResistance("Up", N, node+i, j, MemR, LRowR, LColR), ...
            MemR(node+i, j));
    % Else if traversing upwards in circuit
    elseif node>0
        [Rdown, Rup] = deal(fEquivalentResistance("Down", N, node+i, j, MemR, LRowR, LColR), ...
            MemR(node+i, j));
    else
        [Rdown, Rup] = deal(fEquivalentResistance("Down", N, node+i, j, MemR, LRowR, LColR), ...
            fEquivalentResistance("Up", N, node+i, j, MemR, LRowR, LColR));
    end
    %Calculate Current Factors
    [nCF(node+i, downIdx), nCF(node+i, upIdx)] = fCurrentRatio(Rdown, Rup);
end
end

