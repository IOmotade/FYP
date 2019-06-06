% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates equivalent resistance(to top or bottom/gnd) perceived at memristor's
% bottom node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% type (String) = Either "Up" or "Down" to represent current flow path
% N (Integer) = Order of memristor array i.e. NxN array
% i, j (Integer) = Row, Column of interest
% MemR (NxN Double) = Memristor array values
% LRowR (NxN Double) = Line Row Resistances array values
% LColR (NxN Double) = Line Column Resistances array values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% Req (Double) = Equivalent Resistance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Req = fEquivalentResistance(type, N, i, j, MemR, LRowR, LColR)

if i==1
    if type=="Up"
        Req = Inf;%1e10000 or a really large number
        return
    end
end

if i==N
    if type=="Down"
        Req = LColR(i, j);
        return
    end
end

k = i - (type=="Up") + (type=="Down");

Rn = fEquivalentResistance(type, N, k, j, MemR, LRowR, LColR);
Req = LColR(i, j) + fParallelResistance(MemR(k, j), Rn);

end