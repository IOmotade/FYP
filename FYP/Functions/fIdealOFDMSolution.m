function [corr_fmag] = fIdealOFDMSolution(N, vs_mag, MemR, LRowR, LColR)
%% Generate Accurate Frequency Components
Vs = [vs_mag; zeros(N-1, 1)];
[corr_fmag, ~] = deal(zeros(N));
J = fShiftingMatrix(N);
for i = 1:N
    for j = 1:N
        Circuit = fMacSpiceSim(N, (J^(i-1))*Vs, MemR, LRowR, LColR);
        im = Circuit.IO.value(end, :);
        corr_fmag(i, j) = im(j);
        %             corr_fval(i, j) = fsource(i);
    end
    %         corr_Vbias(i, :) = Circuit.VI.value(i, :);
end
% fmag = corr_fmag;

end