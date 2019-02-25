%%
% load('Circuitfsamp10kHz_nopara.mat')
%Voltage
V = zeros([N N]);
I = zeros([N N]);
val = zeros(size(t));
for idx = 1:N
    for tIdx=1:nsamp
        %     val(tIdx) = Circuit(tIdx).VO.value(end, 1);
        %         val(tIdx) = Circuit(tIdx).VI.value(end, 1);
        %         val(tIdx) = Circuit(tIdx).IO.value(end, 1);
        %     val(tIdx) = Circuit(tIdx).II.value(end, 1);
        val(tIdx) = Circuit(tIdx).VS.value(idx, 1);
        %     val(tIdx) = Circuit(tIdx).IS.value(end, 1);
        %     val(tIdx) = Circuit(tIdx).X.value(end, 1);
    end
    
    hft = haar_fft(val);
    norm_freq = (((1:length(hft))-1)/length(hft));%*fsamp;
    % norm_freq(abs(hft)==max(abs(hft)))*fsamp
    hft(abs(hft)<max(abs(hft))/1.5) = 0; %filter
    hift = haar_ifft(hft);
    
    %Extract frequency components
    val_freqcomp = zeros(1, N);
    f_read = zeros(1, N);
    for i = 1:N
        minIdx = abs(norm_freq-(fsource(i)/fsamp));
        val_freqcomp(i) = hft(minIdx == min(minIdx, [], 2));
        f_read(i) = norm_freq(minIdx == min(minIdx, [], 2));
    end
    val_freqcomp;
    V(idx, :) = val_freqcomp;
end

%Current
for idx = 1:N
    for tIdx=1:nsamp
        %     val(tIdx) = Circuit(tIdx).VO.value(end, 1);
        %         val(tIdx) = Circuit(tIdx).VI.value(end, 1);
                val(tIdx) = Circuit(tIdx).IO.value(end, idx);
        %     val(tIdx) = Circuit(tIdx).II.value(end, 1);
%         val(tIdx) = Circuit(tIdx).VS.value(idx, 1);
        %     val(tIdx) = Circuit(tIdx).IS.value(end, 1);
        %     val(tIdx) = Circuit(tIdx).X.value(end, 1);
    end
    
    hft = haar_fft(val);
    norm_freq = (((1:length(hft))-1)/length(hft));%*fsamp;
    % norm_freq(abs(hft)==max(abs(hft)))*fsamp
%     hft(abs(hft)<max(abs(hft))/1.5) = 0; %filter
    hift = haar_ifft(hft);
    
    %Extract frequency components
    val_freqcomp = zeros(1, N);
    f_read = zeros(1, N);
    for i = 1:N
        minIdx = abs(norm_freq-(fsource(i)/fsamp));
        val_freqcomp(i) = hft(minIdx == min(minIdx, [], 2));
        f_read(i) = norm_freq(minIdx == min(minIdx, [], 2));
    end
    val_freqcomp;
    I(idx, :) = val_freqcomp;
end

% Solve Memristor array
% I = transpose(G)*vs
% I = G' * V;
G = transpose(I*inv(V));
R_est = 1./G;
M = fUnits(MemR, 'Ohm')
fUnits(R_est, 'Ohm')
% Error = fUnits((MemR - R_est), 'Ohm')