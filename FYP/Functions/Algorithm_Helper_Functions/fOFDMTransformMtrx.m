
function H =fOFDMTransformMtrx(N)

R = log2(N);

H = zeros(R+1, N);
for idx = 0:R
    h = [];
    scale = 1;
    while length(h)~=N
        h = [h, scale*ones(1, N/(2^(idx)))];
        scale = scale * -1;
    end
    H(idx+1, :) = h;
end


end
