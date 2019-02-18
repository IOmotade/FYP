function hft = haar_fft(x)

if(mod(length(x), 2)==1)
    x = [x, 0];
end

hft = zeros(size(x));

hft(1) = mean(x, 2);
t = 0:2*pi/(length(x)-1):2*pi;
for idx = 2:length(x)
    h = square(idx*t);
    hft(2) = sum(x./h);
end

end