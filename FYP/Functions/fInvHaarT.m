function h = fInvHaarT(hft)

% if(mod(length(hft), 2)==1)
%     hft = [hft, 0];
% end

h = zeros(size(hft));

t = 0:2*pi/(length(hft)-1):2*pi;
for idx = 1:length(hft)
    h = h + hft(idx)*square((idx-1)*t)/2;
%     plot(square((idx-1)*t));
%     title("Current Figure")
end
% h = h/2;
end