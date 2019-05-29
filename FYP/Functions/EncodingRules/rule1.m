function posMemRValues = rule1(numBitspRes, MemRLowLim, MemRUpLim)

if ~(exist('MemRLowLim', 'var') && exist('MemRUpLim', 'var'))
%     [MemRLowLim, MemRUpLim] = deal(100, 10e3);
    [MemRLowLim, MemRUpLim] = deal(10, 10e6);
end
posMemRValues = MemRLowLim + ((MemRUpLim - MemRLowLim)/(2^numBitspRes-1))*((1:2^numBitspRes)-1);

end