function posMemRValues = rule2(numBitspRes, MemRLowLim, MemRUpLim)

if ~(exist('MemRLowLim', 'var') && exist('MemRUpLim', 'var'))
    [MemRLowLim, MemRUpLim] = deal(10, 10e6);
end
[lgMemRLowLim, lgMemRUpLim] = deal(log10(MemRLowLim), log10(MemRUpLim));
posMemRValues = lgMemRLowLim + ((lgMemRUpLim - lgMemRLowLim)/(2^numBitspRes-1))*((1:2^numBitspRes)-1);
posMemRValues = 10.^posMemRValues;

end