function totalmsglen = fDisplayInternalMessage(msg, totalmsglen)
if ~exist('totalmsglen', 'var')
    totalmsglen = 0;
end

fprintf(strcat(msg, '\n'));

if isunix
    extra = 1;
elseif ispc
    extra = 3;
end

totalmsglen = totalmsglen + numel(msg) + extra;
% totalmsglen = totalmsglen + strlength(msg);

end