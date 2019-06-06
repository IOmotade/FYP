function totalmsglen = fDisplayInternalMessage(msg, totalmsglen)
if ~exist('totalmsglen', 'var')
    totalmsglen = 0;
end

fprintf(strcat(msg, '\n'));
totalmsglen = totalmsglen + numel(msg) + 1;
% totalmsglen = totalmsglen + strlength(msg);

end