function totalmsglen = fDisplayInternalMessage(msg, totalmsglen)
if ~exist('totalmsglen', 'var')
    totalmsglen = 0;
end

fprintf(msg);
totalmsglen = totalmsglen + strlength(msg);

end