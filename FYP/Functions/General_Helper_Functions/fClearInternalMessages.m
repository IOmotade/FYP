function totalmsglen = fClearInternalMessages(totalmsglen)
if ~exist('totalmsglen', 'var')
    totalmsglen = 0;
end

fprintf(strcat(repmat('\b', 1, totalmsglen)));
totalmsglen = totalmsglen - totalmsglen;

end