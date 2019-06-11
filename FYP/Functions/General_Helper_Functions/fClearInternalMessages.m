function totalmsglen = fClearInternalMessages(totalmsglen)
if ~exist('totalmsglen', 'var')
    totalmsglen = 0;
end

% log_file = true;
%% Clear sim.log file
% filename = 'sim.log';
% tmpfilename = ['tmp_', filename];
% if totalmsglen~=0 && log_file
%     try
%         [status, cmdout] = system(sprintf('wc -l %s', filename));
%
%         if(status==0)
%             scanCell = textscan(cmdout,'%u %s');
%             lineCount = double(scanCell{1});
%         else
%             fprintf(1,'Failed to find line count of %s\n', filename);
%             lineCount = -1;
%         end
%
%         extralines = 0;
%         [status, ~] = system(sprintf('head -n%d %s > %s',...
%             lineCount+extralines, filename, tmpfilename));
%
%         if status==0
%             movefile(tmpfilename, filename)
%         end
%
%     catch
%     end
% else
fprintf(strcat(repmat('\b', 1, totalmsglen)));
% end

totalmsglen = totalmsglen - totalmsglen;

end