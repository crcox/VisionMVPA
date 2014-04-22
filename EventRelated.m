function eSeries = EventRelated(rois,nScans,inputTs)
% Take a raw timeseries, and organize it n seconds out per trial

% inputTs = temp_ts;

for k = 1  % # conditions
    for i = 1:length(rois)
        % Get number of scans in this condition
        for j = 1:nScans % Number of sessions / condition
            for trial = 1:34
                eOnset = 5.*(trial-1)+1;
                eSeries{k}{i}{j}(trial,:) = inputTs{k}{i}{j}(eOnset:eOnset+14); % - mean(inputTs{k}{i}{j}(eOnset:eOnset+3)); %14);
            end
%             for trial = 3:34 % Kicking out the first 2 trials, so we can
%             go back in time
%                 eOnset = 5.*(trial-1)+1;
%                 eSeries{k}{i}{j}(trial-2,:) = inputTs{k}{i}{j}(eOnset-10:eOnset+14); % - mean(inputTs{k}{i}{j}(eOnset:eOnset+3)); %14);
%             end
        end
    end
end