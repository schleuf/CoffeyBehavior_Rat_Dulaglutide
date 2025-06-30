function [mT] = interpoweight(mT, plotWeight)
    tags = unique(mT.TagNumber);
    for t = 1:length(tags)
        interp_sessions = char(unique(mT.SessionWeightUpdated(mT.TagNumber == tags(t))));
        interp_sessions = interp_sessions(2:end-1);
        interp_sessions = cell2mat(cellfun(@(s) sscanf(s, '%d,'), cellstr(interp_sessions), 'UniformOutput', false));

        
        % disp(tags(t))
        weights = mT.Weight(mT.TagNumber == tags(t));
        sessions = mT.Session(mT.TagNumber == tags(t));
        
        if plotWeight
            figure
            hold on
            plot(sessions, weights, 'Marker', '.')
            title([tags(t), ' weights'])
            hold off
        end

        temp = weights;
        for s = 1:length(interp_sessions)-1
            first = weights(sessions == interp_sessions(s));
            last = weights(sessions == interp_sessions(s+1));
            if ~isempty(last) & ~isempty(first)
                num = interp_sessions(s+1) - interp_sessions(s) + 1;
                try
                    interp = linspace(first, last, num);
                catch
                    disp('why')
                end
                for n = 1:num
                    temp(sessions == interp_sessions(s) + n - 1) = interp(n);
                end
            end
        end
        mT.Weight(mT.TagNumber == tags(t)) = temp;
        
        if plotWeight
            hold on
            plot(sessions, temp)
            hold off
        end
    end
end