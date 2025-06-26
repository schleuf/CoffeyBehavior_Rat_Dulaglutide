function [mT] = removeExcludedData(mT, mKey)
    % SS added
    RemoveSession = zeros([length(mT.Chamber),1]);
    sub_RemSess = mKey.RemoveSession;
    for sub = 1:length(sub_RemSess)
        sess = sub_RemSess{sub}(2:end-1);
        if ~isempty(sess)
            sess = strsplit(sess, ' ');
            tag = mKey.TagNumber(sub);
            for s = 1:length(sess)
                ind = find((mT.TagNumber == tag) .* (mT.Session == str2double(sess{s})));
                RemoveSession(ind) = 1;
            end
        end      
    end
    mT(find(RemoveSession),:) = [];
end
