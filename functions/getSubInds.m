function subind = getSubInds(tab, subset)
    subind = ones([height(tab), 1]);
    for sb = 1:length(subset)
        tmp = zeros([height(tab), 1]);
        for m = 1:length(subset{sb}{2})
            tmp = tmp | (tab.(subset{sb}{1}) == subset{sb}{2}{m});
        end
        subind = subind .* tmp;
    end
    subind = logical(subind); 
end