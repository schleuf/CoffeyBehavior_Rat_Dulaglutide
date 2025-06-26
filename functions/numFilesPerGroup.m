function [numtab] = numFilesPerGroup(mT)
    numtab = table(); 
    sessnum = sort(unique(mT.Session));
    numtab.Session = sessnum;
    C57F = table;
    C57M = table;
    CD1F = table;
    CD1M = table;
    for sn = 1:numel(sessnum)
        C57F = [C57F; table(numel(find(mT.Session == sessnum(sn) & mT.Strain == 'c57' & mT.Sex == 'Female')))];
        C57M = [C57M; table(numel(find(mT.Session == sessnum(sn) & mT.Strain == 'c57' & mT.Sex == 'Male')))];
        CD1F = [CD1F; table(numel(find(mT.Session == sessnum(sn) & mT.Strain == 'CD1' & mT.Sex == 'Female')))];
        CD1M = [CD1M; table(numel(find(mT.Session == sessnum(sn) & mT.Strain == 'CD1' & mT.Sex == 'Male')))];
    end
    C57F.Properties.VariableNames{'Var1'} = 'C57F';
    C57M.Properties.VariableNames{'Var1'} = 'C57M';
    CD1F.Properties.VariableNames{'Var1'} = 'CD1F';
    CD1M.Properties.VariableNames{'Var1'} = 'CD1M';
    numtab = [numtab, C57F, C57M, CD1F, CD1M];
end