function [valuenames,fnames]=getvaluenames(datasum,needfnames)
cells=fieldnames(datasum);
features=fieldnames(datasum.(char(cells(1))));
    j=0;
    for i=1:length(features)
        featlength = zeros(length(cells),1);
        for c=1:length(cells)
            featlength(c) = length(datasum.(char(cells(c))).(char(features(i))));
        end
        if all(featlength == 1)
            j=j+1;
            valuenames(j)=features(i);
        end
    end
    if exist('needfnames','var') && needfnames>0
        for i=1:length(cells)
            temp=char(cells(i));
            fnames{i}=[temp(6:end),'.mat'];
        end
    end
end

