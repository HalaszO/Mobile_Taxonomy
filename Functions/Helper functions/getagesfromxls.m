% function datasum=getagesfromxls(datasum)
 [fname,path] = uigetfile({'*.xls'});
[num, txt, raw] = xlsread([path,fname]);
agecol=1;
exp=fieldnames(datasum);
while ~strcmp('Age',raw(1,agecol))
    agecol=agecol+1;
end
for i=2:size(raw,1)
    experiment=cell2mat(raw(i,1));
    if ~isnan(experiment)
        for j=1:length(exp)
            expact=char(exp(j));
            if strcmp(experiment, expact(6:5+length(experiment)))
                temp=char(raw(i,agecol));
                datasum.(expact).age=str2num(temp(2:end));
            end
        end
    end
end
clear txt num raw fname path exp experiment temp
% end