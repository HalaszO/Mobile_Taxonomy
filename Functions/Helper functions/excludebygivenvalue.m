function data=excludebygivenvalue(data)
excludeday={'110602','110610','110615','110811','111026','120607'};
valuenames=data(1,:);
rscol=find(strcmp(valuenames,'RS'));
RSmax=4000;
rincol=find(strcmp(valuenames,'Rin'));
Rinmax=1000;
v0col=find(strcmp(valuenames,'vresting'));
v0max=0;
for i=size(data,1):-1:2
    if cell2mat(data(i,rscol))>RSmax
        data(i,:)=[];
    elseif cell2mat(data(i,rincol))>Rinmax
        data(i,:)=[];
    elseif cell2mat(data(i,v0col))>v0max
        data(i,:)=[];
    end
    todelete=0;
    for j=1:length(excludeday)
        if findstr(char(excludeday(j)),char(data(i,1)))
            todelete=1;
        end
    end
    if todelete==1
        data(i,:)=[];
    end

end
end