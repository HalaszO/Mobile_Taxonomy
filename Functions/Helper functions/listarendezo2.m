function listnew=listarendezo2(listold)
header=(listold(1,:));
listnew.pall=listold;
listold(1,:)=[];
p=1;
listnew.(['p',num2str(p)])=[header;listold(1,:)];
listold(1,:)=[];
while size(listold,1)>0
    tmp=char(listnew.(['p',num2str(p)])(end,1));
    tmp=tmp(1:6);
    aas=[];
    for i=1:size(listold,1)
        aa=char(listold(i,1));
        aas{i}=aa(1:6);
    end
    cmpmask=strcmp(aas,tmp);
    cmpmask(cmpmask>1)=0;
    if sum(cmpmask)>0
        lol=sum(cmpmask);
        for i=1:lol
            listnew.(['p',num2str(p)])(size(listnew.(['p',num2str(p)]),1)+1, :)=listold(find(cmpmask==1,1,'last'),:);
            listold(find(cmpmask==1,1,'last'),:)=[];
            cmpmask(find(cmpmask==1,1,'last'))=0;
        end
       p=p+1;
    end
    if size(listold,1)>0
    %listnew.(['p',num2str(p)])(size(listnew,1)+1,:)=listold(1,:);
    listnew.(['p',num2str(p)])=[header;listold(1,:)];
    listold(1,:)=[];
    end
end
end