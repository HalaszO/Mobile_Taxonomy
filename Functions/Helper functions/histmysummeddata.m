function out=histmysummeddata(in,toplots)%% toplots={'Allsum','Pyrsum','Intsum','NGFsum','FSsum'};
%leng=0;
%savefigureasjpg=0;
for i=2:size(in.(char(toplots(1))),2)
    figure;
    temp=[];
    for j=1:length(toplots)
        toplot=in.(char(toplots(j)));
        subplot(length(toplots),1,j);
        alldatatohist=cell2mat(in.Allsum(2:end,i));
        datatohist=cell2mat(toplot(2:end,i));
        temp{j}=datatohist;
        if ~(min(alldatatohist)==max(alldatatohist))
        hist(datatohist,min(alldatatohist):(max(alldatatohist)-min(alldatatohist))/100:max(alldatatohist));
        else
            text(.5,.5,'minden eleme 0 :(')
        end
        title([char(toplots(j)),':  ', char(toplot(1,i)),'  min:',num2str(min(datatohist)),' max:',num2str(max(datatohist)),' mean:',num2str(mean(datatohist)),' mode:',num2str(mode(datatohist))]);

    end
            pause;
        close(gcf);
        
%     if savefigureasjpg==1
%         saveas(gcf, [char(toplot(1,i)),'.jpg']);
%         close(gcf);
%     else
%         [pmannw,hmannw] = ranksum(cell2mat(temp(2)), cell2mat(temp(3)), .05);
%         if hmannw==1
%             leng=leng+1;
%             out{leng,1}= char(toplot(1,i));
%             out{leng,2}= pmannw;
%             [num2str(pmannw),':  ',num2str(mean(cell2mat(temp(2)))),' vs  ',num2str(mean(cell2mat(temp(3))))]
% %             pause
%             close(gcf);
%         else
%             close(gcf);
%         end
%     end
end
end