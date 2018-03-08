load('NSSvsgazeoffset');
f = nssgaze;
list=jet(20);

for i=1:20
    if f{i,2} < 3
        valid(i) = 1;
    else
        valid(i) = 0;
    end
end
figure;
subplot(1,2,1);

for i=1:20
    if ~valid(i)
        scatter(f{i,3}, f{i,2},160,[0.8 0 0],'filled','MarkerEdgeColor',[0 0 0]);
        hold on
    else
        scatter(f{i,3}, f{i,2},160,[0.5 0.5 0.5],'filled','MarkerEdgeColor',[0 0 0]);
    end
end
set(gca,'FontSize',18);
xlabel('NSS score','FontSize',24);
ylabel('Average gaze offset)','FontSize',24);
axis square
nss = cell2mat(f(:,3));
offset = cell2mat(f(:,2));

[r , p] = corr(offset, nss, 'type', 'Spearman');
r =

   -0.6677


p =

    0.0017


subplot(1,2,2);


j = 1;
for i=1:20
    
    if valid(i)
        scatter(f{i,3}, f{i,4},160,[0.5 0.5 0.5],'filled','MarkerEdgeColor',[0 0 0]);
        hold on
        values(j) = f{i,4};
        j = j+1;
    end        
    
    nss(i) = f{i,3};
%     else
%         scatter(axes1, f{i,3}, f{i,4},160,[1 1 1],'filled','MarkerEdgeColor',[0 0 0]);
%         hold on
%     end
    
end

plot(0:7,0:7,'--k','LineWidth',1)
xlabel('NSS score','FontSize',24);
ylabel('NSS score (unrelated)','FontSize',24);
xlim([0 7]);
ylim([0 7]);
set(gca,'xtick',0:2:8, 'ytick',0:2:8);
axis square


legend(le);

