
K= [252574.5022, 18614.45473,12927.53255,12239.44079];
str=categorical({'Bent','Int1','Int2','Open'});

color= ['k','b','g','r'];
figure,
h=bar(str,K,'FaceColor','flat');
h.CData(1,:)=[0 0 0]
h.CData(2,:)=[0 0 1]
h.CData(3,:)=[0 1 0]
h.CData(4,:)=[1 0 0]



%set(b(2),'FaceColor','b');
%b(3).FaceColor='g';
%b(4).FaceColor='r';
xlabel ('Structure', 'fontsize',40);
ylabel ('k (pN/\mum)', 'fontsize',40);
%ylim([0 130])
f=figure(1);
f.Position=[100 600 800 1000]
set (gca, 'fontsize',44, 'linewidth', 8)
%legend boxoff;
%legend('Bent', 'Int1', 'Int2','Open','location','southeast');