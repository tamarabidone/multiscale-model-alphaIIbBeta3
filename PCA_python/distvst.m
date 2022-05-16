figure
plot(bent_d(:,1),bent_d(:,2),'color','k','linewidth',6);
hold on;
plot(int1_d(:,1),int1_d(:,2),'color','b','linewidth',6);
hold on;
plot(int2_d(:,1),int2_d(:,2),'color','g','linewidth',6);
hold on;
plot(open_d(:,1),open_d(:,2),'color','r','linewidth',6);
xlabel ('Time (ns)', 'fontsize',40);
ylabel ('Distance (nm)', 'fontsize',40);
set (gca, 'fontsize',44, 'linewidth', 8);

%xticklabels([0,1,2,3,4,5,6,7,8,9])
f=figure(1);
f.Position=[100 700 900 1000]
axis([0 500 5 22])
legend boxoff;
legend('Bent', 'Int1', 'Int2','Open','location','northeast');


