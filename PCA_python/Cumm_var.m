bent_var= [0.71193387, 0.83293091, 0.91384291, 0.94643914, 0.96252436, 0.97524517, 0.98174247, 0.98412349, 0.98637597, 0.98769788]
int1_var=[0.55650635, 0.94664782, 0.97119666, 0.98618865, 0.98881402, 0.9907781, 0.99244667, 0.99340992, 0.99420716, 0.994867]
int2_var=[0.95996549, 0.97683718, 0.98856066, 0.99333261, 0.99531189, 0.99618633,0.99682186, 0.99734443, 0.99771789, 0.99804068]
open_var=[0.37316537, 0.56889794, 0.71050399, 0.81005526, 0.86707345, 0.90301145, 0.92805465, 0.94089194, 0.95210319, 0.96045898]
modes=[0,1,2,3,4,5,6,7,8,9]



figure(1)
set (gcf, 'color', 'w');
plot(modes, bent_var,'k',  'linewidth',6);
hold on
plot(modes, int1_var, 'b', 'linewidth',6);
hold on
plot(modes, int2_var,'g',  'linewidth',6);
hold on
plot(modes, open_var,'r', 'linewidth',6);
xlabel ('Principal Component', 'fontsize',40);
xticklabels([0,1,2,3,4,5,6,7,8,9])
f=figure(1);
f.Position=[100 300 800 600]
axis([0 9 0.35 1.05])
ylabel ('Cumulative Variance', 'fontsize',40);
set (gca, 'fontsize',44, 'linewidth', 8)
legend boxoff;
legend('Bent', 'Int1', 'Int2','Open','location','southeast');


%print('-dtiff','-r500',gcf)