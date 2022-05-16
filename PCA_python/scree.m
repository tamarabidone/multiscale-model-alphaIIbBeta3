int1_v=[1085462.56302558,  760968.7782987 ,   47882.37572335, 29241.80487156, 5120.77291959,    3830.91293447,3254.54454225,    1878.81261703,    1555.00965831, 1287.02461373];
bent_v= [413836.552431  ,  70333.7759093 ,  47032.94074006,  18947.70360797, 9350.09585828,   7394.41829451,   3776.78403722,   1384.05242144,1309.32989062,    768.40754044];
int2_v=[1.12712657e+07, 1.98096045e+05, 1.37649145e+05, 5.60289918e+04, 2.32393936e+04, 1.02670920e+04, 7.46197678e+03, 6.13564357e+03,4.38488689e+03, 3.78997055e+03];
open_v=[208338.26014155, 109277.51071041,  79058.66993882,  55579.48550782,31833.26171773,  20064.18892489,  13981.62015031,   7167.05861626,6259.24372981,   4665.03651701];
plot(mode,open_v,'-ko','MarkerFaceColor','r','MarkerSize',12,'color','r','linewidth',5);
xlabel ('Principal Components', 'fontsize',40);
ylabel ('Variance', 'fontsize',40);
set (gca, 'fontsize',44, 'linewidth', 4);
f=figure(1);
f.Position=[100 300 1000 600]
axis([0 10 0 210000 ])