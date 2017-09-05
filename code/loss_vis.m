e=linspace(-1,1,51);
l1=abs(e);
l2=e.^2;

e_neg=e;
e_neg(e_neg>=0)=0;
e_neg(e_neg<0)=1;

rho1=e.*(0.25-e_neg);
rho2=e.*(0.5-e_neg);
rho3=e.*(0.75-e_neg);

figure;
plot(e,l1);
hold on;
plot(e,l2);
hold off;
legend('L1','L2','Location','southeast');
title('L_1 and L_2 loss');

figure;
plot(e,rho1);
hold on;
plot(e,rho2);
plot(e,rho3);
hold off;
legend('tau=0.25','tau=0.5','tau=0.75','Location','southeast');
title('Tilted loss');