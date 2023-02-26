%%Exercise1%%
num = 20;
den = [3 8 20];
H = tf(num,den);
display(H);

stepplot(H);


%%Exercise2%%

num1 = [1 8 15];
den1 = [1 9 14 0];
[numZ,denP] = tf2zp(num1,den1);
display(numZ);
display(denP);
