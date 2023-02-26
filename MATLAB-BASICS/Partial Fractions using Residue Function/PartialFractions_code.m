% %Exercise 1
 
% %a=[Denominator Coefficient]
% %b=[Numerator Coefficient]

b1=[2 5 3 6];
a1=[1 6 11 6];

[r1,p1,k1]=residue(b1,a1);
disp(r1);
disp(p1);
disp(k1);

%%%%%%%%%%%==================%%%%%%%%%%%

%Exercise 2

%a=[Denominator Coefficient]
%b=[Numerator Coefficient]

b2=[1 2 3];
a2=[1 3 3 1];

[r2,p2,k2]=residue(b2,a2);
disp(r2);
disp(p2);
disp(k2);
