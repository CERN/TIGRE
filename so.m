clear;clc;
numr11=[1 2 3];
denr11=[1 2.273 14.03 20.15 27.82];

R11=tf(numr11,denr11);

numr12=[3 2 1];

R12=tf(numr12,denr11);

numr21=[4 5 6];

R21=tf(numr21,denr11);

numr22=[7 8 9];

R22=tf(numr22,denr11);

R=[R11,R12;R21,R22]

% sys1= ss(R)
% 
% sys_final=minreal(sys1) 