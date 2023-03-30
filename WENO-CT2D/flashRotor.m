% flashRotor.m

TIME = 2000;

% Q1 = str2num(fileread('Q1.txt'));
% Q2 = str2num(fileread('Q2.txt'));
% Q3 = str2num(fileread('Q3.txt'));
% Q4 = str2num(fileread('Q4.txt'));
% Q5 = str2num(fileread('Q5.txt'));
% Q6 = str2num(fileread('Q6.txt'));
Q1 = str2num(fileread('Q1OT.txt'));
Q2 = str2num(fileread('Q2OT.txt'));
Q3 = str2num(fileread('Q3OT.txt'));
Q4 = str2num(fileread('Q4OT.txt'));
Q5 = str2num(fileread('Q5OT.txt'));
Q6 = str2num(fileread('Q6OT.txt'));
Q1 = Q1(:,1);
Q2 = Q2(:,1);
Q3 = Q3(:,1);
Q4 = Q4(:,1);
Q5 = Q5(:,1);
Q6 = Q6(:,1);
% T = fileread('T.txt');
% T = str2num(T);

STOP = length(T);

N = round(sqrt(length(Q1)/STOP));

Q1 = reshape(Q1,N,N,STOP);
Q2 = reshape(Q2,N,N,STOP);
Q3 = reshape(Q3,N,N,STOP);
Q4 = reshape(Q4,N,N,STOP);
Q5 = reshape(Q5,N,N,STOP);
Q6 = reshape(Q6,N,N,STOP);

QBP = (Q5.^2 + Q6.^2)/2;
QP = (Q4 - 0.5*(Q2.^2 + Q3.^2)./Q1 - 0.5*(Q5.^2 + Q6.^2))*(2/3);
QMach = ((Q2.^2 + Q3.^2)./Q1.^2).^0.5./((5/3*QP./Q1).^0.5);


