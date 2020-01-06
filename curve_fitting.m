clear; clc;
filename = 'C:\Users\Joey\Desktop\�ɶ�\cl1\edge.txt';
% ��ʼ��XYZ����
[oX, oY, oZ] = textread(filename, '%f%f%f', 'delimiter', ',');

% ��ֵ�ͱ�׼�����
mX = mean(oX); mY = mean(oY);  mZ = mean(oZ);
sX = std(oX); sY = std(oY); sZ = std(oZ);

% ת�����X,Y,Z����,XYת��һ��
X = (oY - mY)/sY;
Y = (oX - mX);
Z = (oZ - mZ);

sp = spaps(X, Y, 0.1);
fnplt(sp)