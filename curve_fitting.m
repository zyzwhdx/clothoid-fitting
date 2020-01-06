clear; clc;
filename = 'C:\Users\Joey\Desktop\成都\cl1\edge.txt';
% 初始的XYZ坐标
[oX, oY, oZ] = textread(filename, '%f%f%f', 'delimiter', ',');

% 均值和标准差计算
mX = mean(oX); mY = mean(oY);  mZ = mean(oZ);
sX = std(oX); sY = std(oY); sZ = std(oZ);

% 转换后的X,Y,Z坐标,XY转换一下
X = (oY - mY)/sY;
Y = (oX - mX);
Z = (oZ - mZ);

sp = spaps(X, Y, 0.1);
fnplt(sp)