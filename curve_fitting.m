clc; clear all;
%% 跑真实数据

%% 读取数据
filename = 'C:\Users\Joey\Desktop\成都\cl1\edge_segment1.txt';
% 初始的XYZ坐标
[oX, oY, oZ] = textread(filename, '%f%f%f', 'delimiter', ',');

% 均值和标准差计算
mX = mean(oX); mY = mean(oY);  mZ = mean(oZ);
sX = std(oX); sY = std(oY); sZ = std(oZ);

% 转换后的X,Y,Z坐标,XY转换一下
X = (oX - mX);
Y = (oY - mY);
Z = (oZ - mZ);

coor(:,1) = X;
coor(:,2) = Y;

% 按照Y方向排序
coor = sortrows(coor, 2);
% scatter(X,Y);

x00 = coor(1, 1);
y00 = coor(1, 2);
s_num = size(coor, 1);


%% 生成初始的u
distance_accu = zeros(s_num,1);
distance_ratio = zeros(s_num,1);
for i = 2: s_num
    x1 = coor(i,1); y1 = coor(i,2);
    x2 = coor(i-1,1); y2 = coor(i-1,2);
    distance_accu(i) = distance_accu(i-1) + sqrt((x1-x2).^2+(y1-y2).^2);
end

distance_total = distance_accu(end);

for i = 2: s_num - 1
    distance_ratio(i) = distance_accu(i)/distance_total;
end

distance_ratio(end) = 1;
u = distance_ratio;

%% 优化器
m0 = 1.95; k0 = 0.000; p0 = 0.0000;
curve_length = distance_total;
major_iter_num = 3;

for major_iter = 1 : major_iter_num
    % golden search
    for i = 2 : s_num-1
        R = 0.618034;
        eps = 1e-5;
        a = u(i-1);
        b = u(i+1);
        c = u(i+1) - (u(i+1)-u(i-1))*R;
        d = u(i-1) + (u(i+1)-u(i-1))*R;
        n_iter = ceil(-2.0780869*log(eps/abs(b-a)));

        for j =1:n_iter
            fc = distance_point_clothoid(coor(i, 1), coor(i, 2), c*curve_length, m0, k0, p0, x00, y00);
            fd = distance_point_clothoid(coor(i, 1), coor(i, 2), d*curve_length, m0, k0, p0, x00, y00);
            if fc<fd
                b = d;
            else
                a = c;
                c = b-(b-a)*R;
                d = a+(b-a)*R;
            end
        end
        u(i) = 0.5*(a+b);
    end

    % 待优化的参数
    greek = optimvar('greek',3);

    x_integ = @(t, mu, ka, ps) cos(mu + ka.*t + 0.5.*ps.*t.*t);
    y_integ = @(t, mu, ka, ps) sin(mu + ka.*t + 0.5.*ps.*t.*t);
    diffun1 = @(greek, uu) x00 + integral(@(t)x_integ(t, greek(1), greek(2), greek(3)), 0, uu*curve_length);
    diffun2 = @(greek, uu) y00 + integral(@(t)y_integ(t, greek(1), greek(2), greek(3)), 0, uu*curve_length);

    diffexpr = 0.0;
    for i = 1 : s_num
        diffexpr = diffexpr + (fcn2optimexpr(diffun1, greek, u(i)) - coor(i,1)).^2;
        diffexpr = diffexpr + (fcn2optimexpr(diffun2, greek, u(i)) - coor(i,2)).^2;
    end

    % 迭代初值，曲线参数
    x0.greek = [m0, k0, p0];
    ssqprob = optimproblem('Objective', diffexpr);
    options.Display = 'iter';
    options = optimoptions('lsqlin');
    [sol, fval, exitflag, output] = solve(ssqprob, x0, 'Options', options); 
    sol.greek

    m0 = sol.greek(1); k0 = sol.greek(2); p0 = sol.greek(3);
    
    % 更新曲线长度
    search(:,1) = linspace(curve_length - 20, curve_length + 20, 4001);
    for i = 1:4001
        search(i,2) = abs(integral(@(t)x_integ(t, m0, k0, p0), 0, search(i,1)) + x00 - coor(end,1));
    end
    search = sortrows(search, 2);
    
    if abs(tmp_length - curve_length) <= 20
        curve_length = tmp_length;
    end
end

%% curve generator
s = linspace(0, curve_length, 100);
coor2 = zeros(100, 2);
coor2(1,1) = x00; coor2(1,2) = y00;
mu0 = sol.greek(1); kappa0 = sol.greek(2); psi = sol.greek(3);
for i = 1 : 100
    % 积分
    fun1 = @(t, mu0, kappa0, psi) cos(mu0 + t.*kappa0 + 0.5.*t.*t.*psi);
    fun2 = @(t, mu0, kappa0, psi) sin(mu0 + t.*kappa0 + 0.5.*t.*t.*psi);
    Q1 = integral(@(t)fun1(t,mu0,kappa0,psi),s(1),s(i));
    Q2 = integral(@(t)fun2(t,mu0,kappa0,psi),s(1),s(i));
    coor2(i,1) = Q1 + x00;
    coor2(i,2) = Q2 + y00;
end
plot(coor2(:,1),coor2(:,2),'b')
hold on
scatter(coor(:,1),coor(:,2), 'r')

%% 计算点到回旋曲线距离
function [distance] = distance_point_clothoid(x, y, u, mu, ka, ps, x0, y0)
    fun1 = @(t, mu, ka, ps) cos(mu + t.*ka + 0.5.*t.*t.*ps);
    fun2 = @(t, mu, ka, ps) sin(mu + t.*ka + 0.5.*t.*t.*ps);
    Q1 = integral(@(t)fun1(t,mu,ka,ps), 0, u);
    Q2 = integral(@(t)fun2(t,mu,ka,ps), 0, u);
    
    x_dst = Q1 + x0;
    y_dst = Q2 + y0;
    
    distance = sqrt((x-x_dst).^2 + (y-y_dst).^2);
end