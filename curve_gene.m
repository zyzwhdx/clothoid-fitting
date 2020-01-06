clc;clear all;
%% clothoid generator
% 单位是m
% x0和y0是第一个点对应的横纵坐标
s_start = 0; s_end = 100; s_num = 51; s_interval = s_end - s_start;
s = linspace(s_start, s_end, s_num);
x00 = 750; y00= 1025.0;

% 曲线参数设定
mu0 = pi;
kappa0 = 0.001;
psi = 0.00001;

% 点坐标为行向量
coor = zeros(s_num, 2);
coor(1,1) = x00; coor(1,2) = y00;

for i = 2 : s_num
    % 积分
    fun1 = @(t, mu0, kappa0, psi) cos(mu0 + t.*kappa0 + 0.5.*t.*t.*psi);
    fun2 = @(t, mu0, kappa0, psi) sin(mu0 + t.*kappa0 + 0.5.*t.*t.*psi);
    Q1 = integral(@(t)fun1(t,mu0,kappa0,psi),s(1),s(i));
    Q2 = integral(@(t)fun2(t,mu0,kappa0,psi),s(1),s(i)) + (rand-0.5).*0.3;
    coor(i,1) = Q1 + x00;
    coor(i,2) = Q2 + y00;
end

% 旋转坐标系
% theta = 30.0.*pi./180.0;
% rotation = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% coor_ = coor*rotation;
scatter(coor(:,1),coor(:,2))

% fp = fopen('clothoid1.txt','wt');
% fprintf(fp, '%f,%f\n', coor_);
% fclose(fp);

%% optimizer

% func_int = @(t) (t-1);
% func = @(u) integral(@(t)func_int(t),0,u);
% 
% u0 = 1.3;
% options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');
% options.Display = 'iter';
% [xx, fval, exitflag, output] = fminunc(func, x0, options);

% 迭代初值
u = linspace(0, 1, s_num);
m0 = 3.0; k0 = 0.000; p0 = 0.0000;

for major_iter = 1:5
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
        fc = distance_point_clothoid(coor(i, 1), coor(i, 2), c*s_interval+s_start, m0, k0, p0, x00, y00);
        fd = distance_point_clothoid(coor(i, 1), coor(i, 2), d*s_interval+s_start, m0, k0, p0, x00, y00);
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
diffun1 = @(greek, uu) x00 + integral(@(t)x_integ(t, greek(1), greek(2), greek(3)), 0, uu*s_interval+s_start);
diffun2 = @(greek, uu) y00 + integral(@(t)y_integ(t, greek(1), greek(2), greek(3)), 0, uu*s_interval+s_start);

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
options
sol.greek

m0 = sol.greek(1); k0 = sol.greek(2); p0 = sol.greek(3);

end

%% curve generator
coor2 = zeros(s_num, 2);
coor2(1,1) = x00; coor2(1,2) = y00;
mu0 = sol.greek(1); kappa0 = sol.greek(2); psi = sol.greek(3);
for i = 2 : s_num
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
