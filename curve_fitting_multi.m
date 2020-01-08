clc; clear all;
%% ����������

%% ��ȡ����

filename = 'C:\Users\Joey\Desktop\�ɶ�\cl1\edge_segment1.txt';
[oX, oY, oZ] = textread(filename, '%f%f%f', 'delimiter', ',');

mX = mean(oX); mY = mean(oY);  mZ = mean(oZ);
sX = std(oX); sY = std(oY); sZ = std(oZ);

% ת�����X,Y,Z����
X = (oX - mX); Y = (oY - mY); Z = (oZ - mZ);
total_point_number = size(X,1);
ocoor(:,1) = X; ocoor(:,2) = Y;
coor = sortrows(ocoor,2);

%% �ֶ�
% �ֶ���
segment_number = 1;
% �����ۼӼ���
distance_accu = zeros(total_point_number,1);

for i = 2: total_point_number
    x1 = coor(i,1); y1 = coor(i,2);
    x2 = coor(i-1,1); y2 = coor(i-1,2);
    distance_accu(i) = distance_accu(i-1) + sqrt((x1-x2).^2+(y1-y2).^2);
end

distance_total = distance_accu(end);  % �����ܺ�
segment_length_threshold = 100.0;     % �ֶγ�����ֵ

low = 1; upp = 1;
for i = 2: total_point_number - 1
    low_ = floor(distance_accu(i - 1)/segment_length_threshold);
    hig_ = floor(distance_accu(i)/segment_length_threshold);
    
    if low_ ~= hig_
        high = i - 1;
        segment_label(segment_number, 1) = low;
        segment_label(segment_number, 2) = high;
        low = i;
        segment_number = segment_number + 1;
    end
    segment_label(segment_number, 1) = low;
    segment_label(segment_number, 2) = total_point_number;
end

%% �Ż�
l_bound = zeros(3*segment_number,1); h_bound = zeros(3*segment_number,1);
    for i = 1:segment_number
        l_bound(3*i-2) = 0.0;
        l_bound(3*i-1) = -0.05;
        l_bound(3*i-0) = -0.0025;
        h_bound(3*i-2) = 2*pi;
        h_bound(3*i-1) = 0.05;
        h_bound(3*i-0) = 0.0025;
    end
greek = optimvar('greek', 3*segment_number, 'LowerBound', l_bound, 'UpperBound', h_bound); % Ҫ�����ı�����ÿһ������

diffexpr = 0.0;
x_integ = @(t, mu, ka, ps) cos(mu + ka.*t + 0.5.*ps.*t.*t);
y_integ = @(t, mu, ka, ps) sin(mu + ka.*t + 0.5.*ps.*t.*t);
% ������ֵ

% ����ÿһ��
for segment_iter = 1 : 4
   % ȡ������
   coor_segment = coor(segment_label(segment_iter,1):segment_label(segment_iter,2),:);
   seg_point_number = size(coor_segment, 1);
   distance_accu_seg = zeros(seg_point_number, 1);
   x00 = coor_segment(1,1); y00 = coor_segment(1,2);

    for i = 2 : seg_point_number
        x1 = coor_segment(i,1); y1 = coor_segment(i,2);
        x2 = coor_segment(i-1,1); y2 = coor_segment(i-1,2);
        distance_accu_seg(i) = distance_accu_seg(i-1) + sqrt((x1-x2).^2+(y1-y2).^2);
    end
    distance_ratio = zeros(seg_point_number,1);
    curve_length_seg = distance_accu_seg(end);
    for i = 2 : seg_point_number - 1
        distance_ratio(i) = distance_accu_seg(i)/curve_length_seg;
    end

    distance_ratio(end) = 1;
    u = distance_ratio;
    
    diffun{segment_iter*2-1} = @(greek, uu) x00 + integral(@(t)x_integ(t, ...
    greek(segment_iter*3-2), greek(segment_iter*3-1), greek(segment_iter*3)), 0, uu*curve_length_seg);
    diffun{segment_iter*2} = @(greek, uu) y00 + integral(@(t)y_integ(t, ...
    greek(segment_iter*3-2), greek(segment_iter*3-1), greek(segment_iter*3)), 0, uu*curve_length_seg);
    
    for i = 1 : seg_point_number
        diffexpr = diffexpr + (fcn2optimexpr(diffun{segment_iter*2-1}, greek, u(i)) - coor_segment(i,1)).^2;
        diffexpr = diffexpr + (fcn2optimexpr(diffun{segment_iter*2}, greek, u(i)) - coor_segment(i,2)).^2;
    end
end

    x0.greek = zeros(1, 3*segment_number);
    for i = 1:segment_number
        x0.greek(3*i-2) = 1.5;
        x0.greek(3*i-1) = 0.0;
        x0.greek(3*i-0) = 0.0;
    end
    ssqprob = optimproblem('Objective', diffexpr);
    options = optimoptions('lsqlin');
    options.Display = 'iter';
    showproblem(ssqprob)
    [sol, fval, exitflag, output] = solve(ssqprob, x0, 'Options', options); 
    sol.greek
    
    %% curve generator
    s = linspace(0, curve_length_seg, 100);
    coor2 = zeros(100, 2);
    coor2(1,1) = x00; coor2(1,2) = y00;
    mu0 = sol.greek(10); kappa0 = sol.greek(11); psi = sol.greek(12);
    for i = 1 : 100
        % ����
        fun1 = @(t, mu0, kappa0, psi) cos(mu0 + t.*kappa0 + 0.5.*t.*t.*psi);
        fun2 = @(t, mu0, kappa0, psi) sin(mu0 + t.*kappa0 + 0.5.*t.*t.*psi);
        Q1 = integral(@(t)fun1(t,mu0,kappa0,psi),s(1),s(i));
        Q2 = integral(@(t)fun2(t,mu0,kappa0,psi),s(1),s(i));
        coor2(i,1) = Q1 + x00;
        coor2(i,2) = Q2 + y00;
    end
    plot(coor2(:,1),coor2(:,2),'b')
    hold on
    scatter(coor_segment(:,1),coor_segment(:,2), 'r')

