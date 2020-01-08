clc; clear all;
%% 多段曲线拟合

%% 读取数据

filename = 'C:\Users\Joey\Desktop\成都\cl1\edge_segment1.txt';
[oX, oY, oZ] = textread(filename, '%f%f%f', 'delimiter', ',');

mX = mean(oX); mY = mean(oY);  mZ = mean(oZ);
sX = std(oX); sY = std(oY); sZ = std(oZ);

% 转换后的X,Y,Z坐标
X = (oX - mX); Y = (oY - mY); Z = (oZ - mZ);
total_point_number = size(X,1);
ocoor(:,1) = X; ocoor(:,2) = Y;
coor = sortrows(ocoor,2);

%% 分段
% 分段数
segment_number = 1;
% 距离累加计算
distance_accu = zeros(total_point_number,1);

for i = 2: total_point_number
    x1 = coor(i,1); y1 = coor(i,2);
    x2 = coor(i-1,1); y2 = coor(i-1,2);
    distance_accu(i) = distance_accu(i-1) + sqrt((x1-x2).^2+(y1-y2).^2);
end

distance_total = distance_accu(end);  % 距离总和
segment_length_threshold = 100.0;     % 分段长度阈值

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

%% 优化
greek = optimvar('greek', 3*segment_number); % 要迭代的变量，每一段三个
diffexpr = 0.0;
x_integ = @(t, mu, ka, ps) cos(mu + ka.*t + 0.5.*ps.*t.*t);
y_integ = @(t, mu, ka, ps) sin(mu + ka.*t + 0.5.*ps.*t.*t);
% 迭代初值

% 遍历每一段
for segment_iter = 1 : segment_number
   % 取出坐标
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
    
    for i = 1 : 2
        diffexpr = diffexpr + (fcn2optimexpr(diffun{segment_iter*2-1}, greek, u(i)) - coor_segment(i,1)).^2;
        diffexpr = diffexpr + (fcn2optimexpr(diffun{segment_iter*2}, greek, u(i)) - coor_segment(i,2)).^2;
    end
end

    x0.greek = zeros(1, 12);
    x0.greek(1) = 2.0;x0.greek(4) = 2.0;x0.greek(7) = 2.0;x0.greek(10) = 2.0;
    ssqprob = optimproblem('Objective', diffexpr);
    options = optimoptions('lsqlin');
    [sol, fval, exitflag, output] = solve(ssqprob, x0, 'Options', options); 
    sol.greek

