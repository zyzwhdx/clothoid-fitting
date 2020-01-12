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
        high = i;
        segment_label(segment_number, 1) = low;
        segment_label(segment_number, 2) = high;
        low = i;
        segment_number = segment_number + 1;
    end
    segment_label(segment_number, 1) = low;
    segment_label(segment_number, 2) = total_point_number;
end

%% 优化
greek_result = zeros(segment_number, 3);
seg_length_result = zeros(segment_number, 1);
coor_begin = zeros(segment_number, 2);

for major_iteration = 1 : 2
    % 迭代变量边界
    l_bound = zeros(3*segment_number,1); h_bound = zeros(3*segment_number,1);
    for i = 1:segment_number
        l_bound(3*i-2) = 0.0; l_bound(3*i-1) = -0.05; l_bound(3*i-0) = -0.0025;
        h_bound(3*i-2) = 2*pi; h_bound(3*i-1) = 0.05; h_bound(3*i-0) = 0.0025;
    end
    greek = optimvar('greek', 3*segment_number, 'LowerBound', l_bound, 'UpperBound', h_bound); % 要迭代的变量，每一段三个
    % 变量初值
    x0.greek = zeros(1, 3*segment_number);
    if major_iteration < 2
        for i = 1:segment_number
            x0.greek(3*i-2) = 1.5; x0.greek(3*i-1) = 0.0; x0.greek(3*i-0) = 0.0;
        end
    else
        x0.greek = sol.greek;
    end

    % 积分函数
    diffexpr = 0.0;
    x_integ = @(t, mu, ka, ps) cos(mu + ka.*t + 0.5.*ps.*t.*t);
    y_integ = @(t, mu, ka, ps) sin(mu + ka.*t + 0.5.*ps.*t.*t);

    % 遍历每一段，segment_iter
    for segment_iter = 1 : 2
       % 取出坐标
       coor_segment = coor(segment_label(segment_iter,1):segment_label(segment_iter,2),:);
       seg_point_number = size(coor_segment, 1);
       x00 = coor_segment(1,1); y00 = coor_segment(1,2);
       if major_iteration < 2 %第一轮需要对每一段初始化u
           distance_accu_seg = zeros(seg_point_number, 1);

            % 计算u
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
            u_total{segment_iter} = distance_ratio;
        end
        
        if major_iteration > 1 %除了第一轮都要先更新u
            % 目标函数计算，优化曲线长度和u的目标函数
            diffun_u{1} = @(curve_para, us, i) coor_begin(segment_iter, 1) + integral(@(t)x_integ(t, ...
            curve_para(1), curve_para(2), curve_para(3)), 0, us(1)*us(i));
            diffun_u{2} = @(curve_para, us, i) coor_begin(segment_iter, 2) + integral(@(t)y_integ(t, ...
            curve_para(1), curve_para(2), curve_para(3)), 0, us(1)*us(i));
        
            diffun_uend{1} = @(curve_para, us) coor_begin(segment_iter, 1) + integral(@(t)x_integ(t, ...
            curve_para(1), curve_para(2), curve_para(3)), 0, us(1)*1.0);
            diffun_uend{2} = @(curve_para, us) coor_begin(segment_iter, 2) + integral(@(t)y_integ(t, ...
            curve_para(1), curve_para(2), curve_para(3)), 0, us(1)*1.0);

            lu_bound = zeros(seg_point_number, 1); lu_bound(1) = seg_length_result(segment_iter) - 20; lu_bound(end) = 1;
            uu_bound = ones(seg_point_number, 1); uu_bound(1) = seg_length_result(segment_iter) + 20;
            u_seg = optimvar('u_seg', seg_point_number, 'LowerBound', lu_bound, 'UpperBound', uu_bound); % 要迭代的变量，每一段三个。
            % u_seg(0) = length; u_seg(end) = 1
            diffexpr_u = 0.0;
            curve_greek(1) = greek_result(segment_iter, 1);
            curve_greek(2) = greek_result(segment_iter, 2);
            curve_greek(3) = greek_result(segment_iter, 3);
            
            % 加入约束条件
            const = zeros(seg_point_number, seg_point_number);
            for i = 2 : seg_point_number - 1
                para = u_total{segment_iter};
                para(1) = seg_length_result(segment_iter);
                const(i, i) = 1; const(i, i+1) = -1;

                diffexpr_u = diffexpr_u + (fcn2optimexpr(diffun_u{1}, curve_greek, u_seg, i) - coor_segment(i,1)).^2;
                diffexpr_u = diffexpr_u + (fcn2optimexpr(diffun_u{2}, curve_greek, u_seg, i) - coor_segment(i,2)).^2;
            end
                diffexpr_u = diffexpr_u + 10*(fcn2optimexpr(diffun_uend{1}, curve_greek, u_seg) - coor_segment(end,1)).^2;
                diffexpr_u = diffexpr_u + 10*(fcn2optimexpr(diffun_uend{2}, curve_greek, u_seg) - coor_segment(end,2)).^2;

            ssqprob_u = optimproblem('Objective', diffexpr_u);         
            ssqprob_u.Constraints.cons = const * u_seg <= 0;
            options_u = optimoptions('lsqlin');
            options_u.Display = 'iter';
            options_u.ConstraintTolerance = 1.0000e-4;
            % 赋初值用上次迭代的数据
            x_.u_seg = u_total{segment_iter};
            x_.u_seg(1) = seg_length_result(segment_iter);
            [sol_u, fval_u, exitflag_u, output_u] = solve(ssqprob_u, x_, 'Options', options_u); 
            u_total{segment_iter} = sol_u.u_seg;
            u_total{segment_iter}(1) = 0.0000;
        end
        
        u = u_total{segment_iter};

        % 目标函数计算，优化曲线参数的目标函数
        diffun{segment_iter*2-1} = @(greek, uu) x00 + integral(@(t)x_integ(t, ...
        greek(segment_iter*3-2), greek(segment_iter*3-1), greek(segment_iter*3)), 0, uu*curve_length_seg);
        diffun{segment_iter*2} = @(greek, uu) y00 + integral(@(t)y_integ(t, ...
        greek(segment_iter*3-2), greek(segment_iter*3-1), greek(segment_iter*3)), 0, uu*curve_length_seg);

        for i = 1 : seg_point_number - 1
            diffexpr = diffexpr + (fcn2optimexpr(diffun{segment_iter*2-1}, greek, u(i)) - coor_segment(i,1)).^2;
            diffexpr = diffexpr + (fcn2optimexpr(diffun{segment_iter*2}, greek, u(i)) - coor_segment(i,2)).^2;
        end
            % 此处可以加入权重，每一段最后一点的权重，用于“0阶平滑”
             diffexpr = diffexpr + 10*(fcn2optimexpr(diffun{segment_iter*2-1}, greek, u(seg_point_number)) - coor_segment(seg_point_number,1)).^2;
             diffexpr = diffexpr + 10*(fcn2optimexpr(diffun{segment_iter*2}, greek, u(seg_point_number)) - coor_segment(seg_point_number,2)).^2;    
            
        % 方向平滑约束， “一阶平滑”
        % 前一段结尾的方向与后一段开始的方向一致
        if segment_iter ~= 1
            smoothfun{segment_iter-1} = @(greek, length) cos(greek(segment_iter*3-5)+greek(segment_iter*3-4)*length+0.5*greek(segment_iter*3-3)*length*length)- ...
                cos(greek(segment_iter*3-2));
            diffexpr = diffexpr + 1*((fcn2optimexpr(smoothfun{segment_iter-1}, greek, seg_length_result(segment_iter-1))).^2);
        end

        coor_begin(segment_iter, 1) = x00; coor_begin(segment_iter, 2) = y00;
        seg_length_result(segment_iter, 1) = curve_length_seg;
    end

    % solver生成
    ssqprob = optimproblem('Objective', diffexpr);
    options = optimoptions('lsqlin');
    options.Display = 'iter';
    options_u.ConstraintTolerance = 1.0000e-4;
    [sol, fval, exitflag, output] = solve(ssqprob, x0, 'Options', options); 
    sol.greek

    for segment_iter = 1 : 2
        greek_result(segment_iter, 1) = sol.greek(3*segment_iter-2);
        greek_result(segment_iter, 2) = sol.greek(3*segment_iter-1);
        greek_result(segment_iter, 3) = sol.greek(3*segment_iter-0);
    end
end

%% 结果对比
for segment_iter = 1 : 2
    s = linspace(0, seg_length_result(segment_iter), 500);
    coor2 = zeros(100, 2);
    coor2(1,1) = coor_begin(segment_iter, 1); coor2(1,2) = coor_begin(segment_iter, 2);
    mu0 = sol.greek(3*segment_iter-2); kappa0 = sol.greek(3*segment_iter-1); psi = sol.greek(3*segment_iter);
    for i = 1 : 500
        % 积分
        fun1 = @(t, mu0, kappa0, psi) cos(mu0 + t.*kappa0 + 0.5.*t.*t.*psi);
        fun2 = @(t, mu0, kappa0, psi) sin(mu0 + t.*kappa0 + 0.5.*t.*t.*psi);
        Q1 = integral(@(t)fun1(t,mu0,kappa0,psi),s(1),s(i));
        Q2 = integral(@(t)fun2(t,mu0,kappa0,psi),s(1),s(i));
        coor2(i,1) = Q1 + coor2(1,1);
        coor2(i,2) = Q2 + coor2(1,2);
    end
    plot(coor2(:,1),coor2(:,2),'b')
    hold on
    scatter(coor(segment_label(segment_iter,1):segment_label(segment_iter,2),1),coor(segment_label(segment_iter,1):segment_label(segment_iter,2),2), 'r')
end

abs(sol.greek(1)+seg_length_result(1)*sol.greek(2)+0.5*seg_length_result(1)*seg_length_result(1)*sol.greek(3)-sol.greek(4))*180.0/pi

