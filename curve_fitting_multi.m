clc; clear all;
%% 多段曲线拟合

%% 读取数据

filename = '.\shortedge.txt';
[oX, oY, oZ] = textread(filename, '%f%f%f', 'delimiter', ',');

mX = mean(oX); mY = mean(oY);  mZ = mean(oZ);
sX = std(oX); sY = std(oY); sZ = std(oZ);

% 转换后的X,Y,Z坐标，三个坐标分量重心化
X = (oX - mX); Y = (oY - mY); Z = (oZ - mZ);
total_point_number = size(X,1);  % 总点数
ocoor(:,1) = X; ocoor(:,2) = Y;  
coor = sortrows(ocoor,2);        % 按照Y方向大小排序，存在coor中待用


%% 分段
% 分段数
segment_number = 1; % initialization，共分了多少段
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

% 权重
wt_u_endpoint = 1;        % 更新u时最后一点的权重，约束曲线截止于最后一点，和段内点数成正比
wt_connect_point = 1;     % 更新曲线参数时，约束曲线终止于最后一点，和总点数成正比
wt_connect_direction = 1; % 更新曲线时，约束切向一致，和分段数-1成正比


% 迭代几次大循环
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
            x0.greek(3*i-2) = 1.5; x0.greek(3*i-1) = 0.000; x0.greek(3*i-0) = 0.0;
        end
    else
        x0.greek = sol.greek;
    end

    % 积分函数
    diffexpr = 0.0;
    x_integ = @(t, mu, ka, ps) cos(mu + ka.*t + 0.5.*ps.*t.*t);
    y_integ = @(t, mu, ka, ps) sin(mu + ka.*t + 0.5.*ps.*t.*t);

    % 遍历每一段，segment_iter
    for segment_iter = 1 : segment_number
        % 取出坐标
        coor_segment = coor(segment_label(segment_iter,1):segment_label(segment_iter,2),:);
        seg_point_number = size(coor_segment, 1);
        x00 = coor_segment(1,1); y00 = coor_segment(1,2); % 每一段的起点坐标
        
       %% u的初始化
        % 第一次大循环初始化u
        if major_iteration < 2 
            distance_accu_seg = zeros(seg_point_number, 1);
            % 根据距离初始化u
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
            % 存每一段的u
            u_total{segment_iter} = distance_ratio;
        
        else
            % 除了第一轮都要先更新u，通过计算
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
            % ATTENTION 更新greek_result
            curve_greek(1) = greek_result(segment_iter, 1); curve_greek(2) = greek_result(segment_iter, 2); 
            curve_greek(3) = greek_result(segment_iter, 3);
            
            % 加入约束条件
            constr = zeros(seg_point_number, seg_point_number);
            for i = 2 : seg_point_number - 1
                constr(i, i) = 1; constr(i, i+1) = -1;

                diffexpr_u = diffexpr_u + (fcn2optimexpr(diffun_u{1}, curve_greek, u_seg, i) - coor_segment(i,1)).^2;
                diffexpr_u = diffexpr_u + (fcn2optimexpr(diffun_u{2}, curve_greek, u_seg, i) - coor_segment(i,2)).^2;
            end
                diffexpr_u = diffexpr_u + wt_u_endpoint*seg_point_number*(fcn2optimexpr(diffun_uend{1}, curve_greek, u_seg) - coor_segment(end,1)).^2;
                diffexpr_u = diffexpr_u + wt_u_endpoint*seg_point_number*(fcn2optimexpr(diffun_uend{2}, curve_greek, u_seg) - coor_segment(end,2)).^2;
            
            fprintf('Solve U\n'); fprintf('Major Iter: %d, Segment_iter: %d\n',major_iteration, segment_iter);
            ssqprob_u = optimproblem('Objective', diffexpr_u);         
            ssqprob_u.Constraints.cons = constr * u_seg <= 0;
            options_u = optimoptions('fmincon');   %fmincon problem
            options_u.Display = 'iter';
            options_u.ConstraintTolerance = 1.0000e-4;
            options_u.MaxIterations = 20;
            % 赋初值用上次迭代的数据
            x_.u_seg = u_total{segment_iter};
            x_.u_seg(1) = seg_length_result(segment_iter);
            [sol_u, fval_u, exitflag_u, output_u] = solve(ssqprob_u, x_, 'Options', options_u); 
            u_total{segment_iter} = sol_u.u_seg;
            u_total{segment_iter}(1) = 0.0000; u_total{segment_iter}(end) = 1.0000;
            seg_length_result(segment_iter) = sol_u.u_seg(1);
            curve_length_seg = seg_length_result(segment_iter);
        end
        
        u = u_total{segment_iter};
        
       %%

        % 目标函数计算，优化曲线参数的目标函数
        diffun{1} = @(greek, uu) x00 + integral(@(t)x_integ(t, ...
        greek(segment_iter*3-2), greek(segment_iter*3-1), greek(segment_iter*3)), 0, uu*curve_length_seg);
        diffun{2} = @(greek, uu) y00 + integral(@(t)y_integ(t, ...
        greek(segment_iter*3-2), greek(segment_iter*3-1), greek(segment_iter*3)), 0, uu*curve_length_seg);

        diffexpr_seg = 0.0;
    
        for i = 1 : seg_point_number
            diffexpr = diffexpr + (fcn2optimexpr(diffun{1}, greek, u(i)) - coor_segment(i,1)).^2;
            diffexpr = diffexpr + (fcn2optimexpr(diffun{2}, greek, u(i)) - coor_segment(i,2)).^2;
            
            diffexpr_seg = diffexpr_seg + (fcn2optimexpr(diffun{1}, greek, u(i)) - coor_segment(i,1)).^2;
            diffexpr_seg = diffexpr_seg + (fcn2optimexpr(diffun{2}, greek, u(i)) - coor_segment(i,2)).^2;
        end
        
        % 在第一次迭代中不加入平滑约束
        if major_iteration > 1
            % 此处可以加入权重，每一段最后一点的权重，用于“0阶平滑”
            % 0阶平滑与总的点数成正比
            diffexpr = diffexpr + wt_connect_point*total_point_number*(fcn2optimexpr(diffun{1}, greek, 1.0) - coor_segment(seg_point_number,1)).^2;
            diffexpr = diffexpr + wt_connect_point*total_point_number*(fcn2optimexpr(diffun{2}, greek, 1.0) - coor_segment(seg_point_number,2)).^2;    
        
            if segment_iter > 1         
                % 方向平滑约束， “一阶平滑”，与总的分段数成正比
                % 前一段结尾的方向与后一段开始的方向一致
                smoothfun{segment_iter-1} = @(greek, length) greek(segment_iter*3-5)+greek(segment_iter*3-4)*length+0.5*greek(segment_iter*3-3)*length*length- ...
                greek(segment_iter*3-2);
                % 平滑约束
                diffexpr = diffexpr + wt_connect_direction*(segment_iter-1)*((fcn2optimexpr(smoothfun{segment_iter-1}, ...
                    greek, seg_length_result(segment_iter-1))).^2);
            end
        % 第一次迭代先初始化曲线参数
        else
            fprintf('Init CURVE PARAMETERS segmentwise\n'); fprintf('Major Iter: %d, Segment_iter: %d\n',major_iteration, segment_iter);
            ssqprob_seg = optimproblem('Objective', diffexpr_seg);
            options_seg = optimoptions('fmincon');   % fmincon problem
            options_seg.Display = 'iter';
            options_seg.ConstraintTolerance = 1.0000e-2;
            [sol_seg, fval_seg, exitflag_seg, output_seg] = solve(ssqprob_seg, x0, 'Options', options_seg); 
            x0.greek(3*segment_iter-2) = sol_seg.greek(3*segment_iter-2);
            x0.greek(3*segment_iter-1) = sol_seg.greek(3*segment_iter-1);
            x0.greek(3*segment_iter-0) = sol_seg.greek(3*segment_iter-0);
        end
        coor_begin(segment_iter, 1) = x00; coor_begin(segment_iter, 2) = y00;
        seg_length_result(segment_iter, 1) = curve_length_seg;
    end

    % solver生成
    fprintf('Solve CURVE PARAMETERS greek\n'); fprintf('Major Iter: %d\n',major_iteration);
    ssqprob = optimproblem('Objective', diffexpr);
    options = optimoptions('fmincon');
    options.Display = 'iter';
    options_u.ConstraintTolerance = 1.0000e-4;
    [sol, fval, exitflag, output] = solve(ssqprob, x0, 'Options', options); 
    sol.greek

    for segment_iter = 1 : segment_number
        greek_result(segment_iter, 1) = sol.greek(3*segment_iter-2);
        greek_result(segment_iter, 2) = sol.greek(3*segment_iter-1);
        greek_result(segment_iter, 3) = sol.greek(3*segment_iter-0);
    end
    
    %% 结果的显示 每一个major-iteration显示一次
    for segment_iter = 1 : segment_number
    s = linspace(0, seg_length_result(segment_iter), 500);
    coor2 = zeros(500, 2);
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
    
    if major_iteration == 1
        plot(coor2(:,1),coor2(:,2),'b')
    elseif major_iteration == 2
        plot(coor2(:,1),coor2(:,2),'g') 
    elseif major_iteration == 3
        plot(coor2(:,1),coor2(:,2),'k')
    else
        plot(coor2(:,1),coor2(:,2),'c')
    end
    hold on
    scatter(coor(segment_label(segment_iter,1):segment_label(segment_iter,2),1),coor(segment_label(segment_iter,1):segment_label(segment_iter,2),2), 'r')
    xlim([-100 100])
    ylim([-100 100])
    end
    
end

%% 结果对比
% abs(sol.greek(1)+seg_length_result(1)*sol.greek(2)+0.5*seg_length_result(1)*seg_length_result(1)*sol.greek(3)-sol.greek(4))*180.0/pi

fid = fopen('curve_para.txt', 'w');
for i = 1 : segment_number
   fprintf(fid, '%11.4f,%11.4f,%8.4f,%7.6f,%10.9f,%12.11f\r\n', coor_begin(i,1)+mX,coor_begin(i,2)+mY,seg_length_result(i),sol.greek(3*i-2),sol.greek(3*i-1),sol.greek(3*i-0));
end
fclose(fid);
