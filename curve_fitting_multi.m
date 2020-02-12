clc; clear all;
%% ����������

%% ��ȡ����

filename = '.\shortedge.txt';
[oX, oY, oZ] = textread(filename, '%f%f%f', 'delimiter', ',');

mX = mean(oX); mY = mean(oY);  mZ = mean(oZ);
sX = std(oX); sY = std(oY); sZ = std(oZ);

% ת�����X,Y,Z���꣬��������������Ļ�
X = (oX - mX); Y = (oY - mY); Z = (oZ - mZ);
total_point_number = size(X,1);  % �ܵ���
ocoor(:,1) = X; ocoor(:,2) = Y;  
coor = sortrows(ocoor,2);        % ����Y�����С���򣬴���coor�д���


%% �ֶ�
% �ֶ���
segment_number = 1; % initialization�������˶��ٶ�
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
        high = i;
        segment_label(segment_number, 1) = low;
        segment_label(segment_number, 2) = high;
        low = i;
        segment_number = segment_number + 1;
    end
    segment_label(segment_number, 1) = low;
    segment_label(segment_number, 2) = total_point_number;
end

%% �Ż�
greek_result = zeros(segment_number, 3);
seg_length_result = zeros(segment_number, 1);
coor_begin = zeros(segment_number, 2);

% Ȩ��
wt_u_endpoint = 1;        % ����uʱ���һ���Ȩ�أ�Լ�����߽�ֹ�����һ�㣬�Ͷ��ڵ���������
wt_connect_point = 1;     % �������߲���ʱ��Լ��������ֹ�����һ�㣬���ܵ���������
wt_connect_direction = 1; % ��������ʱ��Լ������һ�£��ͷֶ���-1������


% �������δ�ѭ��
for major_iteration = 1 : 2
    % ���������߽�
    l_bound = zeros(3*segment_number,1); h_bound = zeros(3*segment_number,1);
    for i = 1:segment_number
        l_bound(3*i-2) = 0.0; l_bound(3*i-1) = -0.05; l_bound(3*i-0) = -0.0025;
        h_bound(3*i-2) = 2*pi; h_bound(3*i-1) = 0.05; h_bound(3*i-0) = 0.0025;
    end
    greek = optimvar('greek', 3*segment_number, 'LowerBound', l_bound, 'UpperBound', h_bound); % Ҫ�����ı�����ÿһ������
    % ������ֵ
    x0.greek = zeros(1, 3*segment_number);
    if major_iteration < 2
        for i = 1:segment_number
            x0.greek(3*i-2) = 1.5; x0.greek(3*i-1) = 0.000; x0.greek(3*i-0) = 0.0;
        end
    else
        x0.greek = sol.greek;
    end

    % ���ֺ���
    diffexpr = 0.0;
    x_integ = @(t, mu, ka, ps) cos(mu + ka.*t + 0.5.*ps.*t.*t);
    y_integ = @(t, mu, ka, ps) sin(mu + ka.*t + 0.5.*ps.*t.*t);

    % ����ÿһ�Σ�segment_iter
    for segment_iter = 1 : segment_number
        % ȡ������
        coor_segment = coor(segment_label(segment_iter,1):segment_label(segment_iter,2),:);
        seg_point_number = size(coor_segment, 1);
        x00 = coor_segment(1,1); y00 = coor_segment(1,2); % ÿһ�ε��������
        
       %% u�ĳ�ʼ��
        % ��һ�δ�ѭ����ʼ��u
        if major_iteration < 2 
            distance_accu_seg = zeros(seg_point_number, 1);
            % ���ݾ����ʼ��u
            % ����u
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
            % ��ÿһ�ε�u
            u_total{segment_iter} = distance_ratio;
        
        else
            % ���˵�һ�ֶ�Ҫ�ȸ���u��ͨ������
            % Ŀ�꺯�����㣬�Ż����߳��Ⱥ�u��Ŀ�꺯��
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
            u_seg = optimvar('u_seg', seg_point_number, 'LowerBound', lu_bound, 'UpperBound', uu_bound); % Ҫ�����ı�����ÿһ��������
            % u_seg(0) = length; u_seg(end) = 1
            diffexpr_u = 0.0;
            % ATTENTION ����greek_result
            curve_greek(1) = greek_result(segment_iter, 1); curve_greek(2) = greek_result(segment_iter, 2); 
            curve_greek(3) = greek_result(segment_iter, 3);
            
            % ����Լ������
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
            % ����ֵ���ϴε���������
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

        % Ŀ�꺯�����㣬�Ż����߲�����Ŀ�꺯��
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
        
        % �ڵ�һ�ε����в�����ƽ��Լ��
        if major_iteration > 1
            % �˴����Լ���Ȩ�أ�ÿһ�����һ���Ȩ�أ����ڡ�0��ƽ����
            % 0��ƽ�����ܵĵ���������
            diffexpr = diffexpr + wt_connect_point*total_point_number*(fcn2optimexpr(diffun{1}, greek, 1.0) - coor_segment(seg_point_number,1)).^2;
            diffexpr = diffexpr + wt_connect_point*total_point_number*(fcn2optimexpr(diffun{2}, greek, 1.0) - coor_segment(seg_point_number,2)).^2;    
        
            if segment_iter > 1         
                % ����ƽ��Լ���� ��һ��ƽ���������ܵķֶ���������
                % ǰһ�ν�β�ķ������һ�ο�ʼ�ķ���һ��
                smoothfun{segment_iter-1} = @(greek, length) greek(segment_iter*3-5)+greek(segment_iter*3-4)*length+0.5*greek(segment_iter*3-3)*length*length- ...
                greek(segment_iter*3-2);
                % ƽ��Լ��
                diffexpr = diffexpr + wt_connect_direction*(segment_iter-1)*((fcn2optimexpr(smoothfun{segment_iter-1}, ...
                    greek, seg_length_result(segment_iter-1))).^2);
            end
        % ��һ�ε����ȳ�ʼ�����߲���
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

    % solver����
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
    
    %% �������ʾ ÿһ��major-iteration��ʾһ��
    for segment_iter = 1 : segment_number
    s = linspace(0, seg_length_result(segment_iter), 500);
    coor2 = zeros(500, 2);
    coor2(1,1) = coor_begin(segment_iter, 1); coor2(1,2) = coor_begin(segment_iter, 2);
    mu0 = sol.greek(3*segment_iter-2); kappa0 = sol.greek(3*segment_iter-1); psi = sol.greek(3*segment_iter);
    for i = 1 : 500
        % ����
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

%% ����Ա�
% abs(sol.greek(1)+seg_length_result(1)*sol.greek(2)+0.5*seg_length_result(1)*seg_length_result(1)*sol.greek(3)-sol.greek(4))*180.0/pi

fid = fopen('curve_para.txt', 'w');
for i = 1 : segment_number
   fprintf(fid, '%11.4f,%11.4f,%8.4f,%7.6f,%10.9f,%12.11f\r\n', coor_begin(i,1)+mX,coor_begin(i,2)+mY,seg_length_result(i),sol.greek(3*i-2),sol.greek(3*i-1),sol.greek(3*i-0));
end
fclose(fid);
