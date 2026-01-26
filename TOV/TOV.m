% 解决TOV方程的耦合系统 - 核物质和暗物质 (使用MeV/fm^3单位)

clear all;
close all;
clc;

% 物理常数 (使用国际单位制)
c = 2.99792458e8;     % m/s 国际单位制
G = 6.67430e-11;   % G in m^3/kg/s^2 
Msun_turn = 1.98847e30;  


% 读取EOS数据
eos1 = load('eos1.txt'); % 核物质的EOS [p(MeV/fm^3), rho(MeV/fm^3)]
eos2 = load('2e-12 300.txt'); % 暗物质的EOS [p(MeV/fm^3), rho(MeV/fm^3)]

% 检查EOS数据是否加载成功
if isempty(eos1) || isempty(eos2)
    error('EOS数据加载失败，请检查文件路径和内容');
end

% 从用户获取参数
rho1_min = input('请输入核物质中心能量密度的最小值(MeV/fm^3): ');
rho1_max = input('请输入核物质中心能量密度的最大值(MeV/fm^3): ');
rho1_step = input('请输入核物质中心能量密度的步长(0表示保持为最小值): ');

rho2_min = input('请输入暗物质中心能量密度的最小值(MeV/fm^3): ');
rho2_max = input('请输入暗物质中心能量密度的最大值(MeV/fm^3): ');
rho2_step = input('请输入暗物质中心能量密度的步长(0表示保持为最小值): ');

%转换为国际单位
eos1 = eos1 * 1.602e32;
eos2 = eos2 * 1.602e32;
rho1_min = rho1_min * 1.602e32;
rho1_max = rho1_max * 1.602e32;
rho1_step = rho1_step * 1.602e32;
rho2_min = rho2_min * 1.602e32;
rho2_max = rho2_max * 1.602e32;
rho2_step = rho2_step * 1.602e32;

% 处理步长为0的情况
if rho1_step == 0
    rho1_values = rho1_min;
else
    rho1_values = rho1_min:rho1_step:rho1_max;
end

if rho2_step == 0
    rho2_values = rho2_min;
else
    rho2_values = rho2_min:rho2_step:rho2_max;
end

% 初始化结果存储
results = [];

% 主循环 - 遍历所有rho1和rho2的组合
for rho1_center = rho1_values
    for rho2_center = rho2_values
        fprintf('正在计算 rho1=%.2e MeV/fm^3, rho2=%.2e MeV/fm^3...\n', rho1_center/1.602e32,rho2_center/1.602e32);
        
        [~, idx1] = min(abs(eos1(:,1) - rho1_center));
        p1_center = eos1(idx1, 2); % 提取对应的压强

        % 对 eos2 执行相同操作
        [~, idx2] = min(abs(eos2(:,1) - rho2_center));
        p2_center = eos2(idx2, 2);
        
        % 检查是否找到合理的初始压力
        if isempty(p1_center) || isempty(p2_center)
            warning('无法为rho1=%.2e或rho2=%.2e找到初始压力，跳过', rho1_center, rho2_center);
            continue;
        end

         
        % 调用求解器
        try
            [M1_final, M2_final, R1, R2] = solve_TOV(p1_center, p2_center, eos1, eos2, G,c);
            
            % 转换为太阳质量和千米
            M1_final_Msun = M1_final / Msun_turn; 
            M2_final_Msun = M2_final / Msun_turn;
            R1_km = R1 * 1e-3;  
            R2_km = R2 * 1e-3;
            
            % 存储结果
            results = [results; rho1_center/1.602e32, rho2_center/1.602e32, M1_final_Msun, M2_final_Msun, R1_km, R2_km];
            
            % 显示当前结果
            fprintf('结果: M1=%.3f M☉, M2=%.3f M☉, R1=%.2f km, R2=%.2f km\n', ...
                    M1_final_Msun, M2_final_Msun, R1_km, R2_km);
        catch ME
            warning('计算rho1=%.2e, rho2=%.2e时出错: %s', ...
                    rho1_center, rho2_center, ME.message);
        end
    end
end

% 检查是否有结果
if isempty(results)
    error('没有生成任何结果，请检查输入参数和EOS数据');
end

% 显示结果
disp('结果表格: [rho1_center (MeV/fm^3), rho2_center (MeV/fm^3), M1_final (M☉), M2_final (M☉), R1 (km), R2 (km)]');
disp(results);

% 保存结果到文件
try
    header = {'rho1_center(MeV/fm^3)', 'rho2_center(MeV/fm^3)', 'M1(Msun)', 'M2(Msun)', 'R1(km)', 'R2(km)'};
    fid = fopen('TOV_results.txt', 'w');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', header{:});
    fclose(fid);
    dlmwrite('TOV_results.txt', results, 'delimiter', '\t', 'precision', '%.6e', '-append');
    disp('结果已成功保存到 TOV_results.txt');
catch
    disp('无法写入文件，请在命令行窗口查看结果');
    disp(results);
end

% ========== TOV方程求解函数 ==========
function [M1_final, M2_final, R1, R2,step,p1_new,p2_new] = solve_TOV(p1_center, p2_center, eos1, eos2, G, c)
    % 初始条件 (避免r=0时数值不稳定)
    r_min = 1e-5;       % 很小的非零半径 (m)
    M1_initial = 1; % 很小的非零质量 (kg)
    M2_initial = 1; % 很小的非零质量 (kg)
    
    % 积分参数
    r_max = 1e10;        % 最大半径 (m) 
    dr_initial = 13;    % 初始步长 (m)
    dr_min = 1e-3;      % 最小步长 (m)
    max_steps = 1e10;    % 最大步数
    
    % 初始化变量
    r = r_min;
    M1 = M1_initial;
    M2 = M2_initial;
    p1 = p1_center;
    p2 = p2_center;

    
    % 标志位
    if p2 == 0
        p2_zero = true;
    else
        p2_zero = false;
    end
    if p1 == 0
        p1_zero = true;
    else
        p1_zero = false;
    end
    R1 = 0;%这是啥
    R2 = 0;
    
    % 主循环
    step = 0;
    while (~p1_zero || ~p2_zero)  %当至少一个p1_zero或p2_zero为false时，循环继续；若两个均为true时，循环停止
        step = step + 1;
        if p2 <1e-100 && p1 < 1e-100%当暗物质压强小于10e-50Pa,退出循环
            break;
        end

        %若pcenter一开始就是0，规定m_new就保持边界值
        if p2_zero && step == 1
            M2_new = M2;
        end
        if p1_zero && step == 1
            M1_new = M1;
        end

        if step == 1
            dr = dr_initial;
        end

        if step > 1 && p1 < 1e-30
           dr = r/1000;
        end
        
        % 当前总质量
        M = M1 + M2;
        
        % 从EOS获取能量密度（添加外推保护）原先这里似乎写错了，是根据能量密度获得的压强，进行修改
        rho1 = interp1(eos1(:,2), eos1(:,1), p1, 'spline', 0);
        rho2 = interp1(eos2(:,2), eos2(:,1), p2, 'spline', 0);
        
        % 四阶龙格-库塔法
        if ~p1_zero && ~p2_zero
            [k1_p1, k1_M1] = TOV_eqns(r, M, p1,p2, rho1, G, c);
            [k1_p2, k1_M2] = TOV_eqns(r, M, p2,p1, rho2, G, c);
            [k2_p1, k2_M1] = TOV_eqns(r+dr/2, M, p1+dr/2*k1_p1,p2+dr/2*k1_p2, rho1, G, c);
            [k2_p2, k2_M2] = TOV_eqns(r+dr/2, M, p2+dr/2*k1_p2,p1+dr/2*k1_p1, rho2, G, c);
            [k3_p1, k3_M1] = TOV_eqns(r+dr/2, M, p1+dr/2*k2_p1,p2+dr/2*k2_p2, rho1, G, c);
            [k3_p2, k3_M2] = TOV_eqns(r+dr/2, M, p2+dr/2*k2_p2,p1+dr/2*k2_p1, rho2, G, c);
            [k4_p1, k4_M1] = TOV_eqns(r+dr, M, p1+dr*k3_p1,p2+dr*k3_p2, rho1, G, c);
            [k4_p2, k4_M2] = TOV_eqns(r+dr, M, p2+dr*k3_p2,p1+dr*k3_p1, rho2, G, c);
            
            p1_new = p1 + dr/6*(k1_p1 + 2*k2_p1 + 2*k3_p1 + k4_p1);
            M1_new = M1 + dr/6*(k1_M1 + 2*k2_M1 + 2*k3_M1 + k4_M1);
            p2_new = p2 + (dr/6)*(k1_p2 + 2*k2_p2 + 2*k3_p2 + k4_p2);
            M2_new = M2 + (dr/6)*(k1_M2 + 2*k2_M2 + 2*k3_M2 + k4_M2);
            
            % 检查压力是否变为零
            if p1_new <= 0
                p1_zero = true;
                R1 = r;
                p1_new = 0;
                M1_new = M1; % 压力为零后质量不再变化
            end
             % 检查压力是否变为零
            if p2_new <= 0
                p2_zero = true;
                R2 = r;
                p2_new = 0;
                M2_new = M2; % 压力为零后质量不再变化
            end
        end

        
        
        if ~p1_zero && p2_zero%p2降到0的情况
            [k1_p1, k1_M1] = TOV_eqns(r, M, p1,0, rho1, G, c);
            [k2_p1, k2_M1] = TOV_eqns(r+dr/2, M, p1+dr/2*k1_p1,0, rho1, G, c);
            [k3_p1, k3_M1] = TOV_eqns(r+dr/2, M, p1+dr/2*k2_p1,0, rho1, G, c);
            [k4_p1, k4_M1] = TOV_eqns(r+dr, M, p1+dr*k3_p1,0, rho1, G, c);

            p1_new = p1 + dr/6*(k1_p1 + 2*k2_p1 + 2*k3_p1 + k4_p1);
            M1_new = M1 + dr/6*(k1_M1 + 2*k2_M1 + 2*k3_M1 + k4_M1);
        
            % 检查压力是否变为零
            if p1_new <= 0
                p1_zero = true;
                R1 = r;
                p1_new = 0;
                M1_new = M1; % 压力为零后质量不再变化
            end
            p2_new = 0;
        end
        
        if p1_zero && ~p2_zero%p1降到0的情况
           [k1_p2, k1_M2] = TOV_eqns(r, M, p2,0, rho2, G, c);
           [k2_p2, k2_M2] = TOV_eqns(r+dr/2, M, p2+dr/2*k1_p2,0, rho2, G, c);
           [k3_p2, k3_M2] = TOV_eqns(r+dr/2, M, p2+dr/2*k2_p2,0, rho2, G, c);
           [k4_p2, k4_M2] = TOV_eqns(r+dr, M, p2+dr*k3_p2,0, rho2, G, c); 

           p2_new = p2 + (dr/6)*(k1_p2 + 2*k2_p2 + 2*k3_p2 + k4_p2);
           M2_new = M2 + (dr/6)*(k1_M2 + 2*k2_M2 + 2*k3_M2 + k4_M2);           

        
                       % 检查压力是否变为零
           if p2_new <= 0
               p2_zero = true;
               R2 = r;
               p2_new = 0;
               M2_new = M2; % 压力为零后质量不再变化
           end
           p1_new = 0;
        end
        % 更新变量
        r = r + dr;
        p1 = p1_new;
        p2 = p2_new;
        M1 = M1_new;
        M2 = M2_new;
        
        % 如果两者都达到零压力，退出循环
        %if p1_zero && p2_zero
            %break;
        %end
    end
    
    % 最终值
    M1_final = M1;
    M2_final = M2;
    
    % 如果压力没有降到零，使用当前半径
    if ~p1_zero
        R1 = r;
    end
    if ~p2_zero
        R2 = r;
    end
end

% ========== TOV方程 ==========国际单位制
function [dpdr, dMdr] = TOV_eqns(r, M, p,pother, rho, G, c)
    if r <= 0
        dpdr = 0;
        dMdr = 0;
        return;
    end
    
    % 避免除以零
    if M <= 0
        M = 1e-20;
    end
    
    % 自然单位制下的TOV方程 
    term1 = - (G * M * rho) / (r^2 * c^2) ;
    term2 = (1 + p/ rho );
    term3 = (1 + 4*pi* (r^3) * (p + pother)/ (M * c^2) );
    term4 = (1 - 2*G*M/ (r * c^2) )^(-1);
    
    dpdr = term1 * term2 * term3 * term4;
    dMdr = (4 * pi * r^2 * rho)/c^2;
end