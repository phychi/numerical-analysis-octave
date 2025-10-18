function [result, T_table] = romberg_integration(f, a, b, max_iter, tol)
    % 龙贝格积分法
    % 输入: f - 被积函数句柄
    %       a, b - 积分上下限
    %       max_iter - 最大迭代次数
    %       tol - 容差
    % 输出: result - 积分结果
    %       T_table - T数表

    % 初始化T表
    T_table = zeros(max_iter, max_iter);

    % 计算T(0,0)
    h = b - a;
    T_table(1,1) = h/2 * (f(a) + f(b));

    fprintf('Romberg Integration Table:\n');
    fprintf('k\tT0\t\tT1\t\tT2\t\tT3\t\t...\n');
    fprintf('0\t%.8f\n', T_table(1,1));

    for k = 2:max_iter
        % 计算梯形公式
        n = 2^(k-2);  % 子区间数
        sum_f = 0;
        for i = 1:n
            x = a + (2*i-1)*h/2;
            sum_f = sum_f + f(x);
        end
        T_table(k,1) = 0.5 * T_table(k-1,1) + h/2 * sum_f;

        % Richardson外推
        for m = 2:k
            T_table(k,m) = (4^(m-1)*T_table(k,m-1) - T_table(k-1,m-1)) / (4^(m-1)-1);
        end

        % 输出当前行
        fprintf('%d\t', k-1);
        for m = 1:k
            fprintf('%.8f\t', T_table(k,m));
        end
        fprintf('\n');

        % 检查收敛性
        if k > 1 && abs(T_table(k,k) - T_table(k-1,k-1)) < tol
            break;
        end

        h = h / 2;
    end

    result = T_table(k,k);
end

% 测试函数定义
function y = f1(x)
    % 第一个测试函数: f(x) = x^3
    y = x.^3;
end

function y = f2(x)
    % 第二个测试函数: f(x) = sin(x)/x，处理x=0的情况
    y = zeros(size(x));
    for i = 1:length(x)
        if abs(x(i)) < 1e-10
            y(i) = 1;  % 在x=0处取极限值1
        else
            y(i) = sin(x(i)) / x(i);
        end
    end
end

function y = f3(x)
    % 第三个测试函数: f(x) = sin(x^2)
    y = sin(x.^2);
end

% 主程序
fprintf('=== Romberg Integration Method ===\n\n');

% (1) 计算 ∫_6^100 x^3 dx
fprintf('(1) ∫_6^100 x^3 dx:\n');
a1 = 6; b1 = 100;
exact1 = (b1^4 - a1^4)/4;  % 精确解
[result1, T1] = romberg_integration(@f1, a1, b1, 10, 1e-10);
fprintf('Exact result: %.8f\n', exact1);
fprintf('Romberg result: %.8f\n', result1);
fprintf('Error: %.2e\n\n', abs(result1 - exact1));

% (2) 计算 ∫_0^1 sin(x)/x dx
fprintf('(2) ∫_0^1 sin(x)/x dx:\n');
a2 = 0; b2 = 1;
[result2, T2] = romberg_integration(@f2, a2, b2, 10, 1e-10);
fprintf('Romberg result: %.8f\n\n', result2);

% (3) 计算 ∫_0^1 sin(x^2) dx
fprintf('(3) ∫_0^1 sin(x^2) dx:\n');
a3 = 0; b3 = 1;
[result3, T3] = romberg_integration(@f3, a3, b3, 10, 1e-10);
fprintf('Romberg result: %.8f\n\n', result3);
