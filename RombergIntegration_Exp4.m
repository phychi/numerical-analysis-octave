function romberg_integration_experiment()
    % 龙贝格积分法数值实验
    % 实验1: 计算三个不同积分的近似值

    fprintf('Romberg Integration Experiment\n');
    fprintf('=================================\n\n');

    % 问题1: ∫_6^100 x^3 dx
    fprintf('Problem 1: ∫_6^100 x^3 dx\n');
    fprintf('Exact value: %.6f\n', (100^4 - 6^4)/4);
    T1 = romberg(@f1, 6, 100, 5);
    fprintf('\n');

    % 问题2: ∫_0^1 sin(x)/x dx (在x=0处定义为1)
    fprintf('Problem 2: ∫_0^1 sin(x)/x dx\n');
    T2 = romberg(@f2, 0, 1, 5);
    fprintf('\n');

    % 问题3: ∫_0^1 sin(x^2) dx
    fprintf('Problem 3: ∫_0^1 sin(x^2) dx\n');
    T3 = romberg(@f3, 0, 1, 5);

end

function T = romberg(f, a, b, n)
    % 龙贝格积分法
    % 输入:
    %   f - 被积函数句柄
    %   a - 积分下限
    %   b - 积分上限
    %   n - 迭代次数
    % 输出:
    %   T - 龙贝格T数表

    % 初始化T表
    T = zeros(n+1, n+1);

    % 计算T(0,0) - 梯形公式
    h = b - a;
    T(1,1) = h/2 * (f(a) + f(b));

    fprintf('T-Table:\n');
    fprintf('k\tT(k,0)');
    for j = 1:n
        fprintf('\t\tT(k,%d)', j);
    end
    fprintf('\n');

    fprintf('0\t%.8f', T(1,1));
    for j = 1:n
        fprintf('\t\t');
    end
    fprintf('\n');

    % 龙贝格迭代
    for k = 1:n
        % 计算复合梯形公式 T(k,0)
        sum_val = 0;
        m = 2^(k-1);
        for i = 1:m
            x = a + (2*i-1) * h/2;
            sum_val = sum_val + f(x);
        end
        T(k+1,1) = 0.5 * T(k,1) + h/2 * sum_val;

        % Richardson外推
        for j = 1:k
            T(k+1,j+1) = (4^j * T(k+1,j) - T(k,j)) / (4^j - 1);
        end

        % 输出当前行
        fprintf('%d\t%.6f', k, T(k+1,1));
        for j = 1:k
            fprintf('\t%.6f', T(k+1,j+1));
        end
        for j = k+1:n
            fprintf('\t\t');
        end
        fprintf('\n');

        h = h / 2;
    end

    fprintf('Final result: %.10f\n', T(n+1,n+1));
end

function y = f1(x)
    % 被积函数1: x^3
    y = x.^3;
end

function y = f2(x)
    % 被积函数2: sin(x)/x (在x=0处极限为1)
    y = zeros(size(x));
    for i = 1:length(x)
        if x(i) == 0
            y(i) = 1;  % x=0处的极限值
        else
            y(i) = sin(x(i)) / x(i);
        end
    end
end

function y = f3(x)
    % 被积函数3: sin(x^2)
    y = sin(x.^2);
end
