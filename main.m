%% 《数值分析》非线性方程求解实验
%% 实验目标：比较不同非线性方程求解方法的性能和精度

clear all;
close all;
clc;
format long
fprintf('=== Numerical Analysis: Nonlinear Equations Solving ===\n\n');

%% 1. 二分法求解 sin(x) - x^2/2 = 0
fprintf('1. Bisection Method for sin(x) - x^2/2 = 0\n');
fprintf('--------------------------------------------\n');

% 定义函数 - 使用元素级运算
f1 = @(x) sin(x) - x.^2/2;

% 设置参数
a = 1;
b = 2;
zeta = 0.5 * 10^(-5);
max_iter = 100;

% 二分法实现
fprintf('Initial interval: [%.1f, %.1f]\n', a, b);
fprintf('Tolerance: %.2e\n', zeta);
fprintf('Iteration process:\n');

iter = 0;
while (b - a) > zeta && iter < max_iter
    c = (a + b) / 2;
    fc = f1(c);
    fa = f1(a);

    fprintf('Iter %2d: a=%f, b=%f, c=%f, f(c)=%f\n', ...
            iter, a, b, c, fc);

    if fa * fc < 0
        b = c;
    else
        a = c;
    end

    iter = iter + 1;
end

root_bisection = (a + b) / 2;
fprintf('Final root: %f\n', root_bisection);
fprintf('Function value at root: %e\n\n', f1(root_bisection));

%% 2. 牛顿法求解三个不同方程
fprintf('2. Newton''s Method\n');
fprintf('------------------\n');

% 2.1 xe^x - 1 = 0
fprintf('(1) x*exp(x) - 1 = 0\n');
f2_1 = @(x) x .* exp(x) - 1;  % 使用点乘
df2_1 = @(x) exp(x) + x .* exp(x); % 导数，使用点乘
x0 = 0.5;
tolerance = 1e-8;
max_iter_newton = 50;

fprintf('Initial guess: %.1f\n', x0);
x_current = x0;
for i = 1:max_iter_newton
    fx = f2_1(x_current);
    dfx = df2_1(x_current);
    x_next = x_current - fx / dfx;

    fprintf('Iter %2d: x=%.10f, f(x)=%.2e\n', i, x_current, fx);

    if abs(x_next - x_current) < tolerance
        break;
    end
    x_current = x_next;
end
fprintf('Final root: %.10f\n\n', x_current);

% 2.2 x^3 - x - 1 = 0
fprintf('(2) x^3 - x - 1 = 0\n');
f2_2 = @(x) x.^3 - x - 1;  % 使用点乘
df2_2 = @(x) 3*x.^2 - 1; % 导数
x0 = 1;

fprintf('Initial guess: %.1f\n', x0);
x_current = x0;
for i = 1:max_iter_newton
    fx = f2_2(x_current);
    dfx = df2_2(x_current);
    x_next = x_current - fx / dfx;

    fprintf('Iter %2d: x=%.10f, f(x)=%.2e\n', i, x_current, fx);

    if abs(x_next - x_current) < tolerance
        break;
    end
    x_current = x_next;
end
fprintf('Final root: %.10f\n\n', x_current);

% 2.3 (x-1)^2(2x-1) = 0
fprintf('(3) (x-1)^2(2x-1) = 0\n');
f2_3 = @(x) (x-1).^2 .* (2*x-1);  % 全部使用点乘
df2_3 = @(x) 2*(x-1).*(2*x-1) + 2*(x-1).^2; % 导数

% 测试两个不同的初始值
initial_guesses = [0.45, 0.65];
for guess_idx = 1:length(initial_guesses)
    x0 = initial_guesses(guess_idx);
    fprintf('Initial guess: %.2f\n', x0);

    x_current = x0;
    for i = 1:max_iter_newton
        fx = f2_3(x_current);
        dfx = df2_3(x_current);

        % 避免除零
        if abs(dfx) < 1e-12
            fprintf('Warning: Derivative too small at iteration %d\n', i);
            break;
        end

        x_next = x_current - fx / dfx;

        fprintf('Iter %2d: x=%.10f, f(x)=%.2e\n', i, x_current, fx);

        if abs(x_next - x_current) < tolerance
            break;
        end
        x_current = x_next;
    end
    fprintf('Final root: %.10f\n\n', x_current);
end

%% 3. 割线法求解 xe^x - 1 = 0
fprintf('3. Secant Method for x*exp(x) - 1 = 0\n');
fprintf('--------------------------------------\n');

f3 = @(x) x .* exp(x) - 1;  % 使用点乘
x0 = 0.4;
x1 = 0.6;
tolerance_secant = 1e-8;
max_iter_secant = 50;

fprintf('Initial points: x0=%.1f, x1=%.1f\n', x0, x1);
fprintf('Iteration process:\n');

x_prev = x0;
x_curr = x1;
f_prev = f3(x_prev);

for i = 1:max_iter_secant
    f_curr = f3(x_curr);

    % 割线法公式
    if abs(f_curr - f_prev) < 1e-12
        fprintf('Warning: Function values too close at iteration %d\n', i);
        break;
    end

    x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);

    fprintf('Iter %2d: x=%.10f, f(x)=%.2e\n', i, x_curr, f_curr);

    if abs(x_next - x_curr) < tolerance_secant
        break;
    end

    x_prev = x_curr;
    x_curr = x_next;
    f_prev = f_curr;
end

fprintf('Final root: %.10f\n\n', x_curr);

%% 4. 改进的牛顿法求解 (x-1)^2(2x-1) = 0
fprintf('4. Modified Newton''s Method for (x-1)^2(2x-1) = 0\n');
fprintf('---------------------------------------------------\n');

f4 = @(x) (x-1).^2 .* (2*x-1);  % 使用点乘
df4 = @(x) 2*(x-1).*(2*x-1) + 2*(x-1).^2; % 导数
x0 = 0.55;

fprintf('Initial guess: %.2f\n', x0);
fprintf('Comparison with standard Newton method:\n\n');

% 改进的牛顿法
fprintf('Modified Newton (multiplicity 2):\n');
x_current_mod = x0;
for i = 1:max_iter_newton
    fx = f4(x_current_mod);
    dfx = df4(x_current_mod);

    % 改进的牛顿法公式（假设重根阶数为2）
    x_next = x_current_mod - 2 * fx / dfx;

    fprintf('Iter %2d: x=%.10f, f(x)=%.2e\n', i, x_current_mod, fx);

    if abs(x_next - x_current_mod) < tolerance
        break;
    end
    x_current_mod = x_next;
end
fprintf('Final root (modified): %.10f\n\n', x_current_mod);

% 标准牛顿法比较
fprintf('Standard Newton for comparison:\n');
x_current_std = x0;
for i = 1:max_iter_newton
    fx = f4(x_current_std);
    dfx = df4(x_current_std);

    x_next = x_current_std - fx / dfx;

    fprintf('Iter %2d: x=%.10f, f(x)=%.2e\n', i, x_current_std, fx);

    if abs(x_next - x_current_std) < tolerance
        break;
    end
    x_current_std = x_next;
end
fprintf('Final root (standard): %.10f\n\n', x_current_std);

%% 5. 拟牛顿法求解非线性方程组
fprintf('5. Quasi-Newton Method for Nonlinear System\n');
fprintf('-------------------------------------------\n');

% 定义方程组 - 这里使用标量运算，不需要点乘
F = @(x) [x(1)*x(2) - x(3)^2 - 1;
          x(1)*x(2)*x(3) + x(2)^2 - x(1)^2 - 2;
          exp(x(1)) + x(3) - exp(x(2)) - 3];

% 初始猜测
x0 = [1; 1; 1];
fprintf('Initial guess: [%.1f, %.1f, %.1f]\n', x0);
fprintf('Iteration process:\n');

max_iter_quasi = 50;
tolerance_quasi = 1e-8;

% 使用数值微分计算雅可比矩阵的近似
h = 1e-6; % 数值微分步长

x_current = x0;
n = length(x0);

% 初始近似雅可比矩阵（使用单位矩阵）
B = eye(n);

for iter = 1:max_iter_quasi
    F_current = F(x_current);

    fprintf('Iter %2d: x = [%.8f, %.8f, %.8f], ', iter, x_current);
    fprintf('||F(x)|| = %.2e\n', norm(F_current));

    % 检查收敛
    if norm(F_current) < tolerance_quasi
        fprintf('Converged!\n');
        break;
    end

    % 求解线性系统 B * s = -F(x)
    s = -B \ F_current;

    % 更新解
    x_next = x_current + s;

    % 更新拟牛顿矩阵 (Broyden更新)
    F_next = F(x_next);
    y = F_next - F_current;

    % 避免除零
    if abs(s' * s) > 1e-12
        B = B + (y - B * s) * s' / (s' * s);
    end

    x_current = x_next;

    % 检查最大迭代次数
    if iter == max_iter_quasi
        fprintf('Maximum iterations reached.\n');
    end
end

fprintf('Final solution: [%.8f, %.8f, %.8f]\n', x_current);
fprintf('Final function values: [%.2e, %.2e, %.2e]\n\n', F(x_current));

%% 结果总结
fprintf('=== SUMMARY OF RESULTS ===\n');
fprintf('1. Bisection method found root at: %.8f\n', root_bisection);
fprintf('2. Newton method results:\n');
fprintf('   - x*exp(x)-1=0: converges to solution\n');
fprintf('   - x^3-x-1=0: converges to solution\n');
fprintf('   - (x-1)^2(2x-1)=0: sensitive to initial guess\n');
fprintf('3. Secant method: efficient alternative to Newton\n');
fprintf('4. Modified Newton: better for multiple roots\n');
fprintf('5. Quasi-Newton: effective for nonlinear systems\n');
