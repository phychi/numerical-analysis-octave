% 数值分析实验：曲线拟合的最小二乘法
% 自己实现一次、二次和三次多项式拟合

clear; clc; close all;

% 实验数据
x = [3, 4, 5, 6, 7, 8, 9];
y = [2.01, 2.98, 3.50, 5.02, 5.47, 6.02, 7.05];

% 显示原始数据
fprintf('原始数据:\n');
fprintf('x = '); fprintf('%6.2f ', x); fprintf('\n');
fprintf('y = '); fprintf('%6.2f ', y); fprintf('\n\n');

%% 最小二乘法拟合函数实现
function coeffs = least_squares_fit(x, y, n)
    % 最小二乘法多项式拟合
    % 输入：x, y - 数据点，n - 多项式次数
    % 输出：coeffs - 多项式系数（从高次到低次）
    
    m = length(x);
    
    % 构建法方程组的系数矩阵
    A = zeros(n+1, n+1);
    b = zeros(n+1, 1);
    
    % 计算法方程组的系数
    for i = 1:n+1
        for j = 1:n+1
            A(i,j) = sum(x.^(i+j-2));
        end
        b(i) = sum(y .* x.^(i-1));
    end
    
    % 解法方程组得到系数
    coeffs = A \ b;
    coeffs = flipud(coeffs); % 转换为从高次到低次排列
end

%% 计算拟合值和误差
function [y_fit, error, R_squared] = calculate_fit(x, y, coeffs)
    % 计算拟合值、误差和决定系数R²
    n = length(coeffs) - 1;
    y_fit = polyval(coeffs, x);
    error = y - y_fit;
    SSE = sum(error.^2);  % 残差平方和
    SST = sum((y - mean(y)).^2);  % 总平方和
    R_squared = 1 - SSE/SST;  % 决定系数
end

%% 一次多项式拟合（线性拟合）
fprintf('=== 一次多项式拟合（线性拟合） ===\n');
coeffs1 = least_squares_fit(x, y, 1);
[y_fit1, error1, R_sq1] = calculate_fit(x, y, coeffs1);

fprintf('拟合方程: y = %.4f*x + %.4f\n', coeffs1(1), coeffs1(2));
fprintf('残差平方和 SSE = %.6f\n', sum(error1.^2));
fprintf('决定系数 R² = %.6f\n\n', R_sq1);

%% 二次多项式拟合
fprintf('=== 二次多项式拟合 ===\n');
coeffs2 = least_squares_fit(x, y, 2);
[y_fit2, error2, R_sq2] = calculate_fit(x, y, coeffs2);

fprintf('拟合方程: y = %.4f*x² + %.4f*x + %.4f\n', coeffs2(1), coeffs2(2), coeffs2(3));
fprintf('残差平方和 SSE = %.6f\n', sum(error2.^2));
fprintf('决定系数 R² = %.6f\n\n', R_sq2);

%% 三次多项式拟合
fprintf('=== 三次多项式拟合 ===\n');
coeffs3 = least_squares_fit(x, y, 3);
[y_fit3, error3, R_sq3] = calculate_fit(x, y, coeffs3);

fprintf('拟合方程: y = %.4f*x³ + %.4f*x² + %.4f*x + %.4f\n', ...
        coeffs3(1), coeffs3(2), coeffs3(3), coeffs3(4));
fprintf('残差平方和 SSE = %.6f\n', sum(error3.^2));
fprintf('决定系数 R² = %.6f\n\n', R_sq3);

%% 结果对比表格
fprintf('=== 拟合结果对比 ===\n');
fprintf('次数    SSE          R²          拟合方程\n');
fprintf('一次    %.6f    %.6f    y=%.4fx+%.4f\n', sum(error1.^2), R_sq1, coeffs1(1), coeffs1(2));
fprintf('二次    %.6f    %.6f    y=%.4fx²+%.4fx+%.4f\n', sum(error2.^2), R_sq2, coeffs2(1), coeffs2(2), coeffs2(3));
fprintf('三次    %.6f    %.6f    y=%.4fx³+%.4fx²+%.4fx+%.4f\n', sum(error3.^2), R_sq3, coeffs3(1), coeffs3(2), coeffs3(3), coeffs3(4));


%% 子图1：拟合曲线对比
subplot(2,1,1);
xx = linspace(min(x)-0.5, max(x)+0.5, 100);
plot(x, y, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '原始数据');
hold on;
plot(xx, polyval(coeffs1, xx), 'b-', 'LineWidth', 2.5, 'DisplayName', '一次拟合');
plot(xx, polyval(coeffs2, xx), 'g-', 'LineWidth', 2.5, 'DisplayName', '二次拟合');
plot(xx, polyval(coeffs3, xx), 'm-', 'LineWidth', 2.5, 'DisplayName', '三次拟合');

% 设置坐标轴标签和标题的字体大小
xlabel('x', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 16, 'FontWeight', 'bold');
title('最小二乘法拟合对比', 'FontSize', 18, 'FontWeight', 'bold');

% 设置图例字体大小
legend('show', 'FontSize', 14, 'Location', 'best');
grid on;

% 设置坐标轴刻度字体大小
set(gca, 'FontSize', 14);

% 子图2：残差分析
subplot(2,1,2);
plot(x, error1, 'bo-', 'MarkerSize', 8, 'LineWidth', 2.5, 'DisplayName', '一次拟合残差');
hold on;
plot(x, error2, 'go-', 'MarkerSize', 8, 'LineWidth', 2.5, 'DisplayName', '二次拟合残差');
plot(x, error3, 'mo-', 'MarkerSize', 8, 'LineWidth', 2.5, 'DisplayName', '三次拟合残差');
plot([min(x), max(x)], [0, 0], 'k--', 'LineWidth', 1.5);

% 设置坐标轴标签和标题的字体大小
xlabel('x', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('残差', 'FontSize', 16, 'FontWeight', 'bold');
title('拟合残差对比', 'FontSize', 18, 'FontWeight', 'bold');

% 设置图例字体大小
legend('show', 'FontSize', 14, 'Location', 'best');
grid on;

% 设置坐标轴刻度字体大小
set(gca, 'FontSize', 14);