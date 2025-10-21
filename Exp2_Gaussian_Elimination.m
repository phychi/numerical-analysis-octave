function gaussian_elimination_experiment()
    % 实验1：第一个线性方程组
    fprintf('=== EXPERIMENT 1 ===\n');
    A1 = [1e-8, 2, 3; 
          -1, 3.712, 4.623; 
          -2, 1.072, 5.643];
    b1 = [1; 2; 3];
    
    % 高斯消去法
    fprintf('--- Gaussian Elimination ---\n');
    x_ge1 = gaussian_elimination(A1, b1);
    fprintf('Solution: '); disp(x_ge1');
    fprintf('Residual norm: %e\n', norm(A1*x_ge1 - b1));
    
    % 高斯列主元消去法
    fprintf('--- Gaussian Elimination with Partial Pivoting ---\n');
    x_gepp1 = gaussian_elimination_pp(A1, b1);
    fprintf('Solution: '); disp(x_gepp1');
    fprintf('Residual norm: %e\n', norm(A1*x_gepp1 - b1));
    
    fprintf('\n');
    
    % 实验2：第二个线性方程组
    fprintf('=== EXPERIMENT 2 ===\n');
    A2 = [4, -2, 4;
          -2, 17, 10;
          -4, 10, 9];
    b2 = [10; 3; 7];
    
    % 高斯消去法
    fprintf('--- Gaussian Elimination ---\n');
    x_ge2 = gaussian_elimination(A2, b2);
    fprintf('Solution: '); disp(x_ge2');
    fprintf('Residual norm: %e\n', norm(A2*x_ge2 - b2));
    
    % 高斯列主元消去法
    fprintf('--- Gaussian Elimination with Partial Pivoting ---\n');
    x_gepp2 = gaussian_elimination_pp(A2, b2);
    fprintf('Solution: '); disp(x_gepp2');
    fprintf('Residual norm: %e\n', norm(A2*x_gepp2 - b2));
end

function x = gaussian_elimination(A, b)
    % 高斯消去法
    % 输入：系数矩阵A，右端向量b
    % 输出：解向量x
    
    [n, ~] = size(A);
    Ab = [A, b];  % 增广矩阵
    
    % 前向消元
    for k = 1:n-1
        % 检查主元是否为零
        if abs(Ab(k,k)) < 1e-15
            error('Zero pivot encountered at step %d', k);
        end
        
        for i = k+1:n
            factor = Ab(i,k) / Ab(k,k);
            Ab(i,k:n+1) = Ab(i,k:n+1) - factor * Ab(k,k:n+1);
        end
    end
    
    % 回代求解
    x = zeros(n,1);
    for i = n:-1:1
        x(i) = (Ab(i,n+1) - Ab(i,i+1:n)*x(i+1:n)) / Ab(i,i);
    end
end

function x = gaussian_elimination_pp(A, b)
    % 高斯列主元消去法
    % 输入：系数矩阵A，右端向量b
    % 输出：解向量x
    
    [n, ~] = size(A);
    Ab = [A, b];  % 增广矩阵
    
    % 前向消元（带列主元）
    for k = 1:n-1
        % 列主元选择：找到第k列中绝对值最大的元素
        [~, pivot_row] = max(abs(Ab(k:n, k)));
        pivot_row = pivot_row + k - 1;
        
        % 交换行
        if pivot_row ~= k
            Ab([k, pivot_row], :) = Ab([pivot_row, k], :);
        end
        
        % 消元
        for i = k+1:n
            factor = Ab(i,k) / Ab(k,k);
            Ab(i,k:n+1) = Ab(i,k:n+1) - factor * Ab(k,k:n+1);
        end
    end
    
    % 回代求解
    x = zeros(n,1);
    for i = n:-1:1
        x(i) = (Ab(i,n+1) - Ab(i,i+1:n)*x(i+1:n)) / Ab(i,i);
    end
end

% 运行实验
gaussian_elimination_experiment();