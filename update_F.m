function [FSort, beta, flag] = update_F(S, num_cluster, beta)
%使用特征值分解更新矩阵F
%阈值（MATLAB中极小值）
minVal = 10e-11;
n = size(S, 1);
A = (S + S') / 2;
D= full(sparse(1:n ,1:n, sum(A)+eps));
% L = eye(n)-(D^(-1/2) * S * D^(-1/2));
L = D - A;
[F, V] = eigs(L, num_cluster+1, 'SM');
%对特征值和特征向量排序
% [VSort, index] = sort(abs(diag(V)));
% FSort = F(:, index);
FSort = F(:, 1:num_cluster);
num = sum(diag(V) < minVal);
flag = false;
%更新超参beta
if num > num_cluster
    beta = beta/2;
elseif num < num_cluster
    beta = beta*2;
else
    flag = true;
end
end

