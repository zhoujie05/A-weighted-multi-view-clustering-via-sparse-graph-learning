function S = update_S(w, sim_cell, lambda, gamma, zeta, n, m, F)
%%%更新相似度矩阵
%w:权重
%sim_cell:各个视角的相似度矩阵
%n:样本点个数
%m:视图个数
H = zeros(n, n);
for i = 1:m
    H = H + power(w(i), gamma)*sim_cell{i};
end
V = L2_distance_1(F', F');
sum_w = sum(power(w, gamma));
H = H / sum_w;
H = H - zeta/(2*sum_w)*V;
micro = lambda / (2*sum_w);
% S=(H>micro).*(H-micro)+(H<-micro).*(H+micro);
%保证S非负
S = (H>micro).*(H-micro);
end

