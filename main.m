function C = main(w, sim_cell, lambda, gamma, beta, n, m, num_cluster, S)
%%%计算指标值
times = 0;
while(times <= 50) %达到迭代次数或目标值不再变化
    %固定w和S，更新F
    [F, beta, flag] = update_F(S, num_cluster, beta);
    %相似度矩阵的连通分量个数满足条件时结束循环
    if flag
        break;
    end
    %固定w和F，更新S
    S = update_S(w, sim_cell, lambda, gamma, beta, n, m, F);
    %固定F和S，更新w
    w = update_w(S, sim_cell, gamma, m);
    times = times + 1;
end

%根据相似度矩阵直接得到标签
G = graph(S, 'upper');
[y, ~] = conncomp(G);
C = y';
end