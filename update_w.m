    function w = update_w(S, sim_cell, gamma, m)
%%%更新权重
%S:相似度矩阵
%sim_cell:各个视角的相似度矩阵
frobe = zeros(1, m);
for i = 1:m
    frobe(i) = power(norm(S - sim_cell{i}, 'fro'), 2);
end
w = power(gamma*frobe, 1/(1-gamma));
w = w / sum(w);
end

