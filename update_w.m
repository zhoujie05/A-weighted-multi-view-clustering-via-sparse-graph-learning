    function w = update_w(S, sim_cell, gamma, m)
%%%����Ȩ��
%S:���ƶȾ���
%sim_cell:�����ӽǵ����ƶȾ���
frobe = zeros(1, m);
for i = 1:m
    frobe(i) = power(norm(S - sim_cell{i}, 'fro'), 2);
end
w = power(gamma*frobe, 1/(1-gamma));
w = w / sum(w);
end

