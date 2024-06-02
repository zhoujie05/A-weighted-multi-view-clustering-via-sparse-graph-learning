function C = main(w, sim_cell, lambda, gamma, beta, n, m, num_cluster, S)
%%%����ָ��ֵ
times = 0;
while(times <= 50) %�ﵽ����������Ŀ��ֵ���ٱ仯
    %�̶�w��S������F
    [F, beta, flag] = update_F(S, num_cluster, beta);
    %���ƶȾ������ͨ����������������ʱ����ѭ��
    if flag
        break;
    end
    %�̶�w��F������S
    S = update_S(w, sim_cell, lambda, gamma, beta, n, m, F);
    %�̶�F��S������w
    w = update_w(S, sim_cell, gamma, m);
    times = times + 1;
end

%�������ƶȾ���ֱ�ӵõ���ǩ
G = graph(S, 'upper');
[y, ~] = conncomp(G);
C = y';
end