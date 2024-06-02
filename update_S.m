function S = update_S(w, sim_cell, lambda, gamma, zeta, n, m, F)
%%%�������ƶȾ���
%w:Ȩ��
%sim_cell:�����ӽǵ����ƶȾ���
%n:���������
%m:��ͼ����
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
%��֤S�Ǹ�
S = (H>micro).*(H-micro);
end

