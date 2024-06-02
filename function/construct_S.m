function S = construct_S(X, k)
%%%����K���ں͸�˹������������������ƶ�
%K:���ڸ���
%X:���ݾ���ÿ�д���һ������[n, d]
n = size(X, 1);

%���ݹ�һ��
% A = mapminmax(X', 0, 1);
% A = A';

%�����ڽӾ��� [n, n]
dis = L2_distance_1(X', X');
%Ѱ����С��K��������
[~, I] = sort(dis, 2);
%�������ƶȾ���
S = zeros(n);
for i = 1:n
    id = I(i,2:k+2);
    di = dis(i, id);
    S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
S = (S+S')/2;
end

function d = L2_distance_1(a,b)
% a,b: two matrices. each column is a data
% d:   distance matrix of a and b

if (size(a,1) == 1)
  a = [a; zeros(1,size(a,2))]; 
  b = [b; zeros(1,size(b,2))]; 
end

aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
d = repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab;

d = real(d);
d = max(d,0);
end