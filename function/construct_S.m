function S = construct_S(X, k)
%%%利用K近邻和高斯函数计算样本间的相似度
%K:近邻个数
%X:数据矩阵，每行代表一个样本[n, d]
n = size(X, 1);

%数据归一化
% A = mapminmax(X', 0, 1);
% A = A';

%计算邻接矩阵 [n, n]
dis = L2_distance_1(X', X');
%寻找最小的K个样本点
[~, I] = sort(dis, 2);
%生成相似度矩阵
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