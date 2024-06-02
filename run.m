%最佳参数
function run(fileName)
%读取数据
data = load(fileName);
X = data.X;
Y = data.Y;
m = length(X); %视图个数
n = size(Y, 1); %样本点个数
K = 10; %近邻
num_cluster = length(unique(Y)); %样本点种类数

tic;
%生成相似度矩阵
sim_cell = cell(1, m);
for i = 1:m
    sim_cell{1, i} = construct_S(X{i}, K);
end
toc;

%初始化权重
w = ones(1, m) / m;
%初始化相似度矩阵S
S = zeros(n, n);
for i = 1:m
    S = S + w(i)*sim_cell{1, i};
end

lambda = [0.02, 0.06, 0.1, 0.14, 0.18, 0.22];
gamma = [1.1, 1.3, 1.5, 1.7, 1.9];
beta = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05];

lambda_num = size(lambda, 2);
gamma_num = size(gamma, 2);
beta_num = size(beta, 2);
accs = zeros(lambda_num, gamma_num, beta_num); %记录不同参数的指标

%寻找最佳参数
for i = 1:lambda_num
    for j = 1:gamma_num
        for k = 1:beta_num
            C = main(w, sim_cell, lambda(i), gamma(j), beta(k), n, m, num_cluster, S);
            accs(i, j, k) = cluster_acc(C, Y);
%             disp(["lambda is: ", lambda(i), "gamma is: ", gamma(j), "beta is: ", beta(k), "acc is: ", accs(i, j, k)]);
        end
    end
end

%运行多次取均值
T = 10;
nmi_value = zeros(1, T);
acc_value = zeros(1, T);
ARI_value = zeros(1, T);
[~, position] = max(accs(:));
[lambda_index, gamma_index, beta_index] = ind2sub(size(accs), position);
best_lambda = lambda(lambda_index);
best_gamma = gamma(gamma_index);
best_beta = beta(beta_index);
fprintf('The best lambda is %f\n', best_lambda);
fprintf('The best gamma is %f\n', best_gamma);
fprintf('The best beta is %f\n', best_beta);
%计时
tic;
for i = 1:T
    C = main(w, sim_cell, best_lambda, best_gamma, best_beta, n, m, num_cluster, S);
    nmi_value(i) = nmi(C, Y);
    acc_value(i) = cluster_acc(C, Y);
    ARI_value(i) = rand_index(C,Y,'adjusted');
end
disp("------------------------------------------------");
disp([fileName, "cost is:"]);
toc/T
disp("------------------------------------------------");
outcome = [mean(nmi_value), mean(acc_value), mean(ARI_value); std(nmi_value), std(acc_value), std(ARI_value)];
save(fileName(1:length(fileName)-4), 'outcome');
end