function [ x ] = KPCA(Data, k, gamma,type)
%Data = dlmread('arcene_train.data',' ');
    n = size(Data,1);
    d = size(Data,2);
    nm = mean(Data);
    Data = Data-repmat(nm,n,1);
%    gamma = 0.5;
    for i = 1:n
        K(i,i) = 1;
        for j = 1:i-1
            if strcmp(type,'rbf')
                K(i,j) = exp(-gamma * norm(Data(i,:)-Data(j,:))^2);
            else
                K(i,j) = Data(i,:) * Data(j,:)';
            end
            K(j,i) = K(i,j);
        end
    end
    I = repmat(1/n,n,n);
%    k = 10;
    K = K - 2 * I * K + I * K * I;
    [alphas,lamdas] = eigs(K,k);
    lamdas = diag(lamdas);
    for j = 1:k
        for i = 1:n
            alphas(i,j) = alphas(i,j)/sqrt(lamdas(j,:));
        end
    end
    lamdas = sort(lamdas,1,'descend');
    x = alphas' * K;
    x = x';
end
