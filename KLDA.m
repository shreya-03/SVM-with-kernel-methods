function [y] = KLDA(Data,trainlabel,gamma,type)
%Data = dlmread('arcene_train.data',' ');
    n = size(Data,1);
    d = size(Data,2);
    nm = mean(Data);
    Data = Data-repmat(nm,n,1);
%gamma = 0.5;
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
%k = 50;
    K = K - 2 * I * K + I * K * I;
%trainlabel = dlmread('arcene_train.labels',' ');
    l1 = sum(trainlabel==1);
    l2 = sum(trainlabel==-1);
    K1 = K(:,trainlabel(:,1)== 1);
    K2 = K(:,trainlabel(:,1)== -1);
    I1 = eye(l1);
    I2 = eye(l2);
    l11 = repmat(1/l1,l1,l1);
    l12 = repmat(1/l2,l2,l2);
    N = K1 * (I1-l11) * K1' + K2 * (I2-l12) * K2';
    e = 0.0001;
    N_e = N + e*eye(n);
    M1 = sum(K1,2);
    M2 = sum(K2,2);
    alphas = inv(N_e) * (M2 - M1);
    y = K * alphas;
end