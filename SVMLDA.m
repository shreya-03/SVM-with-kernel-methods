clear all;
Data = dlmread('arcene_train.data',' ');
k = 10;
gamma = 0.5;
type = 'linear';
TrLabels = dlmread('arcene_train.labels',' ');
x = KLDA(Data, TrLabels, gamma,type);
svmobj = svmtrain(x,TrLabels,'kernel_function','rbf','rbf_sigma',2);
Valid_data = dlmread('arcene_valid.data',' ');
ValidLabels = dlmread('arcene_valid.labels',' ');
valid_x = KLDA(Valid_data, ValidLabels, gamma,type);
Labels = svmclassify(svmobj,valid_x);
accuracy = mean(Labels == ValidLabels);