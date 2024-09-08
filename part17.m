clc; close all; clear all;
load('ElecPosXYZ.mat');
load('Interictal.mat')
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];
Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;
R = 9.2; % Radius of outter layer of heads sphere
for beta=1:21
    EP = ElecPos{beta};
    Label_Elec{beta} = num2str(EP.Name) ;
    Electrode_Position(beta,:)= R*EP.XYZ;
end
dipole_patch = [500:518];
dipole_dir = zeros(3,19);
for i=1:19
dipole_dir(:,i) = LocMat(:,dipole_patch(i))./norm(LocMat(:,dipole_patch(i)));
end
Q= zeros(3951,length(Interictal(1,:)));
for i=1:19
    Q(3*dipole_patch(i)-2:3*dipole_patch(i),:) = dipole_dir(:,i)*Interictal(4,:);
end
M = GainMat*Q;
spikes = zeros(21,1);
for i=1:21
 spikes(i,1) = abs((1/35)*(sum(M(i,1785-3:1785+3))+sum(M(i,5452-3:5452+3))+sum(M(i,6322-3:6322+3))+sum(M(i,8325-3:8325+3))+sum(M(i,9439-3:9439+3))));
end
a = 0.6;
IN21 = eye(21);
d = 1; %Resolution
A1 = zeros(length(LocMat(1,:)),length(LocMat(1,:))); % Defining [A1]alpha_beta
for i=1:length(LocMat(1,:))
    for j=1:length(LocMat(1,:))
        if (norm(LocMat(:,i)-LocMat(:,j))) == d
            A1(i,j) = 1/6;
        end
    end
end
p = length(LocMat(1,:));
A0 = inv(diag(A1 * ones(p,1))) * A1;
A = kron(A0,eye(3));
B = (6/(d*d))*(A-eye(3*p));
omega = zeros(p , p);
for beta=1:p
    omega_bb = 0;
    for alpha = 1:21
     omega_bb = omega_bb + GainMat(alpha,3*beta-2:3*beta)*transpose(GainMat(alpha,3*beta-2:3*beta));
    end
    omega(beta,beta) = omega_bb^0.5;
end
W_Loreta = kron(omega,eye(3))*transpose(B) * B *kron(omega,eye(3));
Q_Loreta = inv(transpose(W_Loreta)*W_Loreta)*transpose(GainMat)*inv(GainMat*inv(transpose(W_Loreta)*W_Loreta)*transpose(GainMat)+a*IN21)*M;
max3_Loreta = max(transpose(Q_Loreta)); %calculating maximum amplitude of dipoles in all samples in 3 directions
for i=1:length(LocMat(1,:))
    max_Loreta(i) = (sum(max3_Loreta(3*i -2:3*i).^2))^0.5; %calculating maximum amplitude of each dipole in all samples in one directions 
end
p = length(LocMat(1,:));
val_real = zeros(1,p);
val_real(dipole_patch) = 1;
Threshold_LORETA = 0:0.01:3;
AUC_LORETA = 0;
for i=1:length(Threshold_LORETA)
    val_LORETA = zeros(1,p);
    for j=1:p
        if (max_Loreta(j)>Threshold_LORETA(i))
            val_LORETA(j) = 1;
        end
    end
    Confusion_matrix = confusionmat(val_real,val_LORETA);
    Accuracy(i) = (Confusion_matrix(2,2)+Confusion_matrix(1,1))/(Confusion_matrix(1,2)+Confusion_matrix(1,1)+Confusion_matrix(2,2)+Confusion_matrix(2,1));
    TPR(i) = Confusion_matrix(2,2)/(Confusion_matrix(2,2)+Confusion_matrix(2,1));
    FPR(i) = Confusion_matrix(1,2)/(Confusion_matrix(1,2)+Confusion_matrix(1,1));
    precision(i) = Confusion_matrix(1,1)/(Confusion_matrix(1,1)+Confusion_matrix(1,2));
    asses(i) = (((1-TPR(i))^2)+(FPR(i)^2))^0.5;
    if i~=1
        AUC_LORETA =TPR(i)*(FPR(i-1)-FPR(i)) +AUC_LORETA ;
    end
end
plot(FPR,TPR)
title('ROC Of LORETA Algorithm')
xlabel('1-Specificity')
ylabel('Sensitivity')
[s ,Optimum_Index] = min(asses);
Threshold_Optimum_LORETA = Threshold_LORETA(Optimum_Index);
Accuracy_Optimum_LORETA = Accuracy(Optimum_Index);
precision_optimum_LORETA = precision(Optimum_Index);
FPR_LORETA = FPR(Optimum_Index);
AUC_LORETA