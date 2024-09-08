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
dipole_loc = LocMat(:,1203); % 100th dipole was chosen in previous part
dipole_nav = LocMat(:,1203)/norm(LocMat(:,100));
G_leadfield = GainMat(:,(3*1203)-2:3*1203);
Q = dipole_nav * Interictal(1,:);
M = G_leadfield*Q;
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
[dipole_Loreta dipole_Loreta_index] = max(max_Loreta);
di_nav_Loreta = max3_Loreta(3*dipole_Loreta_index -2:3*dipole_Loreta_index)/dipole_Loreta; %direction of dipole
dipole_position_Loreta = LocMat(:,dipole_Loreta_index);
dipole_real_3 = max(transpose(Q));
dipole_real = (sum(dipole_real_3(3 -2:3).^2))^0.5;
error_position_Loreta = norm(dipole_loc-dipole_position_Loreta);
error_direction_Loreta  = dipole_nav - transpose(di_nav_Loreta);