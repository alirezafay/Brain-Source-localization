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
dipole_loc = LocMat(:,1203); % 1203th dipole was chosen in previous part
dipole_nav = LocMat(:,1203)/norm(LocMat(:,1203));
G_leadfield = GainMat(:,(3*1203)-2:3*1203);
Q = dipole_nav * Interictal(1,:);
M = G_leadfield*Q;
spikes = zeros(21,1);
for i=1:21
 spikes(i,1) = abs((1/35)*(sum(M(i,1785-3:1785+3))+sum(M(i,5452-3:5452+3))+sum(M(i,6322-3:6322+3))+sum(M(i,8325-3:8325+3))+sum(M(i,9439-3:9439+3))));
end
a = 0.6;
IN21 = eye(21);
Q_MNE = transpose(GainMat) * inv(GainMat * transpose(GainMat) + a*IN21)* M ;
IN3 = eye(3);
omega = zeros(length(LocMat(1,:)) , length(LocMat(1,:)));
W = zeros(length(GainMat(1,:)),length(GainMat(1,:)));
omega_bb = 0;
for beta=1:length(LocMat(1,:))
    for alpha = 1:21
     omega_bb = omega_bb + GainMat(alpha,3*beta-2:3*beta)*transpose(GainMat(alpha,3*beta-2:3*beta));
    end
    omega(beta,beta) = omega_bb^0.5;
end
W = kron(omega,IN3);
Q_WMNE = inv(transpose(W)*W)*transpose(GainMat)*inv(GainMat*inv(transpose(W)*W)*transpose(GainMat)+a*IN21)*M;