clc; close all; clear all;
load('ElecPosXYZ');
load('Interictal.mat')
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];
Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;
R = 9.2; % Radius of outter layer of heads sphere
for i=1:21
    EP = ElecPos{i};
    Label_Elec{i} = num2str(EP.Name);
    Electrode_Position(i,:)= R*EP.XYZ;
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
for i=1:21 
figure(1)
    subplot(7,3,i);
    plot(M(i,:));
    title(Label_Elec{i});
end
t_spikes = [1750 5372 8302 9361];
spikes = zeros(21,1);
for i=1:21
 spikes(i,1) = abs((1/28)*(sum(M(i,1750-3:1750+3))+sum(M(i,5372-3:5372+3))+sum(M(i,8302-3:8302+3))+sum(M(i,8325-3:8325+3))+sum(M(i,9361-3:9361+3))));
end
figure(2);
normalized_spikes = (spikes-mean(spikes))/std(spikes);
Display_Potential_3D(R,normalized_spikes)

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

