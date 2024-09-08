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
for i=1:21
    EP = ElecPos{i};
    Label_Elec{i} = num2str(EP.Name);
    Electrode_Position(i,:)= R*EP.XYZ;
end
dipole_loc = LocMat(:,1203); % 1203th dipole was chosen in previous part
dipole_nav = LocMat(:,1203)/norm(LocMat(:,1203));
G_leadfield = GainMat(:,(3*1203)-2:3*1203);
Q = dipole_nav * Interictal(1,:);
M = G_leadfield*Q;
t_spikes_1 = [1785 5452 6322 8325 9439];
spikes = zeros(21,1);
for i=1:21
 spikes(i,1) = abs((1/35)*(sum(M(i,1785-3:1785+3))+sum(M(i,5452-3:5452+3))+sum(M(i,6322-3:6322+3))+sum(M(i,8325-3:8325+3))+sum(M(i,9439-3:9439+3))));
end
S = repmat(50,21,1);
normalized_spikes = (spikes-mean(spikes))/std(spikes);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),' ')
hold on
scatter3(Electrode_Position(:,1),Electrode_Position(:,2),Electrode_Position(:,3),' ')
text(Electrode_Position(:,1),Electrode_Position(:,2),Electrode_Position(:,3),Label_Elec);
hold on
scatter3(Electrode_Position(:,1),Electrode_Position(:,2),Electrode_Position(:,3),S,spikes,'filled')
colorbar
text(Electrode_Position(:,1),Electrode_Position(:,2),Electrode_Position(:,3),Label_Elec)
hold on
plot3([dipole_loc(1,1) dipole_loc(1,1)+dipole_nav(1,1)],[dipole_loc(2,1) dipole_loc(2,1)+dipole_nav(2,1)],[dipole_loc(3,1) dipole_loc(3,1)+dipole_nav(3,1)],' ');
hold off