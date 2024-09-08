clc; close all; clear all;
load('ElecPosXYZ');
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
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),' ')
hold on
scatter3(Electrode_Position(:,1),Electrode_Position(:,2),Electrode_Position(:,3),' ')
text(Electrode_Position(:,1),Electrode_Position(:,2),Electrode_Position(:,3),Label_Elec);
hold off