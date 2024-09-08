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
dipole_patch = [500:518];
dipole_dir = zeros(3,19);
for i=1:19
dipole_dir(:,i) = LocMat(:,dipole_patch(i))./norm(LocMat(:,dipole_patch(i)));
end
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),' ')
hold on
scatter3(Electrode_Position(:,1),Electrode_Position(:,2),Electrode_Position(:,3),' ')
text(Electrode_Position(:,1),Electrode_Position(:,2),Electrode_Position(:,3),Label_Elec);
for i=1:19
    hold on
    scatter3(LocMat(1,dipole_patch(i)),LocMat(2,dipole_patch(i)),LocMat(3,dipole_patch(i)),' ')
    hold on
    plot3([LocMat(1,dipole_patch(i)) LocMat(1,dipole_patch(i))+dipole_dir(1,i)], [LocMat(2,dipole_patch(i)) LocMat(2,dipole_patch(i))+dipole_dir(2,i)], [LocMat(3,dipole_patch(i)) LocMat(3,dipole_patch(i))+dipole_dir(3,i)],' ')
end
