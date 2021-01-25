clear all; clc; close all

data = dlmread('../flareTable_unscaled_Tox245K/chemTab_01.dat');

Z = reshape(data(:,1),[401 501]);
Z = squeeze(Z(1,:));

c = reshape(data(:,2),[401 501]);
c = squeeze(c(:,1));

%
close all

for ii = 3:size(data,2)
    Q = reshape(data(:,ii),[401 501]);

% figure
% plot(Q(:),'.')

figure
contourf(Z,c,Q,20,'EdgeColor','None'); colormap jet;
% xlim([0 0.1])
end 

%%
% clear all; clc; close all

fln = '../flareTable_scaled_Tox245K/flare.tbl';
nchem = dlmread(fln,'\t',[1 0 1 0]);

dims = dlmread(fln,'\t',[2+nchem 0 2+nchem 4]);
NZ = dims(1);
Nc = dims(2);
NgZ = dims(3);
Ngc = dims(4);
NgZc = dims(5);

data = dlmread(fln,'\t',3+nchem,0);
data_scal = data(1:NgZc*Ngc*NgZ*Nc*NZ,:);
data_yis = data(NgZc*Ngc*NgZ*Nc*NZ+2:end,:);

Z = reshape(data_scal(:,1),[NgZc Ngc NgZ Nc NZ]);
Z = squeeze(Z(1,1,1,1,:));

c = reshape(data_scal(:,2),[NgZc Ngc NgZ Nc NZ]);
c = squeeze(c(1,1,1,:,1));

gZ = reshape(data_scal(:,3),[NgZc Ngc NgZ Nc NZ]);
gZ = squeeze(gZ(1,1,:,1,1));

gc = reshape(data_scal(:,4),[NgZc Ngc NgZ Nc NZ]);
gc = squeeze(gc(1,:,1,1,1));

gZc = reshape(data_scal(:,5),[NgZc Ngc NgZ Nc NZ]);
gZc = squeeze(gZc(:,1,1,1,1));

%
for ii = 6:13
    Q = reshape(data_scal(:,ii),[NgZc Ngc NgZ Nc NZ]);

% figure
% plot(Q(:),'.')

Q2d = squeeze(Q(1,2,2,:,:));
figure
contourf(Z,c,Q2d,20,'EdgeColor','None'); colormap jet;
% xlim([0 0.1])
end 
