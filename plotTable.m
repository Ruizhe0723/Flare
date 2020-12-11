clear all; clc; close all

data = dlmread('chemTab_01.dat');

Z = reshape(data(:,1),[401 501]); 
Z = squeeze(Z(1,:));

c = reshape(data(:,2),[401 501]); 
c = squeeze(c(:,1));

%%
close all

Q = reshape(data(:,13),[401 501]); 

figure
plot(Q(:),'.')

figure
contourf(Z,c,Q,100,'EdgeColor','None'); colormap jet;
% xlim([0 0.1])