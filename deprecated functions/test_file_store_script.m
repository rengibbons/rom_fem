clc, clear all, close all 

n = 10;
a.mydata = magic(n);
temp = a.mydata;
save('test_mat.mat','temp')

b = load('test_mat.mat');
b.temp
% rec = b.mydata