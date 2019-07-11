clc, close all

nm1_001 = dlmread('m1_001');
nm1_026 = dlmread('m1_026');
nm1_044 = dlmread('m1_044');
nm3_001 = dlmread('m3_001');
nm3_026 = dlmread('m3_026');
nm3_044 = dlmread('m3_044');

err_001 = norm(nm1_001-nm3_001) / norm(nm1_001) * 100
err_001 = norm(nm1_026-nm3_026) / norm(nm1_026) * 100
err_001 = norm(nm1_044-nm3_044) / norm(nm1_044) * 100

