clc;
clear all;
close all;

pkg load image;

img = phantom();
img_noise = imnoise(img, 'gaussian');

fid = fopen('img/shepp_logan_gt.bin', 'wb');
fwrite(fid, img, 'double');
fclose(fid);

fid = fopen('img/shepp_logan_noise.bin', 'wb');
fwrite(fid, img_noise, 'double');
fclose(fid);

%imwrite(img, 'img/shepp_logan_gt.tiff');
%imwrite(img_noise, 'img/shepp_logan_noise.tiff');

%img_load = imread('img/shepp_logan_gt.tiff');
%img_noise_load = imread('img/shepp_logan_noise.tiff');
