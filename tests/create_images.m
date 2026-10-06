clc;
clear all;
close all;

pkg load image;

%% Shepp-Logan phantom
img = phantom();
img_noise = imnoise(img, 'gaussian', 0, 0.005);

fid = fopen('img/shepp_logan_gt.bin', 'wb');
fwrite(fid, size(img,1), 'double');
fwrite(fid, size(img,2), 'double');
fwrite(fid, img, 'double');
fclose(fid);

fid = fopen('img/shepp_logan_noise.bin', 'wb');
fwrite(fid, size(img_noise,1), 'double');
fwrite(fid, size(img_noise,2), 'double');
fwrite(fid, img_noise, 'double');
fclose(fid);

%% Ramp image
img = zeros(256, 256);
for j = 1:size(img,2)
    img(:,j) = (j-1)/255;
end
for j = 1+64:size(img,2)-64
    img(1+64:end-64,j) = 1 - (j-1 - 64)/127;
end
img_noise = imnoise(img, 'gaussian', 0, 0.005);

fid = fopen('img/ramp_gt.bin', 'wb');
fwrite(fid, size(img, 1), 'double');
fwrite(fid, size(img, 1), 'double');
fwrite(fid, img, 'double');
fclose(fid);

fid = fopen('img/ramp_noise.bin', 'wb');
fwrite(fid, size(img_noise,1), 'double');
fwrite(fid, size(img_noise,2), 'double');
fwrite(fid, img_noise, 'double');
fclose(fid);
%imwrite(img, 'img/shepp_logan_gt.tiff');
%imwrite(img_noise, 'img/shepp_logan_noise.tiff');

%img_load = imread('img/shepp_logan_gt.tiff');
%img_noise_load = imread('img/shepp_logan_noise.tiff');
