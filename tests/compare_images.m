clc;
clear all;
close all;

pkg load image;

fid = fopen('img/shepp_logan_gt.bin');
img = fread(fid, [256, 256], 'double');
fclose(fid);

fid = fopen('img/shepp_logan_noise.bin');
img_noise = fread(fid, [256, 256], 'double');
fclose(fid);

fid = fopen('img/shepp_logan_tgv_denoise.bin');
img_denoised = fread(fid, [256, 256], 'double');
fclose(fid);

#img = imread('img/shepp_logan_gt.tiff');
#img_noise = imread('img/shepp_logan_noise.tiff');
#img_denoised = imread('img/shepp_logan_tgv_denoise.tiff');

#img_denoised(img_denoised(:) > prctile(img_denoised(:),99)) = 0;

figure();
subplot(131); imagesc(img); title('Original image');
subplot(132); imagesc(img_noise); title('Noisy image');
subplot(133); imagesc(img_denoised); title('TGV image');

clear functions;
