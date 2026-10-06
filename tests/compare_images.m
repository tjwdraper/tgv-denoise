clc;
clear all;
close all;

pkg load image;

for example = {'shepp_logan', 'ramp'}
  filename = ['img/', example{1}, '_gt.bin'];
  fid = fopen(['img/', example{1}, '_gt.bin'], 'rb');
  rows = fread(fid, 1, 'double');
  cols = fread(fid, 1, 'double');
  img = fread(fid, [rows, cols], 'double');
  fclose(fid);

  fid = fopen(['img/', example{1}, '_noise.bin'], 'rb');
  rows = fread(fid, 1, 'double');
  cols = fread(fid, 1, 'double');
  img_noise = fread(fid, [rows, cols], 'double');
  fclose(fid);

  fid = fopen(['img/', example{1}, '_tv_denoise.bin'], 'rb');
  rows = fread(fid, 1, 'double');
  cols = fread(fid, 1, 'double');
  img_tv = fread(fid, [rows, cols], 'double');
  fclose(fid);

  fid = fopen(['img/', example{1}, '_tgv_denoise.bin'], 'rb');
  rows = fread(fid, 1, 'double');
  cols = fread(fid, 1, 'double');
  img_tgv = fread(fid, [rows, cols], 'double');
  fclose(fid);

  figure();
  subplot(221); imagesc(img); title('Original image');
  subplot(222); imagesc(img_noise); title('Noisy image');
  subplot(223); imagesc(img_tv); title('TV image');
  subplot(224); imagesc(img_tgv); title('TGV image');
end

clear functions;
