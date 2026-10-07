clc;
clear all;
close all;

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load image;
end

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

  f=figure(Position=[100 100 1000 1000]);
  t=tiledlayout(2,2, TileSpacing='tight', Padding='tight');
  nexttile; imagesc(img); title('Original image'); axis off; axis image;
  nexttile; imagesc(img_noise); title('Noisy image'); axis off; axis image;
  nexttile; imagesc(img_tv); title('TV image'); axis off; axis image;
  nexttile; imagesc(img_tgv); title('TGV image'); axis off; axis image;

  colormap parula;

  exportgraphics(f,[example{1}, '.png'],"Resolution",600);
  
end

clear functions;
