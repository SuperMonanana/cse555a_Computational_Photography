clc;clear;close all;
targetImg = im2double(imread('arch.jpg'));
inputImg = im2double(imread('caustics.png'));
targetImg = imresize(targetImg, 0.6, 'bilinear');
inputImg = imresize(inputImg, 0, 'bilinear');
merge = 1;


tilesize = 64;%!!!!szPatch should not be too big
sizeout=size(targetImg);
overlap = 16;
isdebug = 1;

[imout] = imagequilting(inputImg, targetImg,sizeout, tilesize, overlap)
