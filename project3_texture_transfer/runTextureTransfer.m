clc;clear;close all;

targetImg = im2double(imread('DSC_0676.jpg'));
inputImg = im2double(imread('caustics.png'));

targetImg = imresize(targetImg, 0.55, 'bilinear');
inputImg = imresize(inputImg, 0.8, 'bilinear');

alpha = 0.2;
szPatch = 5;
szOverlap = 4;
ifdebug = 0;
[t1 t2]= texture_transfer6(inputImg, targetImg, alpha, szPatch, szOverlap, ifdebug);

figure(1)
imshow(inputImg);
figure(2)
imshow(targetImg);
figure(3)
imshow(t2, [])

