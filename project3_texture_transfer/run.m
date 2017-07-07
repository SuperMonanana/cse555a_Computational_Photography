texture = imread('logo.jpg');
texture = imresize(texture, 0.25, 'bilinear');

figure

imshow(texture);

figure;

%outsize = size(texture)*3;
blocksize = 80;
overlap = round(blocksize/6);%1/6 of the block according to the paper 
err = 0.1;

%t2 = synthesize(texture,   outsize , blocksize, overlapsize,err);
Y = imagequilt(texture, blocksize, 10, overlap, err)
%Y=basicquilt(texture, blocksize, 10, overlap, err);
figure
imshow(uint8(Y))
