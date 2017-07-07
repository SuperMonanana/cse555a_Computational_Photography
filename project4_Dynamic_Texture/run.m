clear all
clc
n = 50;
nv = 20;
tau=500;


%%read the video
vidObj = VideoReader('Trees in the wind-SD.mp4');
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
%mov = struct('cdata',zeros(vidHeight,vidWidth,'uint8'),'colormap',[]);
mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),'colormap',[]);
k = 1;
while hasFrame(vidObj)
    mov(k).cdata = readFrame(vidObj);
    k = k+1;
end


 %Y = zeros(vidHeight*vidWidth,150);
 Y = zeros(vidHeight*vidWidth*3,250);
 t=700;
 for k = 1:250
     %J = im2double(rgb2gray(mov(k).cdata));
     
     J = im2double((mov(t).cdata));
     Y(:,k) = J(:);
     t=t+1;
 end


[x0, Ymean, Ahat, Bhat, Chat]=dytex(Y,n,nv);
I=synth(x0,Ymean,Ahat, Bhat, Chat,tau);
I(I>1)=1;
I(I<0)=0;

for k = 1: tau
    %Iee(:,:,1,k) = reshape(I(:,k),[vidHeight, vidWidth]);
    Iee(:,:,:,k) = reshape(I(:,k),[vidHeight, vidWidth,3]);     
end

gvideo = VideoWriter('challenge.mp4');
open(gvideo)
writeVideo(gvideo,Iee);
close(gvideo)
