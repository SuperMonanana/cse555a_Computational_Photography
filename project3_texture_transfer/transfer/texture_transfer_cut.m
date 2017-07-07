function [outRGB Y] = texture_transfer_cut(inputImg, tarImg, alpha, szPatch, szOverlap,  merge)

inputImg_rgb = inputImg;
tarImg_rgb = tarImg;
if ndims(inputImg_rgb)==2,
inputImg_rgb = repmat(inputImg_rgb,[1 1 3]);
inputImg = repmat(inputImg,[1 1 3]);
end; 
if ndims(tarImg_rgb)==2,
tarImg_rgb = repmat(tarImg_rgb,[1 1 3]);
tarImg = repmat(tarImg,[1 1 3]);
end; 
tarImg = rgb2gray(tarImg);
inputImg = rgb2gray(inputImg);
sizeout = size(tarImg);
outputImg = zeros(sizeout); 
outRGB = zeros(size(tarImg_rgb));
sizein = size(inputImg);
Y = zeros(size(outRGB));
M = ones(szPatch,szPatch);

%% Main Function

    % prepare top error to calculate SSD
    temp = ones([szOverlap szPatch]);
    errTop = xcorr2(inputImg.^2, temp);
    temp = ones([szPatch szOverlap]);
    errSide = xcorr2(inputImg.^2, temp);
    temp = ones([szPatch-szOverlap szOverlap]);
    errSidesmall = xcorr2(inputImg.^2, temp); 
    temp = ones([szPatch szPatch]);
    errorTarget = xcorr2(inputImg.^2, temp);
    
  

    
for iter = 1 : 1,  
    
   
    for i = 1:szPatch-szOverlap:sizeout(1)-szPatch+1
        
        for j = 1:szPatch-szOverlap:sizeout(2)-szPatch+1
            
%          startI = (i2-1)*(szPatch - szOverlap) + 1;
%          startJ = (j2-1)*(szPatch - szOverlap) + 1;
%          endI = startI + szPatch -1 ;
%          endJ = startJ + szPatch -1;
        
        if (i > 1) && (j > 1), 
            sharedTop = outputImg(i:i+szOverlap-1,j:j+szPatch-1);
            err = errTop - 2 * xcorr2(inputImg, sharedTop) + sum(sharedTop(:).^2);
            err = err(szOverlap:end-szPatch+1,szPatch:end-szPatch+1); 
            sharedSide = outputImg(i+szOverlap:i+szPatch-1,j:j+szOverlap-1);
            err2 = errSidesmall - 2 * xcorr2(inputImg, sharedSide) + sum(sharedSide(:).^2);
            err2=err2(szPatch:end-szPatch+szOverlap+1, szOverlap:end-szPatch+1);
            err = err + err2;
            tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);
            errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); 
            if(iter >1)
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);
                err = alpha*(err + errSyn) + (1-alpha) * errTarget;
            else
            err = err + alpha * errTarget; 
            end
            [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
            c = ceil(rand * length(ibest));
            pos = [ibest(c) jbest(c)]; 
            
            
           
        elseif i > 1 % left part
            sharedTop = outputImg(i:i+szOverlap-1,j:j+szPatch-1);
            err = errTop - 2 * xcorr2(inputImg, sharedTop) + sum(sharedTop(:).^2);
            err = err(szOverlap:end-szPatch+1,szPatch:end-szPatch+1); 
            tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);
            errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); 
            if(iter >1)
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);
                err = alpha*(err + errSyn) + (1-alpha) * errTarget;
            else
            err =  err + alpha * errTarget; 
            end
            [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
            c = ceil(rand * length(ibest));
            pos = [ibest(c) jbest(c)];  
            imout(i:i+tilesize-1,j:j+tilesize-1,:) = imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+tilesize-1,:);
      
      %min cut
      B1overlaph = imout(i:i+overlap-1,j:j+tilesize-1,:); % shared
      B2overlaph = imin(pos(1):pos(1)+overlap-1,pos(2):pos(2)+tilesize-1,:);
      
      errmat = sum((B1overlaph-B2overlaph).^2,3);
      
      fph = zeros(overlap,tilesize);
      pathh = zeros(overlap,tilesize);
      
      fph(:,tilesize) = errmat(:,tilesize);
      
      for k = tilesize-1:-1:1
          for l = 1:overlap
              index =  max(1,l-1):min(overlap,l+1);
              [fph(l,k), temp_index] = min( fph(index,k+1));
              fph(l,k) = fph(l,k) + errmat(l,k);
              pathh(l,k) = index(temp_index);
          end
      end
      
      min_err_boundh = zeros(1,tilesize);
      
      [temp,min_err_boundh(1)] = min(fph(1,:));
      
      for k=2:tilesize
          min_err_boundh(k) = pathh(min_err_boundh(k-1),k-1);
      end
      
      for k = 1:tilesize
          imout(i+min_err_boundh(k):i+tilesize-1,j+k-1,:) = imin(pos(1)+min_err_boundh(k):pos(1)+tilesize-1,pos(2)+k-1,:);
      end
        elseif j > 1 % top part
            sharedSide = outputImg(i:i+szPatch-1,j:j+szOverlap-1);
            err = errSide - 2 * xcorr2(inputImg, sharedSide) + sum(sharedSide(:).^2);
            err = err(szPatch:end-szPatch+1,szOverlap:end-szPatch+1);
            tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);
            errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1);         
            if(iter >1)
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);
                err = alpha*(err + errSyn) + (1-alpha) * errTarget;
            else
            err = err + alpha* errTarget;            
            end
            [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
            c = ceil(rand * length(ibest));
            pos = [ibest(c) jbest(c)];   
            
      
      %imout(i:i+tilesize-1,j:j+tilesize-1) = imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+tilesize-1);
    
      pos = ceil(rand([1 2]) .* (sizein-tilesize+1));
      imout(i:i+tilesize-1,j:j+tilesize-1,:) = imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+tilesize-1,:);
    
        else 
            if(iter >1)
                tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
                errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); 
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);
                err = errSyn + alpha * errTarget;
                [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
                c = ceil(rand * length(ibest));
                pos = [ibest(c) jbest(c)]; 
             else
                pos = ceil(rand([1 2]) .* (sizein-szPatch+1));
            end
               
      
        end
        %min cut
%         if(i>1)
%             E = ( inputImg(pos(1):pos(1)+szOverlap-1, pos(2):pos(2)+szPatch-1) - Y(i:i+szOverlap-1, j:j+szPatch-1) ).^2;
%             C = mincut(E, 1);
%             M(1:szOverlap, 1:end) = M(1:szOverlap, 1:end) .* double(C >= 0);
%         end
%         if(j>1)
%             E = ( inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szOverlap-1) - Y(i:i+szPatch-1, j:j+szOverlap-1) ).^2;
%             C = mincut(E, 0);
%             M(1:end, 1:szOverlap) = double(C >= 0);
%         end
%         if( i == 1 && j == 1 )
%                 Y(i:i+szPatch-1, j:j+szPatch-1, 1:3) = inputImg_rgb(pos(1):pos(1)+szPatch-1, pos(2):pos(2)+szPatch-1, 1:3);
%         else
%                 Y(i:i+szPatch-1, j:j+szPatch-1, :) = filtered_write(Y(i:i+szPatch-1, j:j+szPatch-1, :),inputImg_rgb(pos(1):pos(1)+szPatch-1, pos(2):pos(2)+szPatch-1, :), M);
%         end
%         blend_mask_rgb = repmat(M,[1 1 3]);
        outRGB(i:i+szPatch-1,j:j+szPatch-1,:) = inputImg_rgb(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1,:);
        end
i   
    end
end
end
function A = filtered_write(A, B, M)
for i = 1:3
    A(:, :, i) = A(:,:,i) .* (M == 0) + B(:,:,i) .* (M == 1);
end
end
function C = mincut(X,dir)

if( nargin > 1 && dir == 1 )
    X = X';
end

%Allocate the current cost array, and set first row to first row of X
E = zeros(size(X));
E(1:end,:) = X(1:end,:);

%Starting with the second array, compute the path costs until the end
for i=2:size(E,1),
    E(i,1) = X(i,1) + min( E(i-1,1), E(i-1,2) );
    for j=2:size(E,2)-1,
        E(i,j) = X(i,j) + min( [E(i-1,j-1), E(i-1,j), E(i-1,j+1)] );
    end
    E(i,end) = X(i,end) + min( E(i-1,end-1), E(i-1,end) );

end

%Backtrace to find the cut
C = zeros(size(X));

[cost, idx] = min(E(end, 1:end));
C(i, 1:idx-1) = -1;
C(i, idx) = 0;
C(i, idx+1:end) = +1;

for i=size(E,1)-1:-1:1,
    for j=1:size(E,2),

        if( idx > 1 && E(i,idx-1) == min(E(i,idx-1:min(idx+1,size(E,2))) ) )
            idx = idx-1;
        elseif( idx < size(E,2) && E(i,idx+1) == min(E(i,max(idx-1,1):idx+1)) )
            idx = idx+1;
        end


        C(i, 1:idx-1) = -1;
        C(i, idx) = 0;
        C(i, idx+1:end) = +1;

    end
end

if( nargin > 1 && dir == 1 )
    %E = E';
    C = C';
end
end
function A = filtered_write(A, B, M)
for i = 1:3
    A(:, :, i) = A(:,:,i) .* (M == 0) + B(:,:,i) .* (M == 1);
end
end