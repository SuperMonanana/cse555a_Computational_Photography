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
   
    for i = [1:szPatch-szOverlap:sizeout(1)-szPatch+1, sizeout(1)-szPatch+1]
        for j = [1:szPatch-szOverlap:sizeout(2)-szPatch+1,sizeout(2)-szPatch+1]
        startI = (i-1)*(szPatch - szOverlap) + 1;
        startJ = (j-1)*(szPatch - szOverlap) + 1;
        endI = startI + szPatch -1 ;
        endJ = startJ + szPatch -1;
        
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
                err = err + errSyn + (1-alpha) * errTarget;
            else
            err = alpha *err + (1-alpha) * errTarget; 
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
                err = err + errSyn + (1-alpha) * errTarget;
            else
            err =  alpha*err + (1-alpha) * errTarget; 
            end
            [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
            c = ceil(rand * length(ibest));
            pos = [ibest(c) jbest(c)];                         
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
                err = err + errSyn + (1-alpha) * errTarget;
            else
            err = alpha* err + (1-alpha)* errTarget;            
            end
            [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
            c = ceil(rand * length(ibest));
            pos = [ibest(c) jbest(c)];            
        else 
            if(iter >1)
                tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
                errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); 
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);
                err = errSyn + (1-alpha) * errTarget;
                [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
                c = ceil(rand * length(ibest));
                pos = [ibest(c) jbest(c)]; 
             else
                pos = ceil(rand([1 2]) .* (sizein-szPatch+1));
                end
        end
        if(i>1)
            E = ( inputImg(pos(1):pos(1)+szOverlap-1, pos(2):pos(2)+szPatch-1) - Y(startI:startI+szOverlap-1, startJ:endJ) ).^2;
            C = mincut(E, 1);
            M(1:szOverlap, 1:end) = M(1:szOverlap, 1:end) .* double(C >= 0);
        end
        if(j>1)
            E = ( inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szOverlap-1) - Y(startI:endI, startJ:startJ+szOverlap-1) ).^2;
            C = mincut(E, 0);
            M(1:end, 1:szOverlap) = double(C >= 0);
        end
        if( i == 1 && j == 1 )
                Y(startI:endI, startJ:endJ, 1:3) = inputImg_rgb(pos(1):pos(1)+szPatch-1, pos(2):pos(2)+szPatch-1, 1:3);
        else
                Y(startI:endI, startJ:endJ, :) = filtered_write(Y(startI:endI, startJ:endJ, :),inputImg_rgb(pos(1):pos(1)+szPatch-1, pos(2):pos(2)+szPatch-1, :), M);
        end
        blend_mask_rgb = repmat(M,[1 1 3]);
        outputTemp = outRGB(i:i+szPatch-1,j:j+szPatch-1,:).* blend_mask_rgb + inputImg_rgb(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1,:).*(1- blend_mask_rgb);
        outRGB(i:i+szPatch-1,j:j+szPatch-1,:) = outputTemp;
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