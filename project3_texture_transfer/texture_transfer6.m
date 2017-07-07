function [outputImg outRGB] = texture_transfer6(inputImg, tarImg, alpha, szPatch, szOverlap, niter, isdebug, merge)


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
filt_w = 4;
smooth_filt = binomialFilter(filt_w)*binomialFilter(filt_w)'; 

for iter = 1 : 2,
    temp = ones([szOverlap szPatch]);
    errTop = xcorr2(inputImg.^2, temp);
    temp = ones([szPatch szOverlap]);
    errSide = xcorr2(inputImg.^2, temp);
    temp = ones([szPatch-szOverlap szOverlap]);
    errSidesmall = xcorr2(inputImg.^2, temp); 
    temp = ones([szPatch szPatch]);
    errorTarget = xcorr2(inputImg.^2, temp);
    for i = [1:szPatch-szOverlap:sizeout(1)-szPatch+1, sizeout(1)-szPatch+1]
        for j = [1:szPatch-szOverlap:sizeout(2)-szPatch+1,sizeout(2)-szPatch+1]
       
        blend_mask = logical(zeros([szPatch szPatch]));

        if (i > 1) & (j > 1), 
            sharedTop = outputImg(i:i+szOverlap-1,j:j+szPatch-1);
            err = errTop - 2 * xcorr2(inputImg, sharedTop) + sum(sharedTop(:).^2);
            err = err(szOverlap:end-szPatch+1,szPatch:end-szPatch+1); 

            sharedSide = outputImg(i+szOverlap:i+szPatch-1,j:j+szOverlap-1);
            err2 = errSidesmall - 2 * xcorr2(inputImg, sharedSide) + sum(sharedSide(:).^2);
            err = err + err2(szPatch:end-szPatch+szOverlap+1, szOverlap:end-szPatch+1);

            tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);

            errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); 

            if(iter >1)
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);

                err = err + errSyn + alpha * errTarget;
            else
                err = err + alpha * errTarget; 
            end

            [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
            c = ceil(rand * length(ibest));
            pos = [ibest(c) jbest(c)]; 

            previous = outputImg(i:i+szPatch-1,j:j+szPatch-1);
            newPatch = inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1);

            err_sq = (newPatch - previous).^2;

            blend_mask = logical(zeros(size(err_sq))); 
            blend_mask = dpmain(err_sq,szOverlap); 


        elseif i > 1 
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
                err = err + errSyn + alpha * errTarget;
            else
                err = err + alpha * errTarget; 
            end
            [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
            c = ceil(rand * length(ibest));
            pos = [ibest(c) jbest(c)]; 
            previous = outputImg(i:i+szPatch-1,j:j+szPatch-1);
            newPatch = inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1);
            err_sq = (newPatch - previous).^2;

            blend_mask = logical(zeros(size(err_sq))); 
            blend_mask(1:szOverlap,:) = dp(err_sq(1:szOverlap,:)')';

            
        elseif j > 1 
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

                err = err + errSyn + alpha * errTarget;
            else
                err = err + alpha * errTarget; 
            end
            [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
            c = ceil(rand * length(ibest));
            pos = [ibest(c) jbest(c)]; 
            previous = outputImg(i:i+szPatch-1,j:j+szPatch-1);
            newPatch = inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1);
            err_sq = (newPatch - previous).^2;
            blend_mask = logical(zeros(size(err_sq)));
            blend_mask(:,1:szOverlap) = dp(err_sq(:,1:szOverlap));


            
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

          blend_mask = rconv2(double(blend_mask),smooth_filt);
          blend_mask_rgb = repmat(blend_mask,[1 1 3]);

        outputTemp = outputImg(i:i+szPatch-1,j:j+szPatch-1).* blend_mask + inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1).*(1- blend_mask);   
        outputImg(i:i+szPatch-1,j:j+szPatch-1) = outputTemp;
        outputTemp = outRGB(i:i+szPatch-1,j:j+szPatch-1,:).* blend_mask_rgb + inputImg_rgb(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1,:).*(1- blend_mask_rgb);
        outRGB(i:i+szPatch-1,j:j+szPatch-1,:) = outputTemp;

    end
i   
end


end


end