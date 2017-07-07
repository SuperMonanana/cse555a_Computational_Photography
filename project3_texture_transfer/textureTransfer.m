function [outputImg outRGB] = textureTransfer(inputImg_rgb, tarImg_rgb, alpha, szPatch, szOverlap, isdebug)

%% Input and Config

tarImg = rgb2gray(tarImg_rgb);
inputImg = rgb2gray(inputImg_rgb);


% inputImg = im2double(inputImg);
% tarImg = im2double(tarImg);

sizeout = size(tarImg);

outputImg = zeros(size(tarImg));
sizein = size(inputImg);

outputImg = zeros(sizeout); % output image
outRGB = zeros(size(tarImg_rgb));

% if(size(inputImg,3) ==1)
%     outputImg = zeros(sizeout);
% else
%     outputImg = zeros([sizeout(1:2) size(inputImg,3)]);
% end
% outputImg(1,1,:)=255;


if nargin<5
    isdebug = 0;
end

if isdebug~=0
    h = imshow(uint8(outputImg));
else
    isdebug = floor(isdebug);
end



rng(3290)
for iter = 1 : 1,
temp = ones([szOverlap szPatch]);
    errTop = xcorr2(inputImg.^2, temp);
    temp = ones([szPatch szOverlap]);
    errSide = xcorr2(inputImg.^2, temp);
    temp = ones([szPatch-szOverlap szOverlap]);
    errSidesmall = xcorr2(inputImg.^2, temp); % left bottom rigion

    temp = ones([szPatch szPatch]);
    errorTarget = xcorr2(inputImg.^2, temp);

%% Main Function
for i=1:szPatch-szOverlap:sizeout(1)-szPatch+1,
  for j=1:szPatch-szOverlap:sizeout(2)-szPatch+1,
    %% Synthesis texture given the existing top and left image
    if (i > 1) & (j > 1)
        % TODO#1 Find the existing shared region and target patch
        % Extract shared region and target region
        sharedTop = outputImg(i:i+szOverlap-1,j:j+szPatch-1);
            err = errTop - 2 * xcorr2(inputImg, sharedTop) + sum(sharedTop(:).^2);
            err = err(szOverlap:end-szPatch+1,szPatch:end-szPatch+1); % remove invalid around edge

            sharedSide = outputImg(i+szOverlap:i+szPatch-1,j:j+szOverlap-1);
            err2 = errSidesmall - 2 * xcorr2(inputImg, sharedSide) + sum(sharedSide(:).^2);
            % trim the edge region
            err = err + err2(szPatch:end-szPatch+szOverlap+1, szOverlap:end-szPatch+1);



            % correspond error
            tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);

            errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); % position matched

        % % Total error
        if(iter >1)
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);

                err = err + errSyn + alpha * errTarget;
            else
                err = err + alpha * errTarget; 
            end


        % Find the texture with a reasonable small error
        [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
        c = ceil(rand * length(ibest));
        pos = [ibest(c) jbest(c)];

    elseif i > 1
        % TODO#5 Find the existing shared region and target patch
        % Extract shared region and target region
       % only sharedTop
            sharedTop = outputImg(i:i+szOverlap-1,j:j+szPatch-1);
            err = errTop - 2 * xcorr2(inputImg, sharedTop) + sum(sharedTop(:).^2);
            err = err(szOverlap:end-szPatch+1,szPatch:end-szPatch+1); % remove invalid around edge

            % correspond error
            tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);

            errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); % position matched

        % % Total error
        if(iter >1)
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);

                err = err + errSyn + alpha * errTarget;
        else
                err = err + alpha * errTarget; 
        end


        % Find the texture with a reasonable small error
        [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
        c = ceil(rand * length(ibest));
        pos = [ibest(c) jbest(c)];
    elseif j > 1
        % TODO#9 Find the existing shared region and target patch
        % Extract shared region and target region
        sharedSide = outputImg(i:i+szPatch-1,j:j+szOverlap-1);
            err = errSide - 2 * getxcorr2(inputImg, sharedSide) + sum(sharedSide(:).^2);
            err = err(szPatch:end-szPatch+1,szOverlap:end-szPatch+1);

            % correspond error
            tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);

            errTarget = errorTarget - 2*getxcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); % position matched


        % % Total error
        if(iter >1)
                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);

                err = err + errSyn + alpha * errTarget;
            else
                err = err + alpha * errTarget; 
            end
    
        % Find the texture with a reasonable small error
        [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
        c = ceil(rand * length(ibest));
        pos = [ibest(c) jbest(c)];
    else
        if(iter >1)

                % correspond error
                tarPatch = tarImg(i:i+szPatch-1,j:j+szPatch-1,:);

                errTarget = errorTarget - 2*xcorr2(inputImg,tarPatch) + sum(tarPatch(:).^2);
                errTarget = errTarget(szPatch:end-szPatch+1, szPatch:end-szPatch+1); % position matched

                synthezised_patch = outputImg(i:i+szPatch-1,j:j+szPatch-1,:);
                errSyn = errorTarget - 2*xcorr2(inputImg,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(szPatch:end-szPatch+1, szPatch:end-szPatch+1);

                err = errSyn + alpha * errTarget;

                 % sample a patach from candidates
                 [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
                 c = ceil(rand * length(ibest));
                 pos = [ibest(c) jbest(c)]; 
             else
                pos = ceil(rand([1 2]) .* (sizein-szPatch+1));
            end
    end
     % find minimum error boundary cut 
%             previous = outputImg(i:i+szPatch-1,j:j+szPatch-1);
%             newPatch = inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1);
% 
%             err_sq = (newPatch - previous).^2;
% 
%             blend_mask = logical(zeros(size(err_sq)));  % to record the border
%             blend_mask = dpmain(err_sq,szOverlap); 
    

%             blend_mask = logical(zeros([szPatch szPatch]));
%         
%             
% 
% 
%             blend_mask_rgb = repmat(blend_mask,[1 1 3]);
            
             %outputImg(i:i+szPatch-1,j:j+szPatch-1) = inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1);
%   outputTemp = outputImg(i:i+szPatch-1,j:j+szPatch-1).* blend_mask + inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1).*(1- blend_mask);
        %outputTemp = inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1);
        outputImg(i:i+szPatch-1,j:j+szPatch-1)=inputImg(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1);
        


        % RGB version here
%         outputTemp = outRGB(i:i+szPatch-1,j:j+szPatch-1,:).* blend_mask_rgb + inputImg_rgb(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1,:).*(1- blend_mask_rgb);
        %outputTemp = inputImg_rgb(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1,:);
        outRGB(i:i+szPatch-1,j:j+szPatch-1,:) = inputImg_rgb(pos(1):pos(1)+szPatch-1,pos(2):pos(2)+szPatch-1,:);
    %% Show debug result
    if isdebug~=0
%         figure(3), imshow(outputImg, []);
figure(3);
            subplot(1,2,1);
            imshow(outputImg);
            subplot(1,2,2);
            imshow(outRGB);
    end
    
  end
  i
end
end