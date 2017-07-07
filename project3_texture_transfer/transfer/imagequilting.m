function [imout] = imagequilting(imin, tarImg,sizeout, blocksize, overlap)



imin = double(imin);
tarImg=double(tarImg);

if(size(imin,3) ==1)
    imout = zeros(sizeout);
else
    imout = zeros([sizeout(1:2) size(imin,3)]);
end
imout(1,1,:)=255;

sizein = size(imin);
sizein = sizein(1:2);

for iter=1:5
    alpha = 0.8*(iter-1)/4+0.1;
    temp = ones([overlap blocksize]);
    errtop = getxcorr2(imin.^2, temp);
    temp = ones([blocksize overlap]);
    errside = getxcorr2(imin.^2, temp);
    temp = ones([blocksize-overlap overlap]);
    errsidesmall = getxcorr2(imin.^2, temp);
    %for transfer
    temp = ones([blocksize blocksize]);
    errorTarget = getxcorr2(imin.^2, temp);




for i=1:blocksize-overlap:sizeout(1)-blocksize+1,
  for j=1:blocksize-overlap:sizeout(2)-blocksize+1,

    if (i > 1) & (j > 1),
    % for top 
      shared = imout(i:i+overlap-1,j:j+blocksize-1,:);
      err = errtop - 2 * getxcorr2(imin, shared) + sum(shared(:).^2);
      err = err(overlap:end-blocksize+1,blocksize:end-blocksize+1);

      % for left
      shared = imout(i+overlap:i+blocksize-1,j:j+overlap-1,:);
      err2 = errsidesmall - 2 * getxcorr2(imin, shared) + sum(shared(:).^2);
     
      err = err + err2(blocksize:end-blocksize+overlap+1, overlap:end-blocksize+1);

      %%target picture
            tarPatch = tarImg(i:i+blocksize-1,j:j+blocksize-1,:);
            errTarget = errorTarget - 2*getxcorr2(imin,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(blocksize:end-blocksize+1, blocksize:end-blocksize+1); 
            if(iter >1)
                synthezised_patch = imout(i:i+blocksize-1,j:j+blocksize-1,:);
                errSyn = errorTarget - 2*getxcorr2(imin,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(blocksize:end-blocksize+1, blocksize:end-blocksize+1);
                err = alpha*(err + errSyn) + (1-alpha) * errTarget;
            else
            err = alpha*err + (1-alpha) * errTarget; 
            end
      
      [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
      c = ceil(rand * length(ibest));
      pos = [ibest(c) jbest(c)];
      
      %min cut
          B1overlaph = imout(i:i+overlap-1,j:j+blocksize-1,:); % shared
          B2overlaph = imin(pos(1):pos(1)+overlap-1,pos(2):pos(2)+blocksize-1,:);

          errmat = sum((B1overlaph-B2overlaph).^2,3);

          fph = zeros(overlap,blocksize);
          pathh = zeros(overlap,blocksize);

          fph(:,blocksize) = errmat(:,blocksize);

          for k = blocksize-1:-1:1
              for l = 1:overlap
                  index =  max(1,l-1):min(overlap,l+1);
                  [fph(l,k), temp_index] = min( fph(index,k+1));
                  fph(l,k) = fph(l,k) + errmat(l,k);
                  pathh(l,k) = index(temp_index);
              end
          end
          
          B1overlap = imout(i:i+blocksize-1,j:j+overlap-1,:); % shared
          B2overlap = imin(pos(1):pos(1)+blocksize-1,pos(2):pos(2)+overlap-1,:);

          errmat = sum((B1overlap-B2overlap).^2,3);

          fp = zeros(blocksize,overlap);
          path = zeros(blocksize,overlap);

          fp(blocksize,:) = errmat(blocksize,:);

          for k = blocksize-1:-1:1
              for l = 1:overlap
                  index =  max(1,l-1):min(overlap,l+1);
                  [fp(k,l), temp_index] = min( fp(k+1,index));
                  fp(k,l) = fp(k,l) + errmat(k,l);
                  path(k,l) = index(temp_index);
              end
          end
          
          allerr = fp(1:overlap,1:overlap) + fph(1:overlap,1:overlap);
          
          [tempmin,tempindminclom] = min(allerr);
          [temp, min_bound_indexj] = min(tempmin);
          min_bound_indexi = tempindminclom(min_bound_indexj);
          
          imout(i+ overlap : i+blocksize-1,j + overlap : j+ blocksize-1,:) = ...
              imin(pos(1)+overlap :pos(1)+blocksize-1,pos(2)+overlap :pos(2)+blocksize-1,:);
          
          
      
          %imout(i:i+tilesize-1,j:j+tilesize-1) = imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+tilesize-1);
          
          min_err_bound = zeros(1,blocksize);
          min_err_boundh = zeros(1,blocksize);
      
          min_err_bound(min_bound_indexi) = min_bound_indexj;
          min_err_boundh(min_bound_indexj) = min_bound_indexi;
      
          for k=min_bound_indexi+1 :blocksize
              min_err_bound(k) = path(k-1,min_err_bound(k-1));
          end
          
          for k=min_bound_indexj+1 :blocksize
              min_err_boundh(k) = pathh(min_err_boundh(k-1),k-1);
          end
          
          for k = overlap : blocksize
            imout(i+min_err_boundh(k):i+overlap-1,j+k-1,:) = imin(pos(1)+min_err_boundh(k):pos(1)+overlap-1,pos(2)+k-1,:);
          end
          
          for k = overlap:blocksize
            imout(i+k-1,j+min_err_bound(k):j+overlap-1,:) = imin(pos(1)+k-1,pos(2)+min_err_bound(k):pos(2)+overlap-1,:);
          end
          
          
          for k = min_bound_indexi:overlap-1
              for l = min_bound_indexj:overlap-1
                  if k>=min_err_boundh(l) && l>= min_err_bound(k)
                      imout(i+k,j+l,:) = imin(pos(1)+k,pos(2)+l,:);
                  end
              end
          end
           
    elseif i > 1
      shared = imout(i:i+overlap-1,j:j+blocksize-1,:);
      err = errtop - 2 * getxcorr2(imin, shared) + sum(shared(:).^2);
      
      err = err(overlap:end-blocksize+1,blocksize:end-blocksize+1,:);
     %target 
            tarPatch = tarImg(i:i+blocksize-1,j:j+blocksize-1,:);
            errTarget = errorTarget - 2*getxcorr2(imin,tarPatch) + sum(tarPatch(:).^2);
            errTarget = errTarget(blocksize:end-blocksize+1, blocksize:end-blocksize+1); 
           if(iter >1)
                synthezised_patch = imout(i:i+blocksize-1,j:j+blocksize-1,:);
                errSyn = errorTarget - 2*getxcorr2(imin,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(blocksize:end-blocksize+1, blocksize:end-blocksize+1);
                err = alpha*(err + errSyn) + (1-alpha) * errTarget;
           else 
            err =  alpha*err + (1-alpha) * errTarget; 
           end
      [ibest, jbest] = find(err <= 1.01*1.1*min(err(:)));
      c = ceil(rand * length(ibest));
      pos = [ibest(c) jbest(c)];
      imout(i:i+blocksize-1,j:j+blocksize-1,:) = imin(pos(1):pos(1)+blocksize-1,pos(2):pos(2)+blocksize-1,:);
      
      %min cut
      B1overlaph = imout(i:i+overlap-1,j:j+blocksize-1,:); % shared
      B2overlaph = imin(pos(1):pos(1)+overlap-1,pos(2):pos(2)+blocksize-1,:);
      
      errmat = sum((B1overlaph-B2overlaph).^2,3);
      
      fph = zeros(overlap,blocksize);
      pathh = zeros(overlap,blocksize);
      
      fph(:,blocksize) = errmat(:,blocksize);
      
      for k = blocksize-1:-1:1
          for l = 1:overlap
              index =  max(1,l-1):min(overlap,l+1);
              [fph(l,k), temp_index] = min( fph(index,k+1));
              fph(l,k) = fph(l,k) + errmat(l,k);
              pathh(l,k) = index(temp_index);
          end
      end
      
      min_err_boundh = zeros(1,blocksize);
      
      [temp,min_err_boundh(1)] = min(fph(1,:));
      
      for k=2:blocksize
          min_err_boundh(k) = pathh(min_err_boundh(k-1),k-1);
      end
      
      for k = 1:blocksize
          imout(i+min_err_boundh(k):i+blocksize-1,j+k-1,:) = imin(pos(1)+min_err_boundh(k):pos(1)+blocksize-1,pos(2)+k-1,:);
      end
      
      
      
      
    elseif j > 1
      shared = imout(i:i+blocksize-1,j:j+overlap-1,:);
      err = errside - 2 * getxcorr2(imin, shared) + sum(shared(:).^2);
      
      err = err(blocksize:end-blocksize+1,overlap:end-blocksize+1);
      
      %target
        tarPatch = tarImg(i:i+blocksize-1,j:j+blocksize-1,:);
        errTarget = errorTarget - 2*getxcorr2(imin,tarPatch) + sum(tarPatch(:).^2);
        errTarget = errTarget(blocksize:end-blocksize+1, blocksize:end-blocksize+1); 
            
        if(iter >1)
                synthezised_patch = imout(i:i+blocksize-1,j:j+blocksize-1,:);
                errSyn = errorTarget - 2*getxcorr2(imin,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(blocksize:end-blocksize+1, blocksize:end-blocksize+1);
                err = alpha*(err + errSyn) + (1-alpha) * errTarget;
        else
            err = alpha*err + (1-alpha)* errTarget;    
        end
        
      [ibest, jbest] = find(err <= 1.01*1.1*min(err(:)));
      c = ceil(rand * length(ibest));
      pos = [ibest(c) jbest(c)];
      
      %min cut
      B1overlap = imout(i:i+blocksize-1,j:j+overlap-1,:); % shared
      B2overlap = imin(pos(1):pos(1)+blocksize-1,pos(2):pos(2)+overlap-1,:);
      
      errmat = sum((B1overlap-B2overlap).^2,3);
      
      fp = zeros(blocksize,overlap);
      path = zeros(blocksize,overlap);
      
      fp(blocksize,:) = errmat(blocksize,:);
      
      for k = blocksize-1:-1:1
          for l = 1:overlap
              index =  max(1,l-1):min(overlap,l+1);
              [fp(k,l), temp_index] = min( fp(k+1,index));
              fp(k,l) = fp(k,l) + errmat(k,l);
              path(k,l) = index(temp_index);
          end
      end
      
      min_err_bound = zeros(1,blocksize);
      
      [temp,min_err_bound(1)] = min(fp(1,:));
      
      for k=2:blocksize
          min_err_bound(k) = path(k-1,min_err_bound(k-1));
      end
      
      for k = 1:blocksize
          imout(i+k-1,j+min_err_bound(k):j+blocksize-1,:) = imin(pos(1)+k-1,pos(2)+min_err_bound(k):pos(2)+blocksize-1,:);
      end
      
      %imout(i:i+tilesize-1,j:j+tilesize-1) = imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+tilesize-1);
    else
        if(iter >1)
                tarPatch = tarImg(i:i+blocksize-1,j:j+blocksize-1,:);
                errTarget = errorTarget - 2*getxcorr2(imin,tarPatch) + sum(tarPatch(:).^2);
                errTarget = errTarget(blocksize:end-blocksize+1, blocksize:end-blocksize+1); 
                synthezised_patch = imout(i:i+blocksize-1,j:j+blocksize-1,:);
                errSyn = errorTarget - 2*getxcorr2(imin,synthezised_patch) + sum(synthezised_patch(:).^2);
                errSyn = errSyn(blocksize:end-blocksize+1, blocksize:end-blocksize+1);
                err = alpha*errSyn + (1-alpha) * errTarget;
                [ibest, jbest] = find(err <= 1.1*1.01*min(err(:)));
                c = ceil(rand * length(ibest));
                pos = [ibest(c) jbest(c)]; 
        else
        pos = ceil(rand([1 2]) .* (sizein-blocksize+1));
        imout(i:i+blocksize-1,j:j+blocksize-1,:) = imin(pos(1):pos(1)+blocksize-1,pos(2):pos(2)+blocksize-1,:);
        end
      
    end


  end
  i
end
blocksize = round(blocksize*0.5);
overlap=round(overlap*0.5);

figure(iter)
imshow(imout)
end
