function imout = basicquilt(imin, blocksize, n, overlap, err)
%%Implement the basic algotithm with SSD
%%Inputs:
%%imin: source image
%%blocksize: the size of block defined by user
%%n: the number of blocks used in total; used to calculate destination image
%%size by having destination size = n * blocksize - (n-1) * overlap
%%overlap: the amount of overlap between tiles (by default 1/6 of the tile size)
%%err: an error tolerance (by default 0.1)

%%Initialization
imin = double(imin);
destinationsize = n * blocksize - (n-1) * overlap;
imout = zeros(destinationsize, destinationsize, 3);
simple=0;

for i=1:n,
     for j=1:n,
%%Seed the output image by copying a random tile from the source image to the top-left corner of the destination.
        if i==1&j==1
            sizein = size(imin);
            seed = ceil(rand([1 2]) .* (sizein(1:2)-blocksize+1));
            imout(i:i+blocksize-1,j:j+blocksize-1,:) = imin(seed(1):seed(1)+blocksize-1,seed(2):seed(2)+blocksize-1,:);
        else
            startI = (i-1)*blocksize - (i-1) * overlap + 1;
            startJ = (j-1)*blocksize - (j-1) * overlap + 1;
            endI = startI + blocksize -1 ;
            endJ = startJ + blocksize -1;
            
            %Consider the error as "distance" in the graph problem
            %exhaustive search over every possible patch
            distances = zeros( size(imin,1)-blocksize, size(imin,2)-blocksize );
            for a = 1:size(distances,1)
                v1 = imout(startI:startI+overlap-1, startJ:startJ+overlap-1, 1:3);
                for b = 1:size(distances,2),                 
                    v2 = imin(a:a+overlap-1,b:b+overlap-1, 1:3);
                    distances(a,b) = sum((double((v1(:) > 0)).*(v1(:)-v2(:)) ).^2);
                       
                end
            end
            
            %Find the best candidates for the match
            best = min(distances(:));
            candidates = find(distances(:) <= (1+err)*best);
          
            idx = candidates(ceil(rand(1)*length(candidates)));
                         
            [sub(1), sub(2)] = ind2sub(size(distances), idx);
        
if simple==0
            imout(startI:endI, startJ:endJ, 1:3) = imin(sub(1):sub(1)+blocksize-1, sub(2):sub(2)+blocksize-1, 1:3);
else
end           %Initialize the mask to all ones
             M = ones(blocksize, blocksize);
             
             %We have a left overlap
             if( j > 1 )
                 
                 %Compute the SSD in the border region
                 E = ( imin(sub(1):sub(1)+blocksize-1, sub(2):sub(2)+overlap-1) - imout(startI:endI, startJ:startJ+overlap-1) ).^2;
                 
                 %Compute the mincut array
                 C = mincut(E, 0);
                 
                 %Compute the mask and write to the destination
                 M(1:end, 1:overlap) = double(C >= 0);
                
                 
             end
             
             %We have a top overlap
             if( i > 1 )
                 %Compute the SSD in the border region
                 E = ( imin(sub(1):sub(1)+overlap-1, sub(2):sub(2)+blocksize-1) - imout(startI:startI+overlap-1, startJ:endJ) ).^2;
                 
                 %Compute the mincut array
                 C = mincut(E, 1);
                 
                 %Compute the mask and write to the destination
                 M(1:overlap, 1:end) = M(1:overlap, 1:end) .* double(C >= 0);
                 %Y(startI:endI, startJ:endJ, :) = filtered_write(Y(startI:endI, startJ:endJ, :), ...
                 %    X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :), M); 
             end;
             
             
             if( i == 1 && j == 1 )
                 imout(startI:endI, startJ:endJ, 1:3) = imin(sub(1):sub(1)+blocksize-1, sub(2):sub(2)+blocksize-1, 1:3);
             else
                 %Write to the destination using the mask
                 imout(startI:endI, startJ:endJ, :) = filtered_write(imout(startI:endI, startJ:endJ, :), ...
                     imin(sub(1):sub(1)+blocksize-1, sub(2):sub(2)+blocksize-1, :), M); 
             end;
             

            
            
            
            
            image(uint8(imout));
             drawnow;
        end
        end
     end

end

function A = filtered_write(A, B, M)
for i = 1:3,
    A(:, :, i) = A(:,:,i) .* (M == 0) + B(:,:,i) .* (M == 1);
end
end
