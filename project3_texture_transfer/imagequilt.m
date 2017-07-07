%function imout = imagequilt(imim, blocksize, n, overlap, err)
%Performs the Efros/Freeman Image quilting algorithm on the input

%%Inputs:
%%imin: source image
%%blocksize: the size of block defined by user
%%n: the number of blocks used in total; used to calculate destination image
%%size by having destination size = n * blocksize - (n-1) * overlap
%%overlap: the amount of overlap between tiles (by default 1/6 of the tile size)
%%err: an error tolerance (by default 0.1)


function imout = imagequilt(imin, blocksize, n, overlap, err)

imin = double(imin);
destinationsize = n * blocksize - (n-1) * overlap;


if( length(size(imin)) == 2 )
    imin = repmat(imin, [1 1 3]);
elseif( length(size(imin)) ~= 3 )
    error('Input image must be 2 or 3 dimensional');
end;
    
simple = 0;

if( nargin < 5 )
    err = 0.002;
end;

if( nargin < 4 )
    overlap = round(blocksize / 6);
end;

if( overlap >= blocksize )
    error('Overlap must be less than tilesize');
end;

 

imout = zeros(destinationsize, destinationsize, 3);

for i=1:n,
     for j=1:n,
         startI = (i-1)*blocksize -(i-1)*overlap + 1;
         startJ = (j-1)*blocksize -(j-1)*overlap + 1;
         endI = startI + blocksize -1 ;
         endJ = startJ + blocksize -1;
         
         %Determine the distances from each tile to the overlap region
         %This will eventually be replaced with convolutions
         distances = zeros( size(imin,1)-blocksize, size(imin,2)-blocksize );         

%          K = ones(tilesize, overlap);
%          y = Y(startI:endI,startJ:endJ,1:3);
%          for k=1:3,
%                 
%          end;
         
         
%          M = zeros(tilesize, tilesize);
%          M(1:overlap,:) = 1;
%          M(:,1:overlap) = 1;
%          
%          y = Y(startI:endI,startJ:endJ,1:3);
%          a = (y(:,:,1) .* M) + (y(:,:,2) .* M) + (y(:,:,3) .* M);
%          a2 = sum(sum(sum(a.^2)));
%             
%          %a2 = sum(sum(sum(a.^2)));
%          
%          b2 = filter2(M, b.^2);
%          
%          ab = filter2(a, b);
%      
%          distances = sqrt(a2 + (2*ab + b2));
%          distances = distances(1:size(X,1)-tilesize, 1:size(X,2)-tilesize);
         
%          b2 = filter2(X, zerosones(tilesize, tilesize)
%          
%          a2 = sum(sum(Y(startI:endI, startJ:startJ+overlap-1, 1:3).^2)) + ...
%               sum(sum(Y(startI:startI+overlap-1,startJ+overlap:endJ, 1:3).^2));
%          b2 = sum(sum(
%          

        useconv = 1;
        
        if( useconv == 0 )
            
            %Compute the distances from the template to target for all i,j
            for a = 1:size(distances,1)
                v1 = imout(startI:startI+overlap-1, startJ:startJ+overlap-1, 1:3);
                for b = 1:size(distances,2),                 
                    v2 = imin(a:a+overlap-1,b:b+overlap-1, 1:3);
                    distances(a,b) = myssd( double((v1(:) > 0)) .* (v1(:) - v2(:)) );
                    %distances(a,b) = D;    
                end;
            end;
            
        else
            
            %Compute the distances from the source to the left overlap region
            if( j > 1 )
                distances = ssd( imin, imout(startI:endI, startJ:startJ+overlap-1, 1:3) );    
                distances = distances(1:end, 1:end-blocksize+overlap);
            end;
            
            %Compute the distance from the source to top overlap region
            if( i > 1 )
                Z = ssd( imin, imout(startI:startI+overlap-1, startJ:endJ, 1:3) );
                Z = Z(1:end-blocksize+overlap, 1:end);
                if( j > 1 ) distances = distances + Z;
                else distances = Z;
                end;
            end;
            
            %If both are greater, compute the distance of the overlap
            if( i > 1 && j > 1 )
                Z = ssd( imin, imout(startI:startI+overlap-1, startJ:startJ+overlap-1, 1:3) );
                Z = Z(1:end-blocksize+overlap, 1:end-blocksize+overlap);                   
                distances = distances - Z;
            end;
            
            %distances = distances(1:end-tilesize, 1:end-tilesize);
            
        end;

         %Find the best candidates for the match
         best = min(distances(:));
         candidates = find(distances(:) <= (1+err)*best);
          
         idx = candidates(ceil(rand(1)*length(candidates)));
                         
         [sub(1), sub(2)] = ind2sub(size(distances), idx);
         fprintf( 'Picked tile (%d, %d) out of %d candidates.  Best error=%.4f\n', sub(1), sub(2), length(candidates), best );       
         
         %If we do the simple quilting (no cut), just copy image
         if( simple )
             imout(startI:endI, startJ:endJ, 1:3) = imin(sub(1):sub(1)+blocksize-1, sub(2):sub(2)+blocksize-1, 1:3);
         else
             
             %Initialize the mask to all ones
             M = ones(blocksize, blocksize);
             
             %We have a left overlap
             if( j > 1 )
                 
                 %Compute the SSD in the border region
                 E = ( imin(sub(1):sub(1)+blocksize-1, sub(2):sub(2)+overlap-1) - imout(startI:endI, startJ:startJ+overlap-1) ).^2;
                 
                 %Compute the mincut array
                 C = mincut(E, 0);
                 
                 %Compute the mask and write to the destination
                 M(1:end, 1:overlap) = double(C >= 0);
                 %Y(startI:endI, startJ:endJ, :) = filtered_write(Y(startI:endI, startJ:endJ, :), ...
                 %    X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :), M); 
                 
                 %Y(startI:endI, startJ:endJ, 1:3) = X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, 1:3);
                 
                 %Compute the mask and write to the destination
                 %                  M = zeros(tilesize, tilesize);
                 %                  M(1:end, 1:overlap) = double(C == 0);
                 %                  Y(startI:endI, startJ:endJ, :) = filtered_write(Y(startI:endI, startJ:endJ, :), ...
                 %                      repmat(255, [tilesize, tilesize, 3]), M); 
                 
             end;
             
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
             
         end;
                 
             
         image(uint8(imout));
         drawnow;
     end;
end;

%figure;
%image(uint8(imout));

function y = myssd( x )
y = sum( x.^2 );

function A = filtered_write(A, B, M)
for i = 1:3,
    A(:, :, i) = A(:,:,i) .* (M == 0) + B(:,:,i) .* (M == 1);
end;
