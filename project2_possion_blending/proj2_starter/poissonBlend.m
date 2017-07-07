function im_blend = poissonBlend(im_s, mask_s, im_background);

s=im_s;t=im_background;
[imbh, imbw, nb] = size(s); 
im2var = zeros(imbh, imbw); 
im2var(1:imbh*imbw) = 1:imbh*imbw; 


%initialization
b= zeros(imbh*imbw,3);
A=sparse(imbh*imbw,imbh*imbw);
v=zeros(1,imbh*imbw);
e=1;%which equation

for x=1:imbw
    for y=1:imbh
        if mask_s(y,x) == 0%outside the mask
            A(e, im2var(y,x)) = 1;
            for i=1:3
                b(e, i) = b(e,i)+t(y,x,i);
            end
        else
    

%         elseif mask_s(y,x) == 0&&(mask_s(y+1,x)==1||mask_s(y-1,x)==1||mask_s(y,x+1)==1||mask_s(y,x-1)==1)%inside the mask,not the boundary
            A(e, im2var(y,x)) = 4;
            A(e, im2var(y,x+1))= -1;
            A(e, im2var(y,x-1))= -1;
            A(e, im2var(y+1,x)) = -1;
            A(e, im2var(y-1,x)) = -1;
            
            for i = 1 : 3
                b(e, i) =b(e,i)+4 * s(y, x, i)...
                    - s(y + 1, x, i) ...
                    - s(y - 1, x, i) ...
                    - s(y, x + 1, i) ...
                    - s(y, x - 1, i);
            end
%         else
%             A(e, im2var(y,x)) = 1;
%             for i=1:3
%                 b(e, i) = b(e,i)+s(y,x,i);
%             end
%         
         end
        e=e+1;
    end
    x
end
        


%solve v        
vr = A\b(:,1);
vg = A\b(:,2);
vb = A\b(:,3);  

vr(vr < 0) = 0;
vg(vg < 0) = 0;
vb(vb < 0) = 0;


im_blend = zeros(imbh,imbw,3);

for x=1:imbw
    for y=1:imbh
        im_blend(y,x,1) = vr(im2var(y,x));
        im_blend(y,x,2) = vg(im2var(y,x));         
        im_blend(y,x,3) = vb(im2var(y,x));
        
    end
end
end
    