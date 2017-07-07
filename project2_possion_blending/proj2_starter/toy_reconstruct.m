function [im_out,v] = toy_reconstruct(im)
 
[imh, imw, nb] = size(im); 
im2var = zeros(imh, imw); 
im2var(1:imh*imw) = 1:imh*imw; 
b= [];
A=sparse([]);
v=zeros(1,imh*imw);


%%objective1

    e=0;
    for x=1:imw-1
        for y=1:imh 
            e=e+1
            A(e, im2var(y,x+1))=1; 
            A(e, im2var(y,x))=-1; 
            b(e) = im(y,x+1)-im(y,x); 
            
        end
    end
    
%%objective2
    
    for x=1:imw
        for y=1:imh-1
            e=e+1
            A(e, im2var(y+1,x))=1; 
            A(e, im2var(y,x))=-1; 
            b(e) = im(y+1,x)-im(y,x);           
        end
    end

%%objective3
    e=e+1;
    A(e, im2var(1,1))=1; 
    b(e)=im(1,1); 
        

v = A\b';


im_out = zeros(size(im));
i = 1;
for x=1:imw
    for y=1:imh
        im_out(y,x) = v(im2var(y,x));
    end
end
