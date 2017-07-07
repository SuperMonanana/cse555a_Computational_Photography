clear all;close all;
%%READ DATA
I=imread('01087v.jpg');

%DIVIDE IT INTO 3 EQUAL PARTS
[r c]=size(I);
IB=I(1:r/3,1:c);
IG=I(r/3+1:2*r/3,1:c);
IR=I(2*r/3+1:r,1:c);

VarG=sum(sum( (IB-IG).^2 ));
IGmove=IG;
for i=-20:20
    for j=-20:20
        seG=translate(strel(1),[i,j]);
        IGTrymove=imdilate(IG,seG);
        VarGMOVE=sum(sum( (IB-IGTrymove).^2 ));
        if VarGMOVE<VarG
            VarG=VarGMOVE;
            IGmove=IGTrymove;
            displacementG=[i j];
        end
            
    end
end
IG=IGmove;

VarR=sum(sum( (IB-IR).^2 ));
IRmove=IR;
for a=-20:20
    for b=-20:20
        seR=translate(strel(1),[a,b]);
        IRTrymove=imdilate(IR,seR);
        VarRMOVE=sum(sum( (IB-IRTrymove).^2 ));
        if VarRMOVE<VarR
            VarR=VarRMOVE;
            IRmove=IRTrymove;
            displacementR=[a b];
        end
            
    end
end
IR=IRmove;

figure(2)
imshow(cat(3,IR,IG,IB));


