clear all;close all;
%%READ DATA
I=imread('00458u.tif');


%Construct the pyramid
[r c]=size(I);
Image=cell(5,1);

w=fspecial('gaussian',[3 3]);
Image{1}=I;
for t=2:5
Image{t}=imresize(imfilter(Image{t-1},w),[r/(2^(t-1)) c/(2^(t-1))]);
end


%DIVIDE IT INTO 3 EQUAL PARTS
range=-10:10;
displacementB=[0,0];
displacementG=[0,0];
displacementR=[0,0];
 for x=5:-1:1
   [m,n]=size(Image{x});
   IB=Image{x}(1:m/3,1:n);
   IG=Image{x}(m/3+1:2*m/3,1:n);
   IR=Image{x}(2*m/3+1:m,1:n);
   
%    VarB=sum(sum((IG-IB).^2 ))+sum(sum((IG-IB)'.^2));
%    i0=2*displacementB(1);j0=2*displacementB(2);
%    IBmove=IB;
%     for i=range
%         for j=range
%             seB=translate(strel(1),[i0+i,j0+j]);
%             IBTrymove=imdilate(IB,seB);
%             VarBMOVE=sum(sum ((IG-IBTrymove).^2))+sum(sum((IG-IBTrymove)'.^2));
%             if  VarBMOVE< VarB 
%                 VarB=VarBMOVE;
%                 IBmove=IBTrymove;
%                 displacementB=[i0+i j0+j];
%             end
%         end
%     end
%     IB=IBmove;
% 
%       VarR=sum(sum( (IG-IR).^2))+sum(((IG-IR)'.^2));
%      %VarR=max(dot( IB./abs(IB), IR./abs(IR)));
% 
% 
%     a0=2*displacementR(1);b0=2*displacementR(2);
%     IRmove=IR;
%     for a=range
%         for b=range
%             seR=translate(strel(1),[a0+a,b0+b]);
%             IRTrymove=imdilate(IR,seR);
%             VarRMOVE=sum(sum((IG-IRTrymove).^2))+sum(sum((IG-IRTrymove)'.^2));
% %             VarRMOVE=max(dot( IB./abs(IB), IRTrymove./abs(IRTrymove)));
% 
% 
% 
%             if VarRMOVE<VarR
%                 VarR=VarRMOVE;
%                 IRmove=IRTrymove;
%                 displacementR=[a0+a b0+b];
%             end
%             
%         end
%     end
%     IR=IRmove;
%     x
%  end
 
% VarG=sum(sum( (IB-IG).^2))+sum(((IB-IG)'.^2));
%    i0=2*displacementG(1);j0=2*displacementG(2);
%    IGmove=IG;
%     for i=range
%         for j=range
%             seG=translate(strel(1),[i0+i,j0+j]);
%             IGTrymove=imdilate(IG,seG);
%             VarGMOVE=sum(sum((IB-IGTrymove).^2))+sum(sum((IB-IGTrymove)'.^2));
%             if  VarGMOVE< VarG 
%                 VarG=VarGMOVE;
%                 IGmove=IGTrymove;
%                 displacementG=[i0+i j0+j];
%             end
%         end
%     end
%     IG=IGmove;
% 
%   VarR=sum(sum( (IB-IR).^2))+sum(((IB-IR)'.^2));
%      VarR=max(dot( IB./abs(IB), IR./abs(IR)));
% 
% 
%     a0=2*displacementR(1);b0=2*displacementR(2);
%     IRmove=IR;
%     for a=range
%         for b=range
%             seR=translate(strel(1),[a0+a,b0+b]);
%             IRTrymove=imdilate(IR,seR);
%             VarRMOVE=sum(sum((IB-IRTrymove).^2))+sum(sum((IB-IRTrymove)'.^2));
%             VarRMOVE=max(dot( IB./abs(IB), IRTrymove./abs(IRTrymove)));
% 
% 
% 
%             if VarRMOVE<VarR
%                 VarR=VarRMOVE;
%                 IRmove=IRTrymove;
%                 displacementR=[a0+a b0+b];
%             end
%             
%         end
%     end
%     IR=IRmove;
%     x
%  end

VarG=sum(sum( (IR-IG).^2))+sum(((IR-IG)'.^2));
   i0=2*displacementG(1);j0=2*displacementG(2);
   IGmove=IG;
    for i=range
        for j=range
            seG=translate(strel(1),[i0+i,j0+j]);
            IGTrymove=imdilate(IG,seG);
            VarGMOVE=sum(sum((IR-IGTrymove).^2))+sum(sum((IR-IGTrymove)'.^2));
            if  VarGMOVE< VarG 
                VarG=VarGMOVE;
                IGmove=IGTrymove;
                displacementG=[i0+i j0+j];
            end
        end
    end
    IG=IGmove;


 VarB=sum(sum( (IR-IB).^2))+sum(((IR-IB)'.^2));
   a0=2*displacementB(1);b0=2*displacementB(2);
   IBmove=IB;
    for a=range
        for b=range
            seB=translate(strel(1),[a0+a,b0+b]);
            IBTrymove=imdilate(IB,seB);
            VarBMOVE=sum(sum((IR-IBTrymove).^2))+sum(sum((IR-IBTrymove)'.^2));
            if  VarBMOVE< VarB 
                VarB=VarBMOVE;
                IBmove=IBTrymove;
                displacementB=[i0+i j0+j];
            end
        end
    end
    IB=IBmove;
x
 end
 Isum=cat(3,IR,IG,IB);
 imshow(Isum);