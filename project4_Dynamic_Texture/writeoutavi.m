
load(strcat(num2str(1),'_synth_Result.mat'));
%Y = zeros(120*160,length(mov));
v = VideoWriter('1_test_1.avi');
open(v);
I=[];


for k = 1:length(synth_Result)
%     %I = rgb2gray(mov(k).cdata);
     I(:,:,1,k) = synth_Result(k).frame;
end
    writeVideo(v,I);
    %Y(:,k) = I(:);
% end
close(v)