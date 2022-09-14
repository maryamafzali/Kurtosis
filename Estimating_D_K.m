clear
close all
clc

Delta = 0.023;
delta = 0.012; %[seconds]
Gamma = 2.675987E8; % rad/s/T
D0 = 3000; % [micrometer^2/seconds]

bval = load('20211104_114911LRTE55bmax15000multipledirectionss004a001.bval');
bvec = load('20211104_114911LRTE55bmax15000multipledirectionss004a001.bvec');
S_img = load_untouch_nii('real_4_11_2021_multi_gibbsCorrSubVoxShift_driftCo_TED.nii.gz');
S_img = double(S_img.img);
S_img = S_img(end:-1:1,:,:,:);

mask = load_untouch_nii('real_4_11_2021_multi_gibbsCorrSubVoxShift_driftCo_TED_brain_mask.nii.gz');
mask = double(mask.img);
mask = mask(end:-1:1,:,:,:);
se = strel('disk',2);
mask = imerode(mask,se);

[b_sort, order] = sort(bval);
[C, ia, ic] = unique(b_sort);

q = sqrt(C/10^6/(Delta - delta/3));

p = length(C);
S0 = S_img(:,:,:,bval == 0);

S_shell = zeros(size(S_img,1),size(S_img,2),size(S_img,3),p-1);

S_img_sort = S_img(:,:,:,order);
for kk = 1:p-1
    S_shell(:,:,:,kk) = mean(S_img_sort(:,:,:,ia(kk):ia(kk+1)-1),4).*mask;
end
S_shell(:,:,:,p) = mean(S_img_sort(:,:,:,ia(p):end),4).*mask;
S_shell(S_shell<0) = 10^-3;

S_norm = zeros(size(S_shell,1),size(S_shell,2),size(S_shell,3),size(S_shell,4));

for k = 1:size(S_shell,4)
    S_norm(:,:,:,k) = S_shell(:,:,:,k).*mask./(S_shell(:,:,:,1));
end

S_norm(isnan(S_norm)) = 10^-3;
S_norm(S_norm == inf) = 10^-3;
S_norm(S_norm >1) = 1;
S_norm(S_norm <0) = 10^-3;

[s1, s2, s3, s4] = size(S_norm);

xdata(:,1) = q;
xdata(1,2) = Delta;
xdata(2,2) = delta;

opts = optimset('Display','off');
lb = [0.01, 0.01];
ub = [3, 2];
x1 = [2, 1];

map_d_k = zeros(s1,s2,s3,2);
for i = 1:s1
    for j = 1:s2
        for k = 37 % 1:s3
            if mask(i,j,k)>0
                signal = squeeze(S_norm(i,j,k,:));
                map_d_k(i,j,k,:) = lsqcurvefit(@kurtosis_d, x1, xdata(1:5,:), signal(1:5,:) ,lb,ub,opts);
            end
        end
    end
    i
end
