
names = {'B0'
         'V20'
         'V26'
         'V33'
         'V36'
         'V37'
         'V43'
         'V45'
         'V48'
         'V53'
         'V55'
         'V57'};

imgs = zeros(size(names));
for i = 1:numel(names)
        nm = names{i};
        nm(1) = [];
        imgs(i) = str2double(nm);
end
 
f1 = [0.947
      0.858
     -0.540
     -0.222
     -0.312
      0.805
     -0.141
     -1.101
      0.344
     -0.280
      0.378
      0.061];
  
f2 = [-0.113
       0.788
      -1.117
      -1.096
       1.035
      -0.920
       1.242
       0.574
      -1.645
       1.216
       1.258
      -1.155];
  
f1 = [0,014
-0,020
0,041
0,034
0,004
-0,071
-0,038
0,052
0,051
-0,019
-0,009
0,038
-7,990
];

[a b c d] = size(nii.img);
img1 = zeros(a,b,c);
nii.img = double(nii.img);
for k = 1:numel(imgs)
    img1 = nii.img(:,:,:,imgs(k)+1).*f2(k) + img1;
end

img1 = nii.img(:,:,:,imgs(k)+1).*f1(k+1) + img1;
% img1 = img1*(-1);
nii2 = nii;
nii2.img = img1;
nii2.hdr.dime.dim(5)=1;

