close all

tam = 40; % Image size

% Creating both images
gray_levels = 512;
B0_gray1 = 0.7;
B0_gray2 = 0.5;
V1_gray1 = 0.3;
V1_gray2 = 0.15;

E1 = [1.0000    0.3000    0.4000         0.5         0         0;...
    B0_gray1    0.3000    0.4000         -0.5       0.2         0;...
    B0_gray2    0.2000    0.2000         -0.2      -0.6         0];
[P1] = phantom(tam,E1)*gray_levels;
figure, imshow(P1,[0,511]), axis image, colormap(gray), title('B0')

E2 = [0.8000    0.3000    0.4000         0.5         0         0;...
    V1_gray1    0.3000    0.4000         -0.5       0.2         0;...
    V1_gray2    0.2000    0.2000         -0.2      -0.6         0];
[P2] = phantom(tam,E2)*gray_levels;
figure, imshow(P2,[0, 511]), axis image, colormap(gray), title('V1')

idx_white = find(P1 == gray_levels);
idx_black = find(P1 == 0);
idx_gray1 = find(P1 == B0_gray1*gray_levels);
idx_gray2 = find(P1 == B0_gray2*gray_levels);

% Noise addition
noise_factor = 25;
P1_white = round(P1(idx_white) - rand(size(P1(idx_white)))*noise_factor*3);
P1_black = round(P1(idx_black) + rand(size(P1(idx_black)))*noise_factor*2);
P1_gray1 = round(P1(idx_gray1) + rand(size(P1(idx_gray1)))*noise_factor*4);
P1_gray2 = round(P1(idx_gray2) + rand(size(P1(idx_gray2)))*noise_factor*4);

P1(idx_white) = P1_white;
P1(idx_black) = P1_black;
P1(idx_gray1) = P1_gray1;
P1(idx_gray2) = P1_gray2;
figure(1), imshow(P1,[0, 511]), axis image, colormap(gray), title('B0')

P2_white = round(P2(idx_white) - rand(size(P2(idx_white)))*noise_factor*3);
P2_black = round(P2(idx_black) + rand(size(P2(idx_black)))*noise_factor*2);
P2_gray1 = round(P2(idx_gray1) + rand(size(P2(idx_gray1)))*noise_factor*4);
P2_gray2 = round(P2(idx_gray2) + rand(size(P2(idx_gray2)))*noise_factor*4);

P2(idx_white) = P2_white;
P2(idx_black) = P2_black;
P2(idx_gray1) = P2_gray1;
P2(idx_gray2) = P2_gray2;
figure(2), imshow(P2,[0, 511]), axis image, colormap(gray), title('V1')

% Selecting the pixel intensities
VoxIntensity = cell(3,3);
VoxIntensity{1,1} = 'Gray1';
VoxIntensity{2,1} = 'Gray2';
VoxIntensity{3,1} = 'WhiteM';

VoxIntensity{1,2} = P1_gray1;
VoxIntensity{1,3} = P2_gray1;
VoxIntensity{2,2} = P1_gray2;
VoxIntensity{2,3} = P2_gray2;
VoxIntensity{3,2} = P1_white;
VoxIntensity{3,3} = P2_white;

% Scattering diagram
B0 = [VoxIntensity{1,2}; VoxIntensity{2,2}; VoxIntensity{3,2}];
V1 = [VoxIntensity{1,3}; VoxIntensity{2,3}; VoxIntensity{3,3}];
tag = [repmat({VoxIntensity{1,1}}, length(VoxIntensity{1,2}),1);...
       repmat({VoxIntensity{2,1}}, length(VoxIntensity{2,2}),1);...
       repmat({VoxIntensity{3,1}}, length(VoxIntensity{3,2}),1)];
figure(3);
h1 = gscatter(B0,V1,tag,'rbg','v^o',[],'off');
set(h1,'LineWidth',2)
legend('Gray 1','Gray 2','White','Location','NW')
% axis([0 255 0 255])

% Classification
% [X,Y] = meshgrid(linspace(0,gray_levels),linspace(0,gray_levels));
% X = X(:); Y = Y(:);
% [C,err,P,logp,coeff] = classify([X Y],[B0 V1],tag,'linear');
% 
% % Visualize the classification
% hold on, %grid on
% gscatter(X,Y,C,'rbg','.',1,'off');
% K = coeff(1,2).const;
% L = coeff(1,2).linear; 
% 
% f = sprintf('0 = %g + %g*x + %g*y',K,L);
% h2 = ezplot(f,[0 gray_levels+100 0 gray_levels+100]);
% set(h2,'Color','g','LineWidth',2)
% axis([0 gray_levels 0 gray_levels])
% xlabel('B0')
% ylabel('V1')
% hold off

P1_white
P1_black
P1_gray1
P1_gray2
P2_white
P2_black
P2_gray1
P2_gray2

% Unstandardized coeff
P = 0.012.*P1 + 0.040.*P2 - 18.9;
figure, imagesc(P), axis image, colormap(gray), title('P')
% Fisher coeff
P = 0.012.*P1 + 0.040.*P2 - 18.9;
figure, imshow(P, []), axis image, colormap(gray), title('P')


