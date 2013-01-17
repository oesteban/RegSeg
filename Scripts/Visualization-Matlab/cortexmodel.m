clear all;
close all;
clc;


N = 30; % Radius of the WM ball, in pixels
        % don't go much higher than 50 -- complexity increases to the 6th power or so
M = ceil(N/5); % margin around the brain
T = 2*(floor(M/2)); % thickness of the cortex :> slightly smaller than the margin
L = 2*N+1; % diameter of the WM ball

I2 = ones(L+2*M); % size of one slice
[i,j] = find(I2); % just get the coordinates of every point

block = zeros( L+2*M,L+2*M,L+2*M ); % empty block where the brain will be

disp('Creating WM skeleton');

for depth = 1:(L+2*M)
    I = I2;
    r2 =  N.^2 - (depth-N-M-1).^2;
    I( (i-N-M-1).^2 + (j-N-M-1).^2 > r2 ) = 0; % create disc in slice
    I( (i-1) < N+M-T  & abs(j-N-M-1) < T ) = 0; % cutout
    
    imshow(I);
    block(:,:,depth) = I; % write slice to block
    drawnow;
end

disp('Rounding edges');

%se = strel('disk', T/2, 0);
se = strel_sphere( T/2 );
block = imopen( block, se );   % 3D math. morph. to round edges (giral ridge)

%se = strel('disk', T-1, 0);
se = strel_sphere( T-1 );
block = imclose( block, se );  % 3D MM to round sulcal bottom

disp('Adding cortical GM sheet');

%se = strel('disk', T, 0);
se = strel_sphere( T );
block = 0.5*block + 0.5*imdilate( block, se );  % add a uniformly thick GM layer

%%
disp('Lowpass filtering');

gauss = [1 3 1] / 5;

% slightly lowpass along all individual dimensions, to smoothen mesh
block = convn( shiftdim(block,1), gauss, 'same' );
block = convn( shiftdim(block,1), gauss, 'same' );
block = convn( shiftdim(block,1), gauss, 'same' );

n=4;
m=4;
for k=1:(n*m) % show a few slices
    d = round(linspace( 1, N+M+1, n*m ));
    subplot(n,m,k);
    imshow( block(:,:,d(k)) );
end

figure;
%%
disp('Extracting surfaces');

subplot(121);
FV = isosurface(block, 0.25); % get the pial surface between 0 and 0.5
FV.FaceVertexCData = 1;
p = patch(FV);
isonormals(block, p)
set(p, 'FaceColor', 'flat', 'EdgeColor', 'none');
daspect([1 1 1]);
view(0,-90);
box off;
axis off;
camlight;
lighting gouraud;


subplot(122);
FV = isosurface(block, 0.75); % get the GM/WM interface between 0.5 and 1
FV.FaceVertexCData = 1;
p = patch(FV);
isonormals(block, p)
set(p, 'FaceColor', 'flat', 'EdgeColor', 'none');
daspect([1 1 1]);
view(0,-90);
box off;
axis off;
camlight;
lighting gouraud;