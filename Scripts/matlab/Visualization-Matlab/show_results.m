clear all;
close all;
clc;

addpath ./nifti/

%results = 'results/model1_offset/';
%results = 'results/model1_nooffset_easy/';

%results = 'results/model2_2/';
results = 'results/model1_offset2/';

model = '1';

% load volume
nii = load_nii([results 'deformed2_FA.nii.gz']);

% copy image into block matrix
b2 = nii.img;
%b2 = permute( b2, [2 1 3] );

% compute central slice
dd = ceil(size(b2,3)/2);


n=1;
m=5;
d = round(linspace( 25, 70, n*m ));

bachar = ['b' 'a'];

for ba = 1:2
    % load contours
    if ba == 1
        [vp, fp] = read_vtk([results  'csf_prior.vtk' ]);
        [vi, fi] = read_vtk([results  'wm_prior.vtk' ]);
    elseif ba == 2
        [vp, fp] = read_vtk([results  'deformed2-csf.vtk' ]);
        [vi, fi] = read_vtk([results  'deformed2-wm.vtk' ]);
    end
    
    for k=1:(n*m) % show a few slices
        clf;
        %subplot(n,m,k);
        depth = d(k);
        % show central slice
        imshow( b2(:,:,depth) );
        
        hold on;
        
        % show pial contour / only vertices close to center slice
        v2 = vp+1;
        tmp = v2(:,1);
        v2(:,1) = v2(:,2);
        v2(:,2) = tmp;
        v2( abs(v2(:,3) - depth) > 0.75, 3) = NaN; % threshold based on density
        p = patch('Faces', fp, 'Vertices', v2, 'EdgeColor', 'green', 'FaceColor', 'none');
        set(p, 'LineWidth', 2);
        
        % show gmwm contour
        v2 = vi+1;
        tmp = v2(:,1);
        v2(:,1) = v2(:,2);
        v2(:,2) = tmp;
        v2( abs(v2(:,3) - depth) > 0.75, 3) = NaN;
        p = patch('Faces', fi, 'Vertices', v2, 'EdgeColor', 'red', 'FaceColor', 'none');
        set(p, 'LineWidth', 2);
        
        print('-dpng', '-r300', [results 'model' model 'result_' bachar(ba) '_' num2str(k) '.png']);
        print('-depsc2', '-r300', [results 'model' model 'result_' bachar(ba) '_' num2str(k) '.eps']);
    end
    
end

%%
figure;
tmp = vp(:,1);
vp(:,1) = vp(:,2);
vp(:,2) = tmp;
p = patch('Faces', fp, 'Vertices', vp, 'FaceVertexCData', 1);
%isonormals(b2, p)
set(p, 'FaceColor', 'flat', 'EdgeColor', 'none');
set(p, 'FaceAlpha', 'flat', 'FaceVertexAlphaData', 0.25 );

tmp = vi(:,1);
vi(:,1) = vi(:,2);
vi(:,2) = tmp;
p = patch('Faces', fi, 'Vertices', vi, 'FaceVertexCData', 1);
%isonormals(b2, p)
set(p, 'FaceColor', 'flat', 'EdgeColor', 'none');

for k = 1:2
    v2 = vi;
    v2( abs(v2(:,3) - d(k)) > 0.75, 3) = NaN;
    p = patch('Faces', fi, 'Vertices', v2, 'EdgeColor', 'red', 'FaceColor', 'none');
    set(p, 'LineWidth', 2);
end

% for k = 1:2
%     v2 = vp;
%     v2( abs(v2(:,3) - d(k)) > 0.75, 3) = NaN;
%     p = patch('Faces', fp, 'Vertices', v2, 'EdgeColor', 'green', 'FaceColor', 'none');
%     set(p, 'LineWidth', 2);
% end

colormap(gray);
daspect([1 1 1]);
view(0,-90);
box off;
axis off;

camlight;
lighting gouraud;

print( '-dpng', '-r300', [results 'model' model 'surf.png']);
print( '-depsc2', '-r300', [results 'model' model 'surf.eps']);

