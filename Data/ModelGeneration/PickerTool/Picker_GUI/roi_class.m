
B0 = [VoxIntensity{1,2}; VoxIntensity{2,2}; VoxIntensity{3,2}];
V1 = [VoxIntensity{1,3}; VoxIntensity{2,3}; VoxIntensity{3,3}];
tag = [repmat({VoxIntensity{1,1}}, length(VoxIntensity{1,2}),1);...
       repmat({VoxIntensity{2,1}}, length(VoxIntensity{2,2}),1);...
       repmat({VoxIntensity{3,1}}, length(VoxIntensity{3,2}),1)];
figure;
h1 = gscatter(B0,V1,tag,'rbg','v^o',[],'off');
set(h1,'LineWidth',2)
legend('Gray matter','White matter','CSF','Location','NE')

% Classification
[X,Y] = meshgrid(linspace(100,1100),linspace(0,250));
X = X(:); Y = Y(:);
[C,err,P,logp,coeff] = classify([X Y],[B0 V1],tag,'quadratic');

% Visualize the classification
hold on;
gscatter(X,Y,C,'rb','.',1,'off');
K = coeff(1,2).const;
L = coeff(1,2).linear; 
Q = coeff(1,2).quadratic;
f = sprintf('0 = %g+%g*x+%g*y+%g*x^2+%g*x.*y+%g*y.^2',...
            K,L,Q(1,1),Q(1,2)+Q(2,1),Q(2,2));
h2 = ezplot(f,[4.5 8 2 4]);
set(h2,'Color','m','LineWidth',2)
axis([4.5 8 2 4])
xlabel('Sepal Length')
ylabel('Sepal Width')
title('{\bf Classification with Fisher Training Data}')