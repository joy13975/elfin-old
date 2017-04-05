close all;
figure; 
xyrange = [0,400,0,400];
axis(xyrange); 
hold off;

grid on
box on
[x,y]=getpts();

plot(x,y,'-bx'); 
axis(xyrange)

P=[x,y];

D = (P-circshift(P, 1));
normaliseFactor = 50; % this should be set around mean pair-CoM distance
m = mean(sqrt(sum(D.^2, 2))) / normaliseFactor;

format long g
P=P./m;
P=[P, zeros(size(P,1), 1)]

cd('/Users/joy/src/elfin/bm/fun')

% save using dlmwrite('R3.csv', P*3, ' ')