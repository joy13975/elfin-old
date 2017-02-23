close all;
figure; 
xyrange = [0,400,0,400];
axis(xyrange); 
hold off;

[x,y]=getpts();

plot(x,y,'-bx'); 
axis(xyrange)

P=[x,y];

D = (P-circshift(P, 1));
m = mean(sqrt(sum(D.^2, 2))) / 50;

format long g
P=P./m;
P=[P, zeros(size(P,1), 1)]

