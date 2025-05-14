load('vorticity_18deg.mat')

x = x - min(min(x)) + 0.5 + 5;
y = y - min(min(y)) + 1.5 + 5;

xq = linspace(0, 20, 512)
yq = linspace(0, 20, 512)

[Xq, Yq] = meshgrid(xq,yq);

Vq = interp2(x,y,omega,Xq,Yq,'spline', 0);

figure(1)
contourf(Vq,21)

wc=window(@gausswin,512,5);
wr=window(@gausswin,512,5);
[maskr,maskc]=meshgrid(wr,wc);
w=maskr.*maskc;
figure(2)
imagesc(w)

figure(3)
contourf(w.*Vq,21)

Vq = w .* Vq;
output_matrix(Vq,512,xq-10,yq-10,'omega.bin')
