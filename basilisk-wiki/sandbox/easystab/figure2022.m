figure;
a = 1;
x1s = [-2:0]; y1s =a*x1s;
x1a = [0,3]; y1a = a*x1a;
plot(x1s,y1s,'k','Linewidth',3); hold on;
plot(x1a,y1a,'k--','Linewidth',3); 
x1s = [3:4]; y1s =a*x1s;
plot(x1s,y1s,'k','Linewidth',3); hold on;


x2 = [-1:.01:3];
y2 = 1+sqrt(x2+1);
x2a = [3:.1:4];
y2a = 1+sqrt(x2a+1);
plot(x2,y2,'k','Linewidth',3);
plot(x2a,y2a,'k--','Linewidth',3);

x3a = [-1:.01:0];
y3a = 1-sqrt(x3a+1);
plot(x3a,y3a,'k--','Linewidth',3);
x3 = [0:.01:4];
y3 = 1-sqrt(x3+1);
plot(x3,y3,'k','Linewidth',3);

axis off;

