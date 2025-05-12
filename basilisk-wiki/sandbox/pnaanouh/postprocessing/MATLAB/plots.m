close all

Ls=2*11.4411000000000;
D=1;

%%%%%%%%%% PHI=0.15, MuR=0.4237895, Bo=0.5, RhoR=1.1668599 
[t10,t_bins10,n_int_start10,t_range10,Dist_ref10,vx10,vy10,vt10,DS10,t_cont10,t_min10,Dist_free10,Dist_not_free10,n_min10] = drop_statistics(Ga10,Ls,D,20,121);
[vxrms10,vyrms10,vxp10,vyp10,vxm10,vym10]=V_fluctuations(Ga10,Ga10gl,Ls);
[gav2d10, X,Y] = r_distribution2D(Ga10,Ls,D,3,0.01);
[gav1d10,radii] = r_distribution(Ga10,Ls,D,3,0.01);

figure('Position',[0 0 1080 1080]) 
hold on

plot((Dist_ref10-1)./(2.458861094-1),(n_int_start10./(Ls*Ls*t10(length(t10)))),'DisplayName',"f_c \mu_r=0.42, \phi=0.15, Ga=10",'Color','#0072BD','LineStyle','-','LineWidth',1);
plot((Dist_ref10-1)./(2.458861094-1),n_min10./(2*Ls*Ls*t10(length(t10))),'DisplayName',"\omega_c \mu_r=0.42, \phi=0.15, Ga=10",'Color','#0072BD','LineStyle','--','LineWidth',1);
 
legend('Location','northwest','FontSize',22)
xlabel("$\frac{G-D}{dp-D}$",'Interpreter','latex','FontSize',25)
ylabel("$Frequency (\frac{1}{t L^2})$",'Interpreter','latex','FontSize',25)
grid on
grid minor
xlim([-1 2]);

savefig('./figures/freq_param_1.fig');
h=gcf;
exportgraphics(h,'./figures/freq_param_1.pdf','ContentType','vector','Resolution',300);


B=-1:0.005:1;
x=-1:0.001:1;
figure('Position',[0 0 1080 1080])
histogram(vyp10,B,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1,'DisplayName','Histogram');
hold on
plot(x, exp((-0.5.*((x-mean(vyp10,'all'))./std(vyp10,0,'all')).^2))./(std(vyp10,0,'all')*sqrt(2.*pi)),'LineWidth',1,'DisplayName','Gaussian')
xlabel("$V'_y$",'Interpreter','latex','FontSize',25)
ylabel("$Probability Density$",'Interpreter','latex','FontSize',25)
legend('Location','northeast','FontSize',22)
savefig('./figures/Vp_y_10.fig');
h=gcf;
exportgraphics(h,'./figures/Vp_y_10.pdf','ContentType','vector','Resolution',300);

figure('Position',[0 0 1080 1080])
histogram(vxp10,B,'Normalization','pdf','DisplayStyle','stairs','LineWidth',1,'DisplayName','Histogram');
hold on
plot(x, exp((-0.5.*((x-mean(vxp10,'all'))./std(vxp10,0,'all')).^2))./(std(vxp10,0,'all')*sqrt(2.*pi)),'LineWidth',1,'DisplayName','Gaussian');
xlabel("$V'_x$",'Interpreter','latex','FontSize',25)
ylabel("$Probability Density$",'Interpreter','latex','FontSize',25)
legend('Location','northeast','FontSize',22)
savefig('./figures/Vp_x_10.fig');
h=gcf;
exportgraphics(h,'./figures/Vp_x_10.pdf','ContentType','vector','Resolution',300);


clear x B


figure('Position',[0 0 1920 720]) 
hold on 
plot(t_range10+0.25,t_bins10(:,6)*2,'-k','DisplayName',"\mu_r=0.42, \phi=0.15, Ga=10",'Marker','o')

legend('Location','northeast','FontSize',22)
xlabel("$\frac{t}{\sqrt{D/g}}$",'Interpreter','latex','FontSize',28)
ylabel("$Probability Density$",'Interpreter','latex','FontSize',28)
xlim([0 30])
grid on
grid minor

savefig('./figures/tcont_param_1.fig');
h=gcf;
exportgraphics(h,'./figures/tcont_param_1.pdf','ContentType','vector','Resolution',300);


figure('Position',[0 0 1080 1080]) 
hold on 
plot(t10,vt10./vt10(1),'-k','DisplayName',"\mu_r=0.42, \phi=0.15, Ga=10",'LineWidth',1)

legend('Location','northwest','FontSize',22)
xlabel("$\frac{t}{\sqrt{D/g}}$",'Interpreter','latex','FontSize',25)
ylabel("$Volume/V_{init}$",'Interpreter','latex','FontSize',25)
grid on
grid minor

savefig('./figures/V_cons_param_1.fig');
h=gcf;
exportgraphics(h,'./figures/V_cons_param_1.pdf','ContentType','vector','Resolution',300);

figure('Position',[0 0 1920 720]) 
hold on
plot(t10,vx10-Ga10gl(1:size(Ga10gl,1)-1,6),'-k','DisplayName',"\mu_r=0.42, \phi=0.15, Ga=10",'LineWidth',1)

legend('Location','southeast','FontSize',25)
xlabel("$\frac{t}{\sqrt{D/g}}$",'Interpreter','latex','FontSize',28)
ylabel("$\langle V_x \rangle_p$",'Interpreter','latex','FontSize',28)
grid on
grid minor

savefig('./figures/V_x_part_param_1.fig');
h=gcf;
exportgraphics(h,'./figures/V_x_part_param_1.pdf','ContentType','vector','Resolution',300);

figure('Position',[0 0 1920 720]) 
hold on
plot(t10,vy10-Ga10gl(1:size(Ga10gl,1)-1,7),'-k','DisplayName',"\mu_r=0.42, \phi=0.15, Ga=10",'LineWidth',1)

legend('Location','southeast','FontSize',25)
xlabel("$\frac{t}{\sqrt{D/g}}$",'Interpreter','latex','FontSize',28)
ylabel("$\langle V_y \rangle_p$",'Interpreter','latex','FontSize',28)
grid on
grid minor

savefig('./figures/V_y_part_param_1.fig');
h=gcf;
exportgraphics(h,'./figures/V_y_part_param_1.pdf','ContentType','vector','Resolution',300);

figure('Position',[0 0 1920 720]) 
hold on
plot(t10,Ga10gl(1:size(Ga10gl,1)-1,6),'-k','DisplayName',"\mu_r=0.42, \phi=0.15, Ga=10",'LineWidth',1)

legend('Location','southwest','FontSize',25)
xlabel("$\frac{t}{\sqrt{D/g}}$",'Interpreter','latex','FontSize',28)
ylabel("$Vf_x$",'Interpreter','latex','FontSize',28)
grid on
grid minor

savefig('./figures/V_x_fluid_param_1.fig');
h=gcf;
exportgraphics(h,'./figures/V_x_fluid_param_1.pdf','ContentType','vector','Resolution',300);

figure('Position',[0 0 1920 720]) 
hold on
plot(t10,Ga10gl(1:size(Ga10gl,1)-1,7),'-kx','DisplayName',"\mu_r=0.42, \phi=0.15, Ga=10",'LineWidth',1,'MarkerIndices',1:1000:length(t10))

legend('Location','southwest','FontSize',25)
xlabel("$\frac{t}{\sqrt{D/g}}$",'Interpreter','latex','FontSize',28)
ylabel("$Vf_y$",'Interpreter','latex','FontSize',28)
grid on
grid minor

savefig('./figures/V_y_fluid_param.fig');
h=gcf;
exportgraphics(h,'./figures/V_y_fluid_param.pdf','ContentType','vector','Resolution',300);

figure('Name','Surface','Position',[0 0 600 920])

surf(X,Y,gav2d10,'EdgeColor','none','FaceColor','flat');
view(2);
c=colorbar;
c.Label.String='Pair Distribution Function';
c.Label.FontSize =20;
colormap jet

xlabel("$\frac{X}{D}$",'Interpreter','latex','FontSize',22)
ylabel("$\frac{Y}{D}$",'Interpreter','latex','FontSize',22)
grid on
grid minor

savefig('./figures/pair_disribution_2D_10.fig');
h=gcf;
exportgraphics(h,'./figures/pair_disribution_2D_10.tiff','ContentType','image','Resolution',300);

figure('Position',[0 0 1080 1080]) 
hold on
plot(radii,gav1d10,'-k','DisplayName',"\mu_r=0.42, \phi=0.15, Ga=10",'LineWidth',1);

legend('Location','southeast','FontSize',22)
ylabel("$Radial Distribution$",'Interpreter','latex','FontSize',25)
xlabel("$\frac{r}{D}$",'Interpreter','latex','FontSize',25)
grid on
grid minor

savefig('./figures/pair_disribution_1D_01.fig');
h=gcf;
exportgraphics(h,'./figures/pair_disribution_1D_01.pdf','ContentType','vector','Resolution',300);
