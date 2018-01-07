
function plot_exact_solution3(time)


Globals1D;

%plotting exact solution at all of the nodes of the mesh

x = linspace(a,b,K*(p+1));

x_star = x - time./(1+x.^2);

q0 = @(x)exp(-100*(x-0.3).^2);

plot(x,q0(x_star),'r*','linewidth',1)

title(['$t=',num2str(time),'$'],'interpreter','latex','fontsize',18)
xlabel('$x$','interpreter','latex','fontsize',18)

drawnow




