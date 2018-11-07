clc
clear
close all
tic()

ACC_Settings
lw = 2.2;
color1 = 15;
color2 = 11;
color3 = 9;
coords = [-180 179 -2 3];

X = -180:179;
mu = 1;
dx = 1;
dt = 1;
Nx = length(X);
Lx = Nx*dx;
C = 0.5*mu*dt/dx;
v = 0.5;
Ld = 45;
nsteps = 300;

v1 = (1-4*C^2).*ones(Nx,1);
v2 = (2*C^2-C).*ones(Nx-1,1);
v3 = (C+2*C^2).*ones(Nx-1,1);

M = diag(v1) + diag(v2,1) + diag(v3,-1);
M(1,Nx) = C+2*C^2;
M(Nx,1) = 2*C^2-C;
Q = 0.5.*eye(360);

A_init = sqrwv(X,-180,179);

%% forecast
A_f = mvnrnd(A_init,v);
[Rho,D] = gcorr('gauss',Lx,Ld,Nx,1);
%% 

%% plot of IC
figure()
set(gcf, 'Position', [1, 1, 1280, 720])
h2 = plot(X,A_f,'LineWidth',lw,'Color',Color(:,color2));
hold on
h1 = plot(X,A_init,'LineWidth',lw,'Color',Color(:,color1));
xlabel('x')
ylabel('y')
title('true state vs. forecast state at t = 0 s')
legend([h1(1),h2(1)],'true state','forecast')
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1) + 0.01;
bottom = outerpos(2) + ti(2) + 0.01;
ax_width = outerpos(3) - ti(1) - ti(3) - 0.03;
ax_height = outerpos(4) - ti(2) - ti(4) - 0.02;
ax.Position = [left bottom ax_width ax_height];
axis(coords)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','on');
vidObj = VideoWriter('8a.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);
writeVideo(vidObj, getframe(gcf));
%%

%% plot of time series
A = A_init;
for ii=1:nsteps
    A = M*A;
    A_f = M*A_f;
   
    h2 = plot(X,A_f,'LineWidth',lw,'Color',Color(:,color2));
    hold on
    h1 = plot(X,A,'LineWidth',lw,'Color',Color(:,color1));
    xlabel('x')
    ylabel('y')
    title(['true state vs. forecast state at t = ',num2str(ii),' s'])
    axis(coords)
    legend([h1(1),h2(1)],'true state','forecast')
    hold off
    drawnow()
    writeVideo(vidObj, getframe(gcf));
end
close(vidObj);
%%

save A
save A_f