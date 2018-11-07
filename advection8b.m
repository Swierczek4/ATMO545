clc
clear
close all
tic()

ACC_Settings
lw = 2.2;
sz = 8;
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
P = Q;
R = 0.01;
H = zeros(1,360);
H(10) = 1;

A_init = sqrwv(X,-180,179);

%% forecast
A_fb = mvnrnd(A_init,v);
[Rho,D] = gcorr('gauss',Lx,Ld,Nx,1);
%% 

%% plot of IC
figure()
set(gcf, 'Position', [1, 1, 1280, 720])
plot(X,A_init,'LineWidth',lw,'Color',Color(:,color1))
hold on
plot(X,A_fb,'LineWidth',lw,'Color',Color(:,color2))
xlabel('x')
ylabel('y')
title('true state vs. analysis state at t = 0 s')
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
vidObj = VideoWriter('8b.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);
writeVideo(vidObj, getframe(gcf));
%%

%% plot of time series
Ab = A_init;
for ii=1:nsteps
    Ab = M*Ab;
    A_fb = M*A_fb;
    P = M*P*M'+Q;
    y = H*Ab + normrnd(0,R);
    K = P*H'/(H*P*H'+R);
    A_fb = A_fb + K*(y-H*A_fb);
    
    h1 = plot(X,A_fb,'LineWidth',lw,'Color',Color(:,color2));
    hold on
    h2 = plot(X,Ab,'LineWidth',lw,'Color',Color(:,color1));
    h3 = plot(-170,y,'o','MarkerSize',sz,'Color',Color(:,color3));
    plot(-170,y,'o','MarkerSize',sz-2,'Color',Color(:,color3))
    xlabel('x')
    ylabel('y')
    title(['true state vs. analysis state at t = ',num2str(ii),' s'])
    axis(coords)
    legend([h2(1),h1(1),h3(1)],'true state','analysis','observation')
    hold off
    drawnow()
    writeVideo(vidObj, getframe(gcf));
end
close(vidObj);
%%

save Ab
save A_fb