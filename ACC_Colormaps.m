
%% color map 1: blue to gray to red
n = 30;
m = 4;

c1 = linspace(0,0.9,n).^m;
c2 = linspace(0.4,0.9,n).^m;
c3 = linspace(0.8,0.9,n).^m;

c4 = linspace(0.9,0.8,n).^m;
c5 = linspace(0.9,0,n).^m;
c6 = linspace(0.9,0,n).^m;

colormap1 = [c1',c2',c3';c4',c5',c6'];
%% 

%% color map 2: off white to dark red
n = 60;
m = 2;

c1 = linspace(0.9,0.8,n).^m;
c2 = linspace(0.9,0,n).^m;
c3 = linspace(0.9,0,n).^m;

colormap2 = [c1',c2',c3'];
%%

%% color map 3: blue to dark red
n = 60;
m = 1;

c1 = linspace(0,0.6,n).^m;
c2 = linspace(0.4,0,n).^m;
c3 = linspace(0.8,0,n).^m;

colormap3 = [c1',c2',c3'];
%%

%% color map 4: off white to dark blue
n = 60;
m = 2;

c1 = linspace(0.9,0,n).^m;
c2 = linspace(0.9,0,n).^m;
c3 = linspace(0.9,0.6,n).^m;

colormap4 = [c1',c2',c3'];
%%

%% color map 5: white to dark blue
n = 60;
m = 1;

c1 = linspace(1,0,n).^m;
c2 = linspace(1,0,n).^m;
c3 = linspace(1,0.6,n).^m;

colormap5 = [c1',c2',c3'];
%%

%% color map 1: blue to dark blue to red (good for noticing small deviations from zero)
n1 = 40;
n3 = 100;

m = 1.5;

c1 = linspace(0.3,0.9,n1).^m;
c2 = linspace(0.3,0,n1).^m;
c3 = linspace(0.9,0,n1).^m;

c7 = linspace(0.9,0.97,n3).^m;
c8 = linspace(0,0.97,n3).^m;
c9 = linspace(0,0.97,n3).^m;

colormap6 = 0.9.*[c1',c2',c3';c7',c8',c9'];
%% 

%% clear
clear c1
clear c2
clear c3
clear c4
clear c5
clear c6
clear c7
clear c8
clear c9
clear m
clear n
clear n1
clear n2
clear n3
%%