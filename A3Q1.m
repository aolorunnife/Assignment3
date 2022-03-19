
%part 1
clearvars; clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')  % 'docked' 'normal'
set(0,'DefaultLineLineWidth',1)

%%scattering

%solve for Vth
mo= 9.1093837015E-31;
mn = 0.26*mo;
%may need to an array to do this
l= 200E-9;  %come back into the other side
h= 100E-9;  %bounce back
a= l*h;
%%% ask TA
T= 300;
k= 1.38064852E-23;  %for thermal velocity
Tmn= 0.2E-12;  %mean free path

Vth = sqrt((2*k*T)/(mn));
np = 10000; %number of particles
partWatch = randi(np,7,1);  %an array of indices to select a sample of particles
%mean free path
mfp = Vth*Tmn;

% v(t) = ds/dt
% v(t) = v0 + a*t  %velocity
% s(t) = s0 + v*t + 1/2 a t^2  %position

X = rand(np,1)*l;  %random number between 0 and1 multiplied by limits of silicon
Y = rand(np,1)*h;

X1 = X;
Y1 = Y;
std=sqrt(k*T/mn);

Vx = std*(randn(np,1));  %x component
Vy = std*(randn(np,1));  %y component

avgV = sqrt(mean(Vx.^2 + Vy.^2));  %vx and vy are vectors
semiT = (avgV).^2*mn/(2*k);

dt = h/Vth/100;  %should be 1/100 of the region size


numit=1000;
Xp = X;  %when it hits the top, X remains the same
%when hitting the x max, come back from left so X-width
Yp = Y;   % when it hits the top, direction will just be -Vy



%%1 a) 0.1V across x direction
V = 0.2;
Ex = -V/l;
Ey = 0;

%in x direction

%% 1 b) for on eache electron

q = -1.602*10^-19;
Fex = q*Ex;
Fey = 0;

Area = l * h;
EC = 10E15;

numarea = Area * EC;  %total number of electrons in reggion

numelec = numarea/np;

%%inside the lopp, count particles that go from left to right and how many
%%go fromright to left
%take differece and get value over delta t and caluclate current
%track current then plot it
%current as a function of time (plot in x direction)



%% acceleration and plot

accx = Fex/mn;
accy = Fey/mn;


pathnum = 0;
distancesum = 0;
for i=1:numit


    Xp = X;
    Yp = Y;

    X = X + dt*Vx + 1/2*(accx*(dt^2)) ;
    Y = Y + dt*Vy ;
    Vx = Vx + accx*dt;
    Vy = Vy + accy*dt;

    %left
    ix = X < 0;
    X(ix) = X(ix)+l;
    Xp(ix) = Xp(ix) + l;
    Y = Y;
        
    a = sum(ix);  %number of particles going left

    %right
    ix = X > l;
    X(ix) = X(ix)-l;
    Xp(ix) = Xp(ix)-l;
    Y = Y;

    b = sum(ix); %particles going to the right


    %%drift current
    partpertime = (b - a) / dt ;
    currx(i) = q*EC*partpertime;

    iy = Y<0 | Y > h;
    Vy(iy) = -Vy(iy);

    %scattering, ONLY PLOTS A FEW
    std=sqrt(k*T/mn);
    Pscat = 1 - exp(-dt/Tmn);
    iscat = Pscat > rand(np,1);

    X = X + dt*Vx;
    Y = Y + Vy*dt;
    Vx(iscat) = std*randn(sum(iscat),1);
    Vy(iscat) = std*randn(sum(iscat),1);


    %mean free path
    pathnum = pathnum + sum(iscat);
    distance = sqrt((X1(iscat)-X(iscat)).^2 + ((Y1(iscat)-Y(iscat)).^2));
    distancesum = distancesum + sum(distance);
    X1(iscat) = X(iscat);


    %calculating temperarure using avg velocity

    avgV = sqrt(mean(Vx.^2 + Vy.^2));  %vx and vy are vectors
    semiT = (avgV).^2*mn/(2*k);

    figure(1);
    hold on
    axis([0 l 0 h])

    %titlep = sprintf('Temp = %.3d', semiT);

    title(sprintf('Particle Trajectories, Temp = %.3d', semiT))
    xlabel('x (m)')
    ylabel('y (m)')

    for j=1:length(partWatch)

        plot([Xp(partWatch(j)),X(partWatch(j))]',[Yp(partWatch(j)),Y(partWatch(j))]','SeriesIndex',j)

    end

end


   figure (2)
    hold on
    xlabel('Time (s)');
    ylabel ('Current (A)')
    title('Drift Current over Time');
    plot(currx);


%density plot and %temperature map

binsx = discretize(X, 50); %20 bins
binsy = discretize(Y, 50);

NumPartBin = zeros(50,50);
tempmap = zeros(10,10);
for x = 1:50

    for y= 1:50
        NumPartBin(x,y) =  sum(binsx==x & binsy==y);  %gives u number of particles
        avgV = mean(sqrt(Vx(binsx==x& binsy==y).^2 + Vy(binsx==x & binsy==y).^2)); %vx and vy are vectors
        tempmap(x,y) = (avgV).^2*mn/(2*k);
   
    end

end
figure (3);
surf(tempmap)
xlabel('x')
ylabel('y')
zlabel('Temperature (K)')
title('Average Temperature Over time')

figure (4);
surf(NumPartBin/np)
title('Density Plot')
xlabel('X(m)')
ylabel('Y(m')
zlabel('Electron Count')
title('Electron Density Map')

avgV = mean(sqrt(Vx.^2 + Vy.^2)); %vx and vy are vectors
semiT = (avgV).^2*mn/(2*k);

% 
% current = start off electrons are random so curr = 0
% as we continue, E feild accelerates, more curves
% move particles travelling in different direction, so current increases 
% transient as electrons speed up and gets limited by scattering so it saturates
% after so many iterations, current saturates (constant current) around -12/14 A
