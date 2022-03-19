clearvars; clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')  % 'docked' 'normal'
set(0,'DefaultLineLineWidth',1)

Length = 600E-9; %x direction
W = 400E-9; %y direction
vo = 1;

for f = 1:12


nx = 30;
ny = 20;
COut = 1; %conductivity outside boxs
CIn = 10^-2; %conductivity inside boxs
box = 1; %making box

V = zeros(nx,ny);
F = zeros(nx*ny,1);
G = sparse(nx*ny,nx*ny);  %matrix
CM = zeros(nx,ny);

if box == 1
    b1 = rectangle('Position',[2.0E-7,0,1.8E-7*f,1.8E-7*f]); %box1 position (bottom)
    b2 = rectangle('Position',[2.0E-7,2.2E-7,1.8E-7*f,1.8E-7*f]); %box2 position (top)
end


for j = 1:ny
    for i = 1:nx
        newW = j*(W/ny);
        newL = i*(Length/nx);

        if box == 1 && (newW > 2.2E-7*f || newW < 1.8E-7*f) && box == 1 && newL > 2.0E-7*f && newL < 3.95E-7*f
            CM(i,j) = CIn;
        else
            CM(i,j) = COut;
        end
    end
end


eq = 1;
for j = 1:ny
    for i = 1:nx
        n = j + (i-1)*ny;

        if i == 1
            F(n) = vo;
            G(n,n) = CM(i,j);

        elseif i == nx
            if eq == 1
                F(n) = 0;
                G(n,n) = CM(i,j);
            else
                F(n) = vo;
                G(n,n) = CM(i,j);
            end

        elseif j == 1
            if eq == 1
                F(n) = 0;
                nxm = j + ((i-1)-1)*ny; %(i-1,j)
                nyp = (j+1) + (i-1)*ny; %(i,j+1)
                nxp = j + ((i+1)-1)*ny; %(i+1,j)
                G(n,n) = -(CM(i-1,j) + CM(i+1,j) + CM(i,j+1));
                G(n,nxm) = CM(i-1,j);
                G(n,nxp) = CM(i+1,j);
                G(n,nyp) = CM(i,j+1);
            else
                F(n) = 0;
                G(n,n) = CM(i,j);
            end


        elseif j == ny
            if eq == 1
                F(n) = 0;
                nxm = j + ((i-1)-1)*ny; %(i-1,j)
                nxp = j + ((i+1)-1)*ny; %(i+1,j)
                nym = (j-1) + (i-1)*ny; %(i,j-1)
                G(n,n) = -(CM(i-1,j) + CM(i+1,j) + CM(i,j-1));
                G(n,nxm) = CM(i-1,j);
                G(n,nxp) = CM(i+1,j);
                G(n,nym) = CM(i,j-1);
            else
                F(n) = 0;
                G(n,n) = CM(i,j);
            end

        else
            nxm = j + ((i-1)-1)*ny;
            nxp = j + ((i+1)-1)*ny;
            nym = (j-1) + (i-1)*ny;
            nyp = (j+1) + (i-1)*ny;

            G(n,n) = -(CM(i-1,j) + CM(i+1,j) + CM(i,j-1) + CM(i,j+1));
            G(n,nxm) = CM(i-1,j);
            G(n,nxp) = CM(i+1,j);
            G(n,nym) = CM(i,j-1);
            G(n,nyp) = CM(i,j+1);

        end
    end
end

P = G\F;

for j = 1:ny
    for i = 1:nx
        n = j + (i-1)*ny;
        V(i,j) = P(n);
    end
end


dx = Length/nx;
[Ey, Ex] = gradient(V);
Ey = Ey/dx;
Ex = Ex/dx;


Jx = -CM.*Ex;
Jy = -CM.*Ey;
depth = 1;


    %Current
    A=depth.*W;

    Ix(f)=mean(Jx(1,:)).*A;

end

%plotting current vs bottleneck
plot((1:12), Ix)  
figure (1)
title('Current vs BottleNeck')
xlabel('Distance between Boxes')
ylabel('Current (A)')

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


Area = l * h;
EC = 10E15;

numarea = Area * EC;  %total number of electrons in reggion

numelec = numarea/np;

lx = linspace(0,Length,nx);
ly = linspace(0,W,ny);

[LX,LY]= meshgrid(lx,ly);
% figure
% surf(LX,LY,Ex.');
% interp2(LX,LY,Ex.',2e-7,3e-7);


figure(2);
hold on
axis([0 l 0 h])

box1 = rectangle('Position',[0.8E-7, 0, 0.4E-7, 0.4E-7]);
box2 = rectangle('Position',[0.8E-7, 0.6E-7 ,0.4E-7,0.4E-7]);

InBox = X > 0.8E-7*f & X < 1.2E-7*f & (Y > 0.6E-7*f | Y < 0.4E-7*f);
while sum(InBox)> 0

    X(InBox) = rand(sum(InBox),1)*l;
    Y(InBox) = rand(sum(InBox),1)*h;
    InBox = X > 0.8E-7*f & X < 1.2E-7*f & (Y> 0.6E-7*f | Y < 0.4E-7*f);
end


pathnum = 0;
distancesum = 0;
for i=1:numit

    q = -1.602*10^-19;
    
   Ex_p= interp2(LX,LY,Ex.',Xp, Yp);
   Ey_p= interp2(LX,LY,Ey.',Xp, Yp);
   Ex_p(isnan(Ex_p)) = 0;
   Ey_p(isnan(Ey_p)) = 0;

    Fex = q*Ex_p;
    Fey = q*Ey_p;

    accx = Fex/mn;
    accy = Fey/mn;

    Xp = X;
    Yp = Y;

    X = X + dt*Vx + 1/2*(accx*(dt^2)) ;

    % map particle position to i,j in acc matrix
    %divide into equal parts
    %i = x
    Y = Y + dt*Vy + 1/2*(accy*(dt^2)) ;
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
    currx(i)= q*EC*partpertime;


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


      InBox1 = X > 0.8E-7 & X < 1.2E-7 & Y> 0.6E-7;
    outsidebox = Xp < 0.8E-7 | Xp > 1.2E-7;
    X(InBox1 & outsidebox) = Xp(InBox1 & outsidebox);
    Vx(InBox1 & outsidebox) = -Vx(InBox1 & outsidebox);
    Vy(InBox1 & ~outsidebox) = -Vy(InBox1 & ~outsidebox);
    Y(InBox1 & ~outsidebox) = Yp(InBox1 & ~outsidebox);

    InBox2 = X > 0.8E-7 & X < 1.2E-7 & Y < 0.4E-7;
    outsidebox = Xp < 0.8E-7 | Xp > 1.2E-7;
    X(InBox2 & outsidebox)= Xp(InBox2 & outsidebox);
    Vx(InBox2 & outsidebox) = -Vx(InBox2 & outsidebox);
    Vy(InBox2 & ~outsidebox) = -Vy(InBox2 & ~outsidebox);
    Y(InBox2 & ~outsidebox) = Yp(InBox2 & ~outsidebox);


    %mean free path
    pathnum = pathnum + sum(iscat);
    distance = sqrt((X1(iscat)-X(iscat)).^2 + ((Y1(iscat)-Y(iscat)).^2));
    distancesum = distancesum + sum(distance);
    X1(iscat) = X(iscat);


    %calculating temperarure using avg velocity

    avgV = sqrt(mean(Vx.^2 + Vy.^2));  %vx and vy are vectors
    semiT = (avgV).^2*mn/(2*k);

    figure(2);
    hold on
    axis([0 l 0 h])

    %titlep = sprintf('Temp = %.3d', semiT);

    title(sprintf('Particle Trajectories, Temp = %.3d', semiT))
    xlabel('x (m)')
    ylabel('y (m)')

    title('Curved Trajectories')
    for j=1:length(partWatch)

        plot([Xp(partWatch(j)),X(partWatch(j))]',[Yp(partWatch(j)),Y(partWatch(j))]','SeriesIndex',j)

    end

end



