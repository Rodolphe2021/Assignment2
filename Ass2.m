% Rodolphe Rejouis-101089193
% Assignement 2

%clear all
clearvars
clearvars -GLOBAL
close all
format shortE

%% Parameters for question 1
size = 20; 
W = 2*size;
L = 3*size;
Vo=1;
G1 = sparse(1,L); % Vector for question 1a
G2 = sparse(W,L); % Matrix for question 1b.
G3 = sparse(W,L); % Matrix to solve question 1 b analytical
BC1 = 0; 
BC2 = Vo;
iterations =100; 

%% Part a of question 1:
f1 = figure; 
G1(L)=BC1; % Set the boundary condition for the Vector to be V=0 at x=L 
G1(1)=BC2; % Set the boundary condition for the Vector to be V=Vo at x=0
G1(2:L-1)=deal(0.5*(BC1+BC2)); %populate the rest of the vector with the average of the Boundary Conditions
figure(f1)
for steps =1:iterations
    plot(G1);
    title('1-D finite iteration process')
    pause(0.01);
    for i=2:L-1
        G1(i)= 0.5*(G1(i-1)+G1(i+1));
    end
end

%% Part b of question 1
%Analytical solution for the equation :
% $$V(x,y)=\frac{4Vo}{\pi}\sum_{n=1,3,5...}^\infty\frac{1}{n}\frac{\cosh\left(\frac{n\pi x}{a}\right)}{\cosh\left(\frac{n\pi b}{a}\right)}\sin\left(\frac{n\pi y}{a}\right)$$

Nterms = iterations;
f2 = figure;
figure(f2)
for i = 1:L       
    for j = 1:W
        G3(j,i) = analyticalsol2(i,j,Vo,L,W,Nterms); 
    end
end
subplot(3,1,1);
surf(G3);
title('analytical solution');
pause(0.01);

% V=Vo at x=0,x=L and V=0 at y=0, y=W
[G2(1,1),G2(1,L),G2(W,1),G2(W,L),G2(2:W-1,2:L-1)]=deal( 0.5*(BC1+BC2) );
[G2(1,2:L-1),G2(W,2:L-1)]=deal(BC1);
[G2(2:W-1,1),G2(2:W-1,L)]=deal(BC2);

for step = 1:iterations
    subplot(3,1,2);
    surf(G2);
    title('numerical solution');
    subplot(3,1,3);
    surf(log(abs(G2-G3)))
    set(gca,'ZScale','log')
    title('difference between analytical and numerical')
    pause(0.01);
    for i=2:L-1
        for j=2:W-1
            G2(j,i) = 0.25*(G2(j+1,i)+G2(j-1,i)+G2(j,i+1)+G2(j,i-1));
        end
    end
end



%% Question 2 parameters:
% The length of the bottle neck equal to 1/3 of the total length of the rectangular
% region and the width is about 2/8 of the total width of the
% rectangular region. The conductivity inside the boxes is sig = 10e-2.

cMap2=ones(L,W); % Conductivity mapping
LGap = L/3; %For the length of the bottle neck
WGap =3*W/8; % Width of the boxes enclosing the bottle neck
sig = 10e-2;
[cMap2(floor((L-LGap)/2):floor((L+LGap)/2),1:floor((W-WGap)/2)),...
    cMap2(floor((L-LGap)/2):floor((L+LGap)/2),floor((W+WGap)/2):W)]= deal(sig);


G4 = sparse(L*W,L*W); 
B3 = zeros(1,L*W);

[B3(2:W-1),B3((L-1)*W+1:L*W-1)]=deal(BC2);
[B3(1),B3(W),B3((L-1)*W),B3(L*W)]=deal(0.5*(BC1+BC2)); 
for i=1:L
    for j=1:W
        n= j + (i-1) * W;   
        if i==1 || i == L || j == 1 || j == W
            G4(n,n) = 1;
        else
            nxp = j+i*W;
            nxm = j+(i-2)*W;
            nyp = n+1;
            nym = n-1;
            
            rxp = 0.5* (cMap2(i,j)+cMap2(i+1,j));
            rxm = 0.5* (cMap2(i,j)+cMap2(i-1,j));
            ryp = 0.5* (cMap2(i,j)+cMap2(i,j+1));
            rym = 0.5* (cMap2(i,j)+cMap2(i,j-1));
            
            G4(n,n) = -(rxp+rxm+ryp+rym);
            G4(n,nxp)=rxp;
            G4(n,nxm)=rxm;
            G4(n,nyp)=ryp;
            G4(n,nym)=rym;
        end    
    end
end

V2Vec=G4\(B3'); % The matrix neccessary to obtained our V(x,y) from the conductivity map

V2Map= zeros(L,W);
for i=1:L
    for j=1:W
        V2Map(i,j)=V2Vec(j+(i-1)*W);
    end
end
%% Process to obtain the Electric Field values for the vector format of the Electric Field 

eX= zeros(L,W);
eY= zeros(L,W);
for i=1:L
    for j=1:W
        if i == 1
            eX(i,j) = V2Map(i+1,j)-V2Map(i,j);
        elseif i == L
            eX(i,j) = V2Map(i,j)-V2Map(i-1,j);
        else
            eX(i,j) = (V2Map(i+1,j)-V2Map(i-1,j))*0.5;
        end
        if j == 1
            eY(i,j) = V2Map(i,j+1)-V2Map(i,j);
        elseif j == W
            eY(i,j) = V2Map(i,j)-V2Map(i,j-1);
        else
            eY(i,j) = (V2Map(i,j+1)-V2Map(i,j-1))*0.5;
        end
    end
end

eX=-eX;
eY=-eY;
% Code responsible for creating the Current Density Vectors
% The following code lines are responsible for making the Current Density
% Vector, by first assigning the value of the current Density to their
% matching position with the Electric Field (since J = sig*E). 
Jx=cMap2.*eX;% This line of code create the X value of the Current Density vector
Jy=cMap2.*eY;% This line of code create the Y value of the Current Density vector

pointSampling = 3;
[x,y]= meshgrid(1:floor(W/pointSampling),1:floor(L/pointSampling));
x=pointSampling*x-pointSampling+1;
y=pointSampling*y-pointSampling+1;

eXSample = zeros(floor(L/pointSampling),floor(W/pointSampling));
eYSample = zeros(floor(L/pointSampling),floor(W/pointSampling));
JxSample = zeros(floor(L/pointSampling),floor(W/pointSampling));
JySample = zeros(floor(L/pointSampling),floor(W/pointSampling));

for i=1:floor(L/pointSampling)
    for j = 1:floor(W/pointSampling)
        eXSample(i,j)=eX(pointSampling*i,pointSampling*j);
        eYSample(i,j)=eY(pointSampling*i,pointSampling*j);
        JxSample(i,j)=Jx(pointSampling*i,pointSampling*j);
        JySample(i,j)=Jy(pointSampling*i,pointSampling*j);
    end
end
%% Plot of the question 2.a
% This is plots for the conductivity sig(x,y), V(x,y), E(x,y), and J(x,y).
% This will be standard plot that will compare to the following figures,
% where the setting are modified to investigate the changes.
f3=figure;
figure(f3)
subplot(2,2,1)
surface(cMap2');
title('Conductivity');
subplot(2,2,2)
surface(V2Map');
title('Voltage Vx,Vy');
subplot(2,2,3)
quiver(y,x,eXSample,eYSample,10);
title('Electric Field, Ex,Ey');
subplot(2,2,4)
quiver(y,x,JxSample,JySample,10);
title('Current Density, Jx,Jy');

%% Answer of Question 2.b

% This is with meshing size being double of the standard set:
Plotconductivity(40,10e-2,3);

% This is with meshing size being triple of the standard set:
Plotconductivity(80,10e-2,3);

% This is with meshing size being half of the standard set:
Plotconductivity(10,10e-2,3);

%% Answer of Question 2.c:
% Here are the different version of the a set with different bottle neck size
%
% Compariso:
Plotconductivity(100,10e-2,3);
% This is one with longer bottle neck:
Plotconductivity2(100,10e-2,3,0.5,0.375);
% This is one with narrower bottle neck:
Plotconductivity2(100,10e-2,3,0.3333,0.5);

%% Answer of Question 2.d:
% Comparing the results with different resistivity inside the box:

% Higher resistivity
Plotconductivity(20,10,3);

% Lower resistivity
Plotconductivity(20,10e-5,3);