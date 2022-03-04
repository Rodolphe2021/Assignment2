function [] = Plotconductivity2(size,sig,ps,L1,W1)
%This is condensed code for question 2
%   This code is a repeat of what was used to solve question 2.a) But this
%   time there will be some inputs to allow for change to occur
BC1 = 0;
BC2 = 1;
L = 3*size; 
W = 2*size;
cMap2=ones(L,W);
LGap = L*L1;
WGap = W*W1;

[cMap2(floor((L-LGap)/2):floor((L+LGap)/2),1:floor((W-WGap)/2)),...
    cMap2(floor((L-LGap)/2):floor((L+LGap)/2),floor((W+WGap)/2):W)]= deal(sig);


G4 = sparse(L*W,L*W);
B3 = zeros(1,L*W);
%equations will be 
%4T_{i,j} -T_{i-1,j} - T_{i+1,j}-T_{i,j-1} - T_{i,j+1}= 0
%T_{1,j} = V0=BC2
%T_{L,j} = V0=BC2
%T_{i,1} = 0=BC1
%T_{i,W} = 0=BC1
%T_{1,1}=T_{L,1}=T_{1,W}=T_{L,W}=V0/2=(BC1+BC2)/2


[B3(2:W-1),B3((L-1)*W+1:L*W-1)]=deal(BC2);
[B3(1),B3(W),B3((L-1)*W),B3(L*W)]=deal(0.5*(BC1+BC2));


for i=1:L
    for j=1:W
        n= j + (i-1) * W;   
        if i==1 | i == L | j == 1 | j == W
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

V2Vec=G4\(B3');

V2Map= zeros(L,W);
for i=1:L
    for j=1:W
        V2Map(i,j)=V2Vec(j+(i-1)*W);
    end
end

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

Jx=cMap2.*eX;
Jy=cMap2.*eY;

pointSampling = ps;
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



f3=figure;
figure(f3)
subplot(1,2,1)
surface(V2Map');
title('Voltage Vx,Vy');
subplot(1,2,2)
quiver(y,x,JxSample,JySample,10);
title('Current Density, Jx,Jy');
end