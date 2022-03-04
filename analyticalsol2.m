function V = analyticalsol2(xInd,yInd,Vo,xMax,yMax,nMax)
%This is the function to calculate the analytical solution of the equation
%given in the assignment
%   We assign V1 to be the constant since it is not depedent on the n-term
%   of the equation, and it's defined by 4*Vo/pi. Then we take the sum in
%   respect of the odd n term, as the following equation:
%   (1/n)*(cosh(npix/a)/cosh(npib/a))*sin(npiy/a). a is defined by the
%   maximum width, b is equal to the maximum length divided by 2, x and y
%   defined by the index xInd and yInd. 
V1 = (4*Vo)/pi;
a=yMax;
b=xMax/2;
x=xInd-b;
y=yInd;
nList = 2*(0:round((nMax-1)*0.5))+1;
V = V1*sum((sin(nList*(pi*y/a))./(nList)).*(...
    cosh(nList*(pi*(x)/yMax))./cosh(nList*(pi*b/a))));
end