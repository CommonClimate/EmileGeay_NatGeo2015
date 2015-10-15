function [P,Err,Xhat,Yhat] = fitTLS(XData, YData)

% function [P,Err,Xhat,Yhat] = fitTLS(XData, YData)
%
% Orthogonal linear regression method in 2D for model: y = a + bx   
%
% Input parameters:
%  - XData: input data block -- x: axis
%  - YData: input data block -- y: axis
% Return parameters:
%  - Err: error - sum of orthogonal distances 
%  - P: vector of model parameters [b-slope, a-offset] 
%
% Authors:
% Ivo Petras (ivo.petras@tuke.sk) 
% Dagmar Bednarova (dagmar.bednarova@tuke.sk)
%
% adapted from fit_2D_data.m by J.E.G, USC, 22-Jan-2015

kx=length(XData);
ky=length(YData);;
if kx ~= ky
   disp('Incompatible X and Y data.');
   close all;
end
n=size(YData,2);
sy=sum(YData)./ky;
sx=sum(XData)./kx;
sxy=sum(XData.*YData);
sy2=sum(YData.^2);
sx2=sum(XData.^2);
B=0.5.*(((sy2-ky.*sy.^2)-(sx2-kx.*sx.^2))./(ky.*sx.*sy-sxy));
b1=-B+(B.^2+1).^0.5;
b2=-B-(B.^2+1).^0.5;
a1=sy-b1.*sx;
a2=sy-b2.*sx;
R=corrcoef(XData,YData);
if R(1,2) > 0 
    P=[b1 a1];
    Yhat = XData.*b1 + a1;
    Xhat = ((YData-a1)./b1);
end
if R(1,2) < 0
    P=[b2 a2];
    Yhat = XData.*b2 + a2;
    Xhat = ((YData-a2)./b2);
end   
alpha = atan(abs((Yhat-YData)./(Xhat-XData)));
d=abs(Xhat-XData).*sin(alpha);
%Err=sum(abs(d));
Err=sum(d.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%