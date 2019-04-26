function [R,G,L,C,Z0] = trasmpar(epsR,muR,tand,ro,b,a,f)
%Calculation of transmission line parameters
%   Function takes relative permittivity, relative permeability, loss tangent, resistivity, outer conductor
%   diameter, inner conductor diameter, frequency to calculate transmission line
%   distributed parameters
eps0=8.854187817*1e-12;
eps=eps0*epsR*(1-1i*tand);
mu0=4*pi*1e-7;
mu=mu0*muR;
sk=sqrt(pi*f*ro*mu);
R=sk/(pi*b)+sk/(pi*a);    %[Ohm/m]
G=4*pi^2*f*imag(eps)/log(b/a);  %[S/m]
L=mu/(2*pi)*(log(b/a));         %[H/m]
C=2*pi*eps0*epsR/log(b/a);      %[F/m]
Z0=60*log(b/a)/(epsR);          %[Ohm]
end

