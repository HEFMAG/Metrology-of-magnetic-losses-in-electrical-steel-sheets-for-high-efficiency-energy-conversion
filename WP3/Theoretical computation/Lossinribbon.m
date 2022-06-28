function [Loss,mu,Weddy,mul,ha]=Lossinribbon(alpha,Nyy,d,Ms,sigma,hoo,Nzz,Jp,A,f)
% variables in input: alpha is the damping constant; Nyy is the demagnetizing factor along the
% longitudinal direction of the ribbon; d is the tickness of the ribbon;
% Ms is the saturation magnetization of the given magnetic material (Ampere
% meter^-1); sigma is the conductivity (Siemens meter^-1); 
% hoo is the inverse of the magnetic susceptibility [non dimensional] (https://www.electropedia.org/iev/iev.nsf/display?openform&ievref=121-12-37)
% Nzz is the demagnetizing factor in the direction perpendicular to the ribbon; 
% A is the exchange constant (Joule meter^-1); Jp is the peak
% induction(s) at which the loss is computed; f is the frequency(es).
%
% Variables in output:
% Loss is the specific total loss per cycle [Joule meter^-3]; mu is the complex permeability; 
% Weddy is the specific loss per cycle [Joule meter^-3] due to eddy
% currents;mul is the complex local permeability; ha is the complex applied
% field.

mu0=4e-7*pi; %vacuum permeability
gamma = 1.761e11*mu0; % electron gyromagnetic ratio x mu0 [seconds^-1 Ampere^-1 meter]
l=sqrt(2*A/(mu0*Ms^2)); % exchange length [meter]
jp=Jp/Ms; %reduced induction
omega=2*pi*f/(gamma*Ms); % reduced angular frequency
ho=hoo+1i*omega*alpha;
mul=mu0*(1+1./(ho+Nyy-omega.^2./(ho+Nzz)));
chi=1./(ho+Nyy-omega.^2./(ho+Nzz));
omega=omega*gamma*Ms; % absolute angular frequency
% the solution of the diffusion equation start here
a1=-(1./(l^2*chi)+1i.*omega*sigma*mu0);
a2=1i.*omega*sigma.*mul./(l^2*chi);
s1=sqrt(-a1/2+sqrt(a1.^2/4-a2));
s2=sqrt(-a1/2-sqrt(a1.^2/4-a2));

dh1=1i.*omega.*sigma*d/2.*jp;
dh2=-omega.^2.*sigma^2*d/2*mu0.*jp;

K1=1./(s1.*s2.*(s2.^2-s1.^2)).*((s2.^3).*dh1+(-s2).*dh2);
K2=1./(s1.*s2.*(s2.^2-s1.^2)).*((-s1.^3).*dh1+(s1).*dh2);

% Wdamp=0.5*Ms^2.*2*pi.*1i*(mul).*2/d.* ...
% 0.5.*(abs(K1).^2./abs(sinh(s1*d/2)).^2.*(1./(s1+conj(s1)).*sinh((s1+conj(s1))*d/2)+1./(s1-conj(s1)).*sinh((s1-conj(s1))*d/2))+ ...
% +(K1).*conj(K2)./(sinh(s1*d/2).*sinh(conj(s2)*d/2)).*(1./(s1+conj(s2)).*sinh((s1+conj(s2))*d/2)+1./(s1-conj(s2)).*sinh((s1-conj(s2))*d/2))+ ...
% +(K2).*conj(K1)./(sinh(s2*d/2).*sinh(conj(s1)*d/2)).*(1./(s2+conj(s1)).*sinh((s2+conj(s1))*d/2)+1./(s2-conj(s1)).*sinh((s2-conj(s1))*d/2))+ ...
% +abs(K2).^2./abs(sinh(s2*d/2)).^2.*(1./(s2+conj(s2)).*sinh((s2+conj(s2))*d/2)+1./(s2-conj(s2)).*sinh((s2-conj(s2))*d/2)));
% Wdamp=real(Wdamp);
%

Weddy=0.5/sigma*Ms^2.*2/d.*0.5.* ...
(abs(K1).^2.*abs(s1).^2./abs(sinh(s1*d/2)).^2.*(1./(s1+conj(s1)).*sinh((s1+conj(s1))*d/2)-1./(s1-conj(s1)).*sinh((s1-conj(s1))*d/2))+ ...
+(K1).*conj(K2).*s1.*conj(s2)./(sinh(s1*d/2).*sinh(conj(s2)*d/2)).*(1./(s1+conj(s2)).*sinh((s1+conj(s2))*d/2)-1./(s1-conj(s2)).*sinh((s1-conj(s2))*d/2))+ ...
+(K2).*conj(K1).*s2.*conj(s1)./(sinh(s2*d/2).*sinh(conj(s1)*d/2)).*(1./(s2+conj(s1)).*sinh((s2+conj(s1))*d/2)-1./(s2-conj(s1)).*sinh((s2-conj(s1))*d/2))+ ...
+abs(K2).^2.*abs(s2).^2./abs(sinh(s2*d/2)).^2.*(1./(s2+conj(s2)).*sinh((s2+conj(s2))*d/2)-1./(s2-conj(s2)).*sinh((s2-conj(s2))*d/2)))./f;
Weddy=real(Weddy);
%

ha=1./(s1.*s2.*(s2.^2-s1.^2)).*((s2.^3.*coth(s1*d/2)-s1.^3.*coth(s2*d/2)).*dh1+(-s2.*coth(s1*d/2)+s1.*coth(s2*d/2)).*dh2);

have=1./(s1.*s2.*(s2.^2-s1.^2)).*((s2.^3./(s1*d/2)-s1.^3./(s2*d/2)).*dh1+(-s2./(s1*d/2)+s1./(s2*d/2)).*dh2);
%S=0.5*Ms^2*1i.*omega.*jp.*conj(ha);
s1=conj(s1);
s2=conj(s2);

S1=pi^2/6*sigma*d^2*jp.^2*Ms^2.*f.^2.*3./(s2.^2-s1.^2).*((s2.^2.*coth(s1*d/2)./(s1*d/2)-s1.^2.*coth(s2*d/2)./(s2*d/2))+1i.*omega*sigma*mu0.*(coth(s1*d/2)./(s1*d/2)-coth(s2*d/2)./(s2*d/2)));

%Wclc=pi^2/6*sigma*d^2*(Jp).^2.*f;
mu=-1i*omega./(2*conj(S1))*Ms^2.*jp.^2;
mu=mu/mu0;
mul=mul/mu0;
Loss=real(S1)./f;
%Wdamp=Loss-Weddy;
end
