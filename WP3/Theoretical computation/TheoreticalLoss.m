% Example of theoretical computation of complex permeability, specific loss, and loss separation
% in VITROPERM 500F 
% see: https://vacuumschmelze.com/Nanocrystalline-Material
%      https://vacuumschmelze.com/03_Documents/Brochures/EMC%20Products%20based%20on%20Nanocrystalline%20VITROPERM.pdf
clear
% Input
Jp = 5e-3; % chosen peak induction
load mu_vitroperm_500 VNA % experimental data are loaded 
f = VNA(:,1); % frequencies
mur_exp = VNA(:,2); % experimental permeability real part
mui_exp = VNA(:,3); % experimental permeability imaginary part
%
% Physical parameters
Is = 1.24; % the polarization at saturation
Ms = Is / (4e-7*pi); % the magnetization at saturation
chiDC = 98847; % the susceptibility
rho = 115e-8; % the resistivity
sigma = 1/rho; % the conductivity
% other physical parameters obtained by fitting of the experimental results.
d = 21.95e-6; %tickness of the ribbon
alpha = 0.0270; % damping constant
A = 1.8327e-11; % exchange constant
Nzz = 1; % demagnetizing factor along the direction perpendicular to the ribbon
Nyy = 0; % demagnetizing factor along the longitudinal (in plane) direction. Not used.
%
% Theoretical computation
ho = 1/chiDC;
[Loss, mu, Weddy,mul,ha]=Lossinribbon(alpha,Nyy,d,Ms,sigma,ho,Nzz,Jp,A,f);

% In the following the output is displayed
mur=real(mu); % Theoretical real permeability
mui=imag(mu); % Theoretical imaginary permeability
nu=1./mu; % Theoretical reluctivity
nu_exp=1./(mur_exp-1i*mui_exp); % Experimental reluctivity
figure
subplot(2,1,1)
semilogx(f,mur,'b-','LineWidth',1.5)
hold on
semilogx(f,mur_exp,'b.', 'MarkerSize',8)
semilogx(f,-mui,'r-','LineWidth',1.5)
semilogx(f,mui_exp,'r.','MarkerSize',8)
xlabel('Frequency (Hz)')
ylabel('Permeability')
legend('Theoretical permeability real part','Experimental permeability real part','Theoretical permeability imaginary part','Experimental permeability imaginary part')
title('VITROPERM 500F: permeability behaviour')
subplot(2,1,2)
loglog(f,imag(nu)*pi*Jp^2/(4e-7*pi),'b-','LineWidth',3)
hold on
loglog(f,imag(nu_exp)*pi*Jp^2/(4e-7*pi),'r.','MarkerSize',6)
loglog(f,imag(nu_exp)*pi*Jp^2/(4e-7*pi)- Weddy,'.k')
xlabel('Frequency (Hz)')
ylabel('Loss (J/m^3)')
legend({'Theoretical total loss','Experimental total loss','Excess loss component (total - eddy current loss)'}, 'Location','northwest')
title('Loss behavior at 5mT peak induction')