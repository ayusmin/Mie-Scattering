%% About this code
% Code for Mie scattering calculation of coted spheres
% author: Ayusmin Panda and B. R. K. Nanda
% Date: 15/07/2024
% doi: https://doi.org/10.1002/adpr.202300339
% contact: pandaayusmin97@gmail.com and nandab@iitm.ac.in
% for more information visit https://www.cmtcl-iitm.com
%%
clear
clc
close all

% Inputs
WL = 400:0.1:750;  %Wavelength
r1 = 300;   %Radius of inner sphere
r2 = 400;   %Radius of outer sphere
n1 = 2.8;   %Real part of refractive index of inner sphere
n2 = 1.2;   %Real part of refractive index of outer sphere
n3 = 1;     %Real part of refractive index of outer mediun
k1 = 0;     %Imaginary part of refractive index of inner sphere
k2 = 0;     %Imaginary part of refractive index of outer sphere
k3 = 0;     %Imaginary part of refractive index of outer medium
m0 = 1;     %Refractive index of surounding

%Outputs
Csca = zeros(1,length(WL));%Scattering cross section
Cext = Csca; %Extinction cross section
Cabs = Csca; %absorption cross section

% Calculating the variation of cross-section coefficients for a given range of wavelength 
for i=1:length(WL)
    [Csca(i), Cext(i), Cabs(i), x0] = Cal_Mie_coated(WL(i), r1,r2, n1, k1, n2, k2, n3, k3, m0);
    %[Csca, Cext, Cabs, x0] = Cal_Mie_coated(lambda, a, b, RI_re1, RI_im1,RI_re2, RI_im2, RI_b)
end


% Create plot
figure(1)
x0 = 10;
y0 = 10;
width = 1200;
height = 1000;
set(gcf,'position',[x0,y0,width,height])

hold on
plot(WL,Csca,'DisplayName',strcat('r=nm'),'LineWidth',2,'Color','b');
hold off

% Create axes
ax = gca;
% Set the axes properties
ax.Box = 'on';
grid on
ax.GridLineStyle = '-';
ax.GridLineWidth = 1;
ax.LineWidth = 2;
ax.FontSize = 22;
ax.XLim = [400,750];
%ax.YLim = [];

% Create title
title('Mie Scattering in coated spheres')
% Create xlabel
xlabel('wavelength (nm)', 'Interpreter', 'latex' );
% Create ylabel
ylabel('$\sigma_{abs}$', 'Interpreter', 'latex' );
