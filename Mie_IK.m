%% About this code
% Code for Mie scattering calculation of pristine spheres
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
r  = 2300;          % radius
WL = 400:0.1:460;   % wavelength range
n  = 2.4;           % real refractive index of sphere
k  = 0;             % imaginary refractive index of sphere
m0 = 1;             % Refractive index of surounding

% Outputs 
Csca = zeros(1,length(WL)); % Scattering cross section
Cext = Csca; % Extinction cross section
Cabs = Csca; % Absorption cross section

% Calculating the variation of cross-section coefficients for a given range of wavelength 
for i = 1:length(WL)
    [Csca(i), Cext(i), Cabs(i)] = Cal_Mie(WL(i), r, n, k, m0);
end

% Create plot
figure(1)

x0=10;
y0=10;
width=1200;
height=1000;
set(gcf,'position',[x0,y0,width,height])

hold on
plot(WL,Cext,'DisplayName',strcat('r=nm'),'LineWidth',2,'Color','b');
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
ax.XLim = [400,460];
%ax.YLim = [];

% Create title
title('Mie scattering in pristine spheres')
% Create xlabel
xlabel('wavelength (nm)', 'Interpreter', 'latex' );
% Create ylabel
ylabel('$\sigma_{abs}$', 'Interpreter', 'latex' );
