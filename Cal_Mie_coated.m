function[Csca, Cext, Cabs, x0] = Cal_Mie_coated(lambda, a, b, RI_re1, RI_im1,RI_re2, RI_im2, RI_re3, RI_im3, RI_b)

% The function generates the scattering, absorption and extinction coefficients for a given
% incident wavelength, when the boundary conditions of a coated core-shell sphere has been
% satisfied. 
%%
%  [Csca, Cext, Cabs, x0] = Cal_Mie_coated(lambda, a, b, RI_re1, RI_im1,RI_re2, RI_im2, RI_re3, RI_im3, RI_b)
%  
%  Description of inputs and outputs of the function Cal_Mie_coated
%  Input: 
%      lambda : Incident Wavelength.
%      a : Radius of inner sphere.
%      b : Radius of outer sphere.
%      RI_re1 : Real part of refractive index of inner sphere.
%      RI_im1 : Imaginary part of refractive index of inner sphere.
%      RI_re2 : Real part of refractive index of outer sphere.
%      RI_im2 : Imaginary part of refractive index of outer sphere.
%      RI_re3 : Real part of refractive index of outer medium.
%      RI_im3 : Imaginary part of refractive index of outer medium.
%      RI_b : Refractive index of surrounding medium.
%  Output:
%      Csca : Scattering Cross-section.
%      Cext : Extinction Cross-section.
%      Cabs : Absorption Cross-section.
%   Example:
%      [Csca, Cext, Cabs, x0] = Cal_Mie_coated(650, 300, 400, 2.03, 0 , 2.7, 0, 1, 0, 1)
%% Mie scattering in a coated sphere

n1=(RI_re1+1i*RI_im1)/RI_b;  % Effective refractive index of inner sphere
n2=(RI_re2+1i*RI_im2)/RI_b;  % Effective refractive index of outer sphere
n3=(RI_re3+1i*RI_im3)/RI_b;  % Effective refractive index of outer medium

index = RI_b; 
k0 = 2*pi*index/lambda; %progation constant

%Size parameters for inner and outer sphere
k1 = k0*n1;
k2 = k0*n2;
z0 = k0*b;
z1 = k2*a;
z2 = k2*b;
z3 = k1*a;
z = 1;
nmax = 3; % maximum number of modes
scaele = 0;
extele = 0;

% Boundary conditions for resonance to occur in coated sphere

for n = 1:nmax
    
    %Speherical Bessel functions
    jnz0 = sqrt(pi/(2*z0))*besselj(n+0.5,z0);
    jnminz0 = sqrt(pi/(2*z0))*besselj(n-0.5,z0);
    
    jnz1 = sqrt(pi/(2*z1))*besselj(n+0.5,z1);
    jnminz1 = sqrt(pi/(2*z1))*besselj(n-0.5,z1);
    
    jnz2 = sqrt(pi/(2*z2))*besselj(n+0.5,z2);
    jnminz2 = sqrt(pi/(2*z2))*besselj(n-0.5,z2);
    
    jnz3 = sqrt(pi/(2*z3))*besselj(n+0.5,z3);
    jnminz3 = sqrt(pi/(2*z3))*besselj(n-0.5,z3);
    
    %Speherical Hankel functions
    ynz1 = sqrt(pi/(2*z1))*bessely(n+0.5,z1);
    ynminz1 = sqrt(pi/(2*z1))*bessely(n-0.5,z1);
    
    ynz2 = sqrt(pi/(2*z2))*bessely(n+0.5,z2);
    ynminz2 = sqrt(pi/(2*z2))*bessely(n-0.5,z2);
    
    %Speherical Neumann functions
    hnz0 = sqrt(pi/(2*z0))*besselh(n+0.5,z0);
    hnminz0 = sqrt(pi/(2*z0))*besselh(n-0.5,z0);
    
    %Differentiation of spherical Bessel function
    jnz0diff = z0*jnminz0 - n*jnz0;
    jnz1diff = z1*jnminz1 - n*jnz1;
    jnz2diff = z2*jnminz2 - n*jnz2;
    jnz3diff = z3*jnminz3 - n*jnz3;
   
    %Differentiation of spherical Hankel function
    ynz1diff = z1*ynminz1 - n*ynz1;
    ynz2diff = z2*ynminz2 - n*ynz2;
   
    %Differentiation of spherical Neumann function
    hnz0diff = z0*hnminz0 - n*hnz0;
    
    %Assigning variables
    A1 = (jnz3diff./z3).*((ynz1.*jnz2)-(jnz1.*ynz2));
    B1 = jnz3.*(((ynz1diff./z1).*jnz2)-((jnz1diff./z1).*ynz2));
    C1 = jnz3.*(((ynz1diff./z1).*(jnz2diff./z2))-((jnz1diff./z1).*(ynz2diff./z2)));
    D1 = (jnz3diff./z3).*((ynz2.*(jnz2diff./z2))-(jnz1.*(ynz2diff./z2)));
 
    x1 = n1.*A1-n2.*B1;
    x2 = n2.*A1-n1.*B1;
    x3 = n2.*C1-n1.*D1;
    x4 = n1.*C1-n2.*D1;
    
  
    %Scattering coeffcients
    an = (index.*z.*((jnz0diff./z0)*x1)+((z*jnz0)+1i*index*(jnz0diff./z0))*n2*x3)./(index.*z.*((hnz0diff/z0)*x1)+(z*hnz0+1i*index*(hnz0diff/z0))*n2*x3);
    bn = ((index.*z.*jnz0*x4)+(z*(jnz0diff./z0)-1i*index*jnz0)*n2*x2)./((index.*z.*hnz0*x4)+(z*(hnz0diff/z0)-1i*index*hnz0)*n2*x2);
   
    scaele = scaele + (2*n+1)*(an*conj(an)+bn*conj(bn));
    extele = extele + (2*n+1)*real(an+bn);
    x0 = index.*z.*((hnz0diff/z0)*x1)+(z*hnz0+1i*index*(hnz0diff/z0))*n2*x3;
    % x0=(m.^2*jnmx*xhnxdiff-hnx*mxjnmxdiff);
    %x0(n) = (cn*m3xjnm3x+m3xhnm3x)*(m3*m2xjnm2xdiff)-(m2*m2xhnm2x)*(cn*m3xjnm3xdiff+m3xhnm3xdiff);
end

Csca = 2*pi/(k0.^2)*scaele;
Cext = 2*pi/(k0.^2)*extele;
Cabs = Cext-Csca;
