function[Csca, Cext, Cabs, x0] = Cal_Mie(lambda, radius, RI_re, RI_im, RI_b)

% The function generates the scattering, absorption and extinction coefficients for a given
% incident wavelength, when the boundary conditions of a pristine sphere has been
% satisfied. 
%%
%  [Csca, Cext, Cabs, x0] = Cal_Mie(lambda, radius, RI_re, RI_im, RI_b)
%
%  Description of inputs and outputs of the function Cal_Mie
%  Input: 
%      lambda : Incident Wavelength.
%      radius : Radius of sphere.
%      RI_re : Real part of refractive index of sphere.
%      RI_im : Imaginary part of refractive index of sphere.
%      RI_b : Refractive index of surrounding medium.
%  Output:
%      Csca : Scattering Cross-section.
%      Cext : Extinction Cross-section.
%      Cabs : Absorption Cross-section.
%   Example:
%      [Csca, Cext, Cabs, x0] = Cal_Mie(650,600,2.03,0,1)
%% Mie scattering in a pristine sphere

m = (RI_re+1i*RI_im)/RI_b; % Effective refractive index
index = RI_b;
k = 2*pi*index./lambda; % Propagartion constant
x = k.*radius; % Size parameter
mx = m.*x;
nmax = 100; % Maximum number of modes             
scaele = 0;
extele = 0;

% Boundary conditions for resonance to occur in pristine spheres

for n = 1:nmax

    % Spherical Bessel function 
    jnx = sqrt(pi./(2.*x)).*besselj(n+0.5,x); 
    jnminx = sqrt(pi./(2.*x)).*besselj(n-0.5,x);
    jnmx = sqrt(pi./(2.*mx)).*besselj(n+0.5,mx);
    jnminmx = sqrt(pi./(2.*mx)).*besselj(n-0.5,mx);

    % Spherical Neumann function
    hnx = sqrt(pi./(2.*x)).*besselh(n+0.5,x);  
    hnminx = sqrt(pi./(2.*x)).*besselh(n-0.5,x);

    %Differentiation of spherical bessel function
    xjnxdiff = x.*jnminx - n.*jnx; 
    mxjnmxdiff = mx.*jnminmx - n.*jnmx;

    %Differentiation of spherical Neumann function
    xhnxdiff = x.*hnminx - n.*hnx;

    %Scattering coefficients
    an =(m.^2.*jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(m.^2.*jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    bn =(jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(jnmx.*xhnxdiff-hnx.*mxjnmxdiff);


    scaele = scaele + (2*n+1).*(an.*conj(an)+bn.*conj(bn));
    extele = extele + (2*n+1).*real(an+bn);


    x0=(m.^2.*jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
end


Csca = 2*pi./(k.^2).*scaele;
Cext = 2*pi./(k.^2).*extele;
Cabs = Cext-Csca;