%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the viscosity from rotational diffusion for different sized %%
%% anisotropic nanoparticles. The expression of the friction factor depends 
%% on the exact shape of the spheroid. Reference: Martchenko et al. 2011 %%
%% Wang et al. 2014 (Using the discrete dipole approximation and         %%
%% holographic microscopy to measure rotational dynamics of non-spherical%%
%% colloidal particles)                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nr] = getViscosity(D,r, partType, T)
    p = r(1)/r(2);
    L = r(1)*10^-9;
    switch partType
        case 'Bipyramid'
             W = ((p^4 -1)/(p^4))/((((2*p^2 - 1)*log(p + sqrt(p^2 -1)))/(p*sqrt(p^2 -1)))-1);
        case 'Rod'
             a = [13.04468, -62.6084, 174.0921, -218.8365, 140.26992, -33.27076];
             for i = 1:6
                 factor(i) = a(i)*(p^(i/4));
             end
             W = (pi*D*L^3)/(3*(log(p) + 2*log(2) - 11/6 + (log(2)/log(1*p))*((1/3) - 2*log(2) + 11/6 - sum(a)) + sum(factor)));  

        otherwise
            error('Unknown dim, dim needs to be provided as 1D 2D or 3D')
    end

    nr = (1.380649*10^(-23)*T)/((2/3)*pi*D*L^3*W)*1000;%in cp
end
