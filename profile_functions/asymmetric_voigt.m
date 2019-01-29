function [ asymmetric_voigt ] = asymmetric_voigt( x,x_peak,gamma1,gamma2,sigma1,sigma2 )

asymmetric_voigt=(x>=x_peak).*myvoigt(x,x_peak,gamma1,sigma1)+(x<x_peak).*myvoigt(x,x_peak,gamma2,sigma2);
if(any(isnan(asymmetric_voigt)))
    keyboard
end
end

