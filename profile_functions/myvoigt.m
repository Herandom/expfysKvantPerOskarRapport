function vp = myvoigt(xx, peakVal, gam, sig)
% MYVOIGT returns Voigt profile data.
% xx            Vector of x values.
% peakVal       Peak center value.
% gammaVal      Gamma value.
% sig           Sigma value.
z = arrayfun(@(q) ((q-peakVal)+1i*gam)/(sqrt(2)*sig), xx);
z_peak = arrayfun(@(q) ((q-peakVal)+1i*gam)/(sqrt(2)*sig), peakVal);
%vp = (1/(sig*sqrt(2*pi))) * real(Faddeeva_w(z)); % Get Voigt from Faddeeva fn.
vp = (1/(sig*sqrt(2*pi))) * real(fadf(z)); % Get Voigt from Faddeeva fn.
vp_peak=(1/(sig*sqrt(2*pi))) * real(fadf(z_peak));
vp = vp./max(vp_peak);
    if(gam<0)
        vp=1e38;
    end
end