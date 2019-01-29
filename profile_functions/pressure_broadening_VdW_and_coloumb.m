function pressure_broadened_profile = pressure_broadening_VdW_and_coloumb(omega,c,c2,a_p ,omega_peak, gamma, sigma, coloumb_width)
omega=omega(:);

%We calculate the profile over a bigger range and then return only the
%points corresponding to the omega given as input. This is to avoid edge
%effects on the returned profile caused by using fft and ifft to perform
%the convolution
omega_extension_width=3e3;
omega_lower_extension=[omega(1)+omega_extension_width:omega(2)-omega(1):omega(1)-(omega(2)-omega(1))]';
omega_upper_extension=[omega(end)+(omega(2)-omega(1)):omega(2)-omega(1):omega(end)-omega_extension_width]';

omega_tmp=[omega_lower_extension;omega;omega_upper_extension];

%Create the pure pressure broadened profile
dist=makedist('stable','alpha',0.5,'beta',1,'gam',c,'delta',c);
dist2=makedist('stable','alpha',0.5,'beta',-1,'gam',c2,'delta',-c2);

pure_pressure_broadening=pdf(dist,omega_tmp-omega_peak)+a_p*pdf(dist2,omega_tmp-omega_peak);

%Create the pressure broadened profile due to coloumb interactions
coloumb_profile=coloumb_broadened_profile(omega_tmp,omega_peak,coloumb_width);

%Create the voight profile centered at omega_peak
voigt_profile=myvoigt(omega_tmp,omega_peak,gamma,sigma);

%For some reason the convolution calculated by fft and ifft is rotated.
%Find out how rotated it is by calculating a convolution with a delta (ish)
%function and rotate back. Convolving with a delta function should not change the initial profile and therefore the 
%convolution should have its peak at the exact same point as before the convolution
[~,peak_index0]=findpeaks(voigt_profile);
if(isempty(peak_index0))
    [~,peak_index0]=findpeaks(circshift(voigt_profile,ceil(length(voigt_profile)/2)));
    peak_index0=peak_index0-ceil(length(voigt_profile)/2);
    
end

delta=zeros(size(omega_tmp));
delta(peak_index0)=1;
%delta(1)=1;
[~,peak_index1]=findpeaks(abs(ifft(fft(delta).*fft(voigt_profile))));
if(isempty(peak_index1))
    [~,peak_index1]=findpeaks(circshift(abs(ifft(fft(delta).*fft(voigt_profile))),round(length(voigt_profile)/2)));
    peak_index1=peak_index1-round(length(voigt_profile)/2);
end

pressure_broadened_profile=ifft(fft(pure_pressure_broadening).*fft(voigt_profile));
%pressure_broadened_profile=circshift(pressure_broadened_profile,peak_index0-peak_index1)/max(pressure_broadened_profile);
%For some reason the convolution calculated by fft and ifft is rotated.
%Find out how rotated it is by calculating a convolution with a delta (ish)
%function and rotate back. Convolving with a delta function should not change the initial profile and therefore the 
%convolution should have its peak at the exact same point as before the convolution
[~,peak_index0]=findpeaks(pressure_broadened_profile);
if(isempty(peak_index0))
    [~,peak_index0]=findpeaks(circshift(pressure_broadened_profile,ceil(length(pressure_broadened_profile)/2)));
    peak_index0=peak_index0-ceil(length(pressure_broadened_profile)/2);
    
end

delta=zeros(size(omega_tmp));
delta(peak_index0)=1;
%delta(1)=1;
[~,peak_index1]=findpeaks(abs(ifft(fft(delta).*fft(pressure_broadened_profile))));
if(isempty(peak_index1))
    [~,peak_index1]=findpeaks(circshift(abs(ifft(fft(delta).*fft(pressure_broadened_profile))),round(length(pressure_broadened_profile)/2)));
    peak_index1=peak_index1-round(length(pressure_broadened_profile)/2);
end

pressure_broadened_profile=ifft(fft(pressure_broadened_profile).*fft(coloumb_profile));
pressure_broadened_profile=circshift(pressure_broadened_profile,peak_index0-peak_index1)/max(pressure_broadened_profile);
pressure_broadened_profile=pressure_broadened_profile(length(omega_lower_extension)+1:end-length(omega_upper_extension))/max(pressure_broadened_profile);
pressure_broadened_profile(isnan(pressure_broadened_profile))=inf;
end