function coloumb_broadened_profile=coloumb_broadened_profile(omega,omega_peak,omega_width)
omega=omega(:);
x=linspace(0,20,2000);
dx=x(1)-x(2);
coloumb_broadened_profile=sum(cos((omega-omega_peak)/omega_width*x).*exp(-(repmat(x,size(omega)).^3).^(1/2)),2)*dx;
% f_tmp=@(x)cos((omega-omega_peak)/omega_width*x)*exp(-(x^3).^(1/2));
% coloumb_broadened_profile=integral(f_tmp,0,inf,'ArrayValued',true);
end