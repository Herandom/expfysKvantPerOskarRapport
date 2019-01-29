%addpath kalibrering
%addpath grunduppgift
addpath heliumbös
addpath extra

close all
c=299792458;

%Settings
log_intensity=1;
lin_intensity=~log_intensity;

fit_baseline_to_full_spectrum=1;
include_respons_function_from_datasheet=0;

%Guess for the ratio of the gaussian width of the peaks and the lorentzian width of the peaks
%Usually between 1 and 5 but as big as 10 can be needed in some cases.
%Profiles with thin peaks and long tails corresponds to lower value of
%GLratio and fat profiles with not much tails correspond to higher values.
GLratio=1;
%Guess for the ratio of the voight width of the peaks and the pressure broadening width of the peaks
%Usually between 1 and 5 but as big as 10 can be needed in some cases, as with GLratio.
%Big skewness corresponds to lower values.
VPBratio=1;
%Guess for the ratio of the contribution from thw coloumb broadening and 
%the Van der Waals broadening to the pressure broadened profile
CVdWRatio=0.3;
%Guess for the ratio of the height of the right skewed and the left skewed
%part of the pressure broadening
PLRratio=3;

%The width gets a little too wide due to the sampling step size, compensate
%for this when guessing
width_guess_scale_factor=0.7;

%The estimate of the peak width parameters are the width where the peak has
%decreased a factor with width_definition_decrease relative to peak value
width_definition_decrease=0.6;

%The range around the peak included in the fit is the width where the peak has
%decreased a factor with range_definition_decrease relative to peak value
range_definition_decrease=1/100;

%Choose a model for the fitting
%1,2 or 3 gives symmetric voigt profile with 1,2 or 3 peaks respectively
%4 and 5 gives symmetric voigt profile with 1 and 2 peaks respectively, 
%devided by a 1st degree ploynomial that compensates for a linearized model
%of the spectrometer respons curve, and added to a constant that removes 
%the last of the baseline. Suitable for calibration with sodium duoble peak
%6 gives asymmetric voigt top, with different decay parameters on each side
%of the peak. 
%7 and 8 fits one and two peaks respectively to a model that
%convolves a voigt profile and a profile that modells pressure broadening
%by only taking into account interaction to the closest nearby atoms.
%9 and 10 fits one and two peaks respectively to a model that convolves a
%voigt profile with a general stable distribution used for modelling
%pressure broadening.
%11 and 12 convolves a voigt profile with two Levy profiles, representing
%transitions with an energy difference that is decreased and increased by 
%pressure broadening, respectively
%13 is the same as 11 but convolved with an aditional profile for
%broadening due to coloumb interactions with charged particles
model=11;

%CVdWRatio is not approperiate for other models than 13
if(model==13)
    CVdWRatio=0;
end

%Load data
spectrum=load('He4.txt');
lambda=spectrum(1,:);
omega=2*pi*c./lambda;

%Do we want to plot the intensity with a linear or logarithic scale?
if(lin_intensity)
    intensity=spectrum(2,:);
    
elseif(log_intensity)
    intensity=log(spectrum(2,:));
end

%Adjust calibration so that the calibration is made with the same model as
%the one the peaks are fitted to
calibrated_with_kal13=1;
if(calibrated_with_kal13)
    
    switch model
        case 4  % if you fit single peaks with model 4, the corresponding
            %calibration model should be 5 since model 5 is the same as
            %model 4 but fits 2 peaks instead of 1
            %peaks from kal13 with model=5
            sodium_peak1=588.976756236704;
            sodium_peak2=589.560943817678;
        case 7
            
            %peaks from kal13 with model=8
            sodium_peak1= 589.0458853188291  ;
            sodium_peak2=589.5966234463626 ;
            
        case 9
            %peaks from kal13 with model=10
            sodium_peak1= 589.0143770839022 ;
            sodium_peak2=589.5594530255535 ;
            
        case 11
            %peaks from kal13 with model=12
%             sodium_peak1= 589.0211596876717  ;old
%             sodium_peak2=589.6015693034709 ;
            
            sodium_peak1=589.0143747624688;
            sodium_peak2=589.5967923761523;
        case 13
              %peaks from kal13 with model=12
%             sodium_peak1= 589.0211596876717  ;old
%             sodium_peak2=589.6015693034709 ;
            
            sodium_peak1=589.0143747624688;
            sodium_peak2=589.5967923761523;
    end
    real_sodium_peak1=588.9950;%From NIST
    real_sodium_peak2=589.5924;%From NIST
    
    fine_calibration_lambda_shift_kal13=((real_sodium_peak1-sodium_peak1)+(real_sodium_peak2-sodium_peak2))/2;
    lambda=lambda+fine_calibration_lambda_shift_kal13;
end
calibrated_with_kal14=0;
if(calibrated_with_kal14)
    switch model
        case 7
            %peaks from kal14 with model=8
            sodium_peak1= 589.0458853188291  ;
            sodium_peak2=589.5966234463626 ;
        case 9
            %peaks from kal14 with model=10
            sodium_peak1= 589.0143813357297  ;
            sodium_peak2=589.5880025592129 ;
        case 11
            %peaks from kal14 with model=12
            sodium_peak1= 588.9906837408184  ;
            sodium_peak2=589.6081256471466 ;
        case 13
            %peaks from kal14 with model=12
            sodium_peak1= 588.9906837408184  ;
            sodium_peak2=589.6081256471466 ;            
    end
    real_sodium_peak1=588.9950;%From NIST
    real_sodium_peak2=589.5924;%From NIST
    
    fine_calibration_lambda_shift_kal14=((real_sodium_peak1-sodium_peak1)+(real_sodium_peak2-sodium_peak2))/2;
    lambda=lambda+fine_calibration_lambda_shift_kal14;
end
%Plot spectrum and pause to enable zooming
plot(lambda,intensity)
pause

%Take input for lower and upper wavelength range, and also for the peak
%threshold (saying that all peaks are over this threshold)
lambda_range_and_threshold=ginput(3);
lambda_range=[lambda_range_and_threshold(1,1) lambda_range_and_threshold(2,1)]
peak_threshold=lambda_range_and_threshold(3,2)

%From now on, only use the range specified by the user
intensity=intensity(lambda>lambda_range(1)&lambda<lambda_range(2));
lambda=lambda(lambda>lambda_range(1)&lambda<lambda_range(2));
%Plot the selected range
plot(lambda,intensity)

%Find peaks using findpeaks. Create a temporary intesity variable where we set all elements below
%peak_threshold to zero. This is needed to make sure that findpeaks does
%not find all the very smal "peaks" in the noise
intensity_tmpcpy=intensity;
if(log_intensity)
    intensity_tmpcpy(intensity_tmpcpy<peak_threshold)=min(intensity_tmpcpy);
elseif(lin_intensity)
    intensity_tmpcpy(intensity_tmpcpy<peak_threshold)=0;
end
[~,index_peaks]=findpeaks(intensity_tmpcpy);

%Initialize arrays where we store the wavelength indicies where the intensity has
%decreased by half the peak value (width_param_index) and to 1/20 of the
%peak value (lambda_intervall_index_length, obs!!! number of indicies from peak)
width_param_index=zeros(1,length(index_peaks));
width_param_index_upper=zeros(1,length(index_peaks));
width_param_index_lower=zeros(1,length(index_peaks));
lambda_intervall_index_length_m=zeros(1,length(index_peaks));
lambda_intervall_index_length_p=zeros(1,length(index_peaks));

%If requested, fit a 2nd degree polynomial baseline and subtract that baseline. We fit
%the baseline to the non-logarithmic intensity regardles of the value of
%log_intensity
if(fit_baseline_to_full_spectrum)
    baseline=polyfit(lambda(intensity<peak_threshold),exp(intensity(intensity<peak_threshold)),2);
    %intensity=exp(intensity)-(baseline(1)*lambda.^3+baseline(2)*lambda.^2+baseline(3)*lambda+baseline(4));
    intensity=exp(intensity)-(baseline(1)*lambda.^2+baseline(2)*lambda+baseline(3));
    intensity(intensity<0)=min(intensity(intensity>0));
    if(include_respons_function_from_datasheet)
        intensity=log(intensity./hamamatsu_R375_respons(lambda));
    else
        intensity=log(intensity);
    end
end

%Find index for the estimate of the width parameter and the range within 
%wich we want to make a fit to the peaks
for i=1:length(index_peaks)   
    if(lin_intensity)
        width_param_index_upper(i)=find(intensity(index_peaks(i):end)./intensity(index_peaks(i))<width_definition_decrease,1)+index_peaks(i);
        width_param_index_lower(i)=find(intensity(1:index_peaks(i))./intensity(index_peaks(i))<width_definition_decrease,1,'last');
        lambda_intervall_index_length_p(i)=find(intensity(index_peaks(i):end)-intensity(index_peaks(i))<log(range_definition_decrease),1);
        lambda_intervall_index_length_m(i)=find(fliplr(intensity(1:index_peaks(i)))-intensity(index_peaks(i))<log(range_definition_decrease),1);
    elseif(log_intensity)
        width_param_index_upper(i)=find(intensity(index_peaks(i):end)-intensity(index_peaks(i))<log(width_definition_decrease),1)+index_peaks(i);
        width_param_index_lower(i)=find(intensity(1:index_peaks(i))-intensity(index_peaks(i))<log(width_definition_decrease),1,'last');
        lambda_intervall_index_length_p(i)=find(intensity(index_peaks(i):end)-intensity(index_peaks(i))<log(range_definition_decrease),1);
        lambda_intervall_index_length_m(i)=find(fliplr(intensity(1:index_peaks(i)))-intensity(index_peaks(i))<log(range_definition_decrease),1);
    end
end
%Maximum locations from raw data
lambda_peaks=lambda(index_peaks);
omega_peaks=omega(index_peaks);

%Find estimate of half width expressed in angular frequency, since thats
%what the voigt profile is a function of
HWHF=1/2*(2*pi*c./lambda(width_param_index_lower)-2*pi*c./lambda(width_param_index_upper));

%Calculate guesses for gamma and sigma knowing that the full width is
%related to the Gaussin and Lorenzian width according to the formulas below
%(approximations taken from wikipedia)
gamma_guesses=HWHF/(0.5346+sqrt(0.2166+GLratio^2))*width_guess_scale_factor;
sigma_guesses=HWHF/(0.5346/GLratio+sqrt(0.2166/GLratio^2+1))/sqrt(2*log(2))*width_guess_scale_factor;

%Calculate guesses for the width of the pressure broadened profile(s) and
%change gamma_guesses and sigma_guesses accordingly, assuming the total
%width is the sum of sigma, gamma and c.
c_guesses=HWHF/(VPBratio+1)*width_guess_scale_factor*2/1.1*1/(1+CVdWRatio);%Calculated for model 7, seems to work ok for model 9 as well
gamma_guesses=gamma_guesses*VPBratio/(1+VPBratio);
sigma_guesses=sigma_guesses*VPBratio/(1+VPBratio);
coloumb_width_guesses=HWHF/(VPBratio+1)*width_guess_scale_factor*CVdWRatio/(1+CVdWRatio);


%Calculate guesses for the skewnes parameter b. Guess 0.7 if intensity is
%skewed to the left as a function of lambda and -0.7 if skewed to the right
%NOTE that the profile is a function of omega and that a left skewness in 
%lambda corresponds to a skewness to the right in omega and vice versa
b_guesses=-0.7+1.4*(intensity(width_param_index_upper)<intensity(2*index_peaks-width_param_index_upper));

%Plot the found peaks and the width estimates so the user can check that
%everything worked
plot(lambda,intensity)
hold on

plot(lambda(index_peaks),intensity(index_peaks),'*')
plot(lambda(width_param_index_upper),intensity(width_param_index_upper),'x')
plot(lambda(width_param_index_lower),intensity(width_param_index_lower),'x')
hold off
%%
%addpath Voightfit
%addpath Faddeeva
%addpath extra
addpath profile_functions

%Create a string defining the desired model and prealocate for the output and
%set fit options for the desired model. Hardcoded values for 'upper' and
%'lower' seems to make the fitting robust so far
    switch model
        case 1
            function_string=['A*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda,gamma,sigma)'];
            fitted_peaks=zeros(length(index_peaks),4);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 1e6 1e6 1e6],'lower',[0 0 0 0] );
        case 2
            function_string=['A*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda,gamma,sigma)+A2*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda2,gamma2,sigma2)'];
            %TODO include fitting a linearized respons curve and extra
            %baseline subtracting in this model
            fitted_peaks=zeros(length(index_peaks),8);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 100 1e6 1e6 1e6 1e6 1e6 1e6],'lower',[0 0 0 0 0 0 0 0] );
        case 3
            function_string=['A*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda,gamma,sigma)+A2*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda2,gamma2,sigma2)+A3*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda3,gamma3,sigma3)'];
            fitted_peaks=zeros(length(index_peaks),12);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 100 100 1e6 1e6 1e6 1e6 1e6 1e6 1e6 1e6],'lower',[0 0 0 0 0 0 0 0 0 0 0] );
        case 4
            function_string=['A*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda,gamma,sigma)./(1+t*x)+z'];
            fitted_peaks=zeros(length(index_peaks),6);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 1e6 1e6 1e6 1 1],'lower',[0 0 0 0 -1 -1] );
       case 5
           function_string=['(A*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda,gamma,sigma)+A2*myvoigt(2*pi*299792485./x,2*pi*299792485/lambda2,gamma2,sigma2))./(1+t*x)+z'];
            fitted_peaks=zeros(length(index_peaks),10);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 100 1e6 1e6 1e6 1e6 1e6 1e6 1 1],'lower',[0 0 0 0 0 0 0 0 -1 -1] );
        case 6
            function_string=['A*asymmetric_voigt(2*pi*299792485./x,2*pi*299792485/lambda,gamma1,gamma2,sigma1,sigma2)'];
            fitted_peaks=zeros(length(index_peaks),6);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1e6 1e6 1e6 1e6 1e6],'lower',[0 0 0 0 0] );
        case 7
            if(b_guesses(1)<0)
                function_string=['pressure_broadened_profile(2*pi*299792485./x,c,3,2*pi*299792485/lambda,gamma,sigma,-1)'];
            else
                function_string=['pressure_broadened_profile(2*pi*299792485./x,c,3,2*pi*299792485/lambda,gamma,sigma,1)'];
            end
            fitted_peaks=zeros(length(index_peaks),4);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[500 500 1e4 1e4],'lower',[0 0 0 0] );  
        case 8
            if(b_guesses(1)<0&&b_guesses(2)<0)
                function_string=['pressure_broadened_profile(2*pi*299792485./x,c,3,2*pi*299792485/lambda,gamma,sigma,-1)+A2*pressure_broadened_profile(2*pi*299792485./x,c2,3,2*pi*299792485/lambda2,gamma2,sigma2,-1)'];
            elseif(b_guesses(1)>0&&b_guesses(2)<0)
                function_string=['pressure_broadened_profile(2*pi*299792485./x,c,3,2*pi*299792485/lambda,gamma,sigma,1)+A2*pressure_broadened_profile(2*pi*299792485./x,c2,3,2*pi*299792485/lambda2,gamma2,sigma2,-1)'];
            elseif(b_guesses(1)<0&&b_guesses(2)>0)
                function_string=['pressure_broadened_profile(2*pi*299792485./x,c,3,2*pi*299792485/lambda,gamma,sigma,-1)+A2*pressure_broadened_profile(2*pi*299792485./x,c2,3,2*pi*299792485/lambda2,gamma2,sigma2,1)'];
            else
                function_string=['pressure_broadened_profile(2*pi*299792485./x,c,3,2*pi*299792485/lambda,gamma,sigma,1)+A2*pressure_broadened_profile(2*pi*299792485./x,c2,3,2*pi*299792485/lambda2,gamma2,sigma2,1)'];
            end
            fitted_peaks=zeros(length(index_peaks),9);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 500 500 500 500 1e4 1e4 1e4 1e4],'lower',[0 0 0 0 0 0 0 0 0] ); 
        case 9
            function_string=['pressure_broadening_from_stable_dist(2*pi*299792485./x,a,b,c,2*pi*299792485/lambda,omega_shift,gamma,sigma)'];
            fitted_peaks=zeros(length(index_peaks),7);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[2 1 1000 1000 1000 1000 1000],'lower',[ 0 -1 0 200 0 -1000 0] ); 
        case 10
            function_string=['pressure_broadening_from_stable_dist(2*pi*299792485./x,a,b,c,2*pi*299792485/lambda,omega_shift,gamma,sigma)+A2*pressure_broadening_from_stable_dist(2*pi*299792485./x,a2,b2,c2,2*pi*299792485/lambda2,omega_shift2,gamma2,sigma2)'];
            fitted_peaks=zeros(length(index_peaks),15);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[2 2 2 1 1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000],'lower',[ 0 0 0 -1 -1 0 0 200 200 0 0 -1000 -1000 0 0] );    
        case 11
            function_string=['pressure_broadening_from_stable_dist2(2*pi*299792485./x,c,c2,a_p,2*pi*299792485/lambda,gamma,sigma)'];
            fitted_peaks=zeros(length(index_peaks),6);
            if(b_guesses(1)>0)
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1 1000 1000 1000 1200 1000],'lower',[ 0 0 0 0 200 0] );    
            else
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 1000 1000 1000 1200 1000],'lower',[ 1 0 0 0 200 0] );    
            end   
        case 12
            function_string=['N1*pressure_broadening_from_stable_dist2(2*pi*299792485./x,c,c2,a_p,2*pi*299792485/lambda,gamma,sigma)+N2*pressure_broadening_from_stable_dist(2*pi*299792485./x,c21,c22,a_p2,2*pi*299792485/lambda2,gamma2,sigma2)'];
            fitted_peaks=zeros(length(index_peaks),14);
            if(b_guesses(1)>0&&b_guesses(2)>0)
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1.5 1.5 1 1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000],'lower',[ 0.5 0.5 0 0 0 0 0 0 0 0 200 200 0 0] );    
            elseif (b_guesses(1)<0&&b_guesses(2)>0)
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1.5 1.5 100 1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000],'lower',[ 0.5 0.5 1 0 0 0 0 0 0 0 200 200 0 0] );    
            elseif (b_guesses(1)<0&&b_guesses(2)<0)
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1.5 1.5 100 100 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000],'lower',[ 0.5 0.5 1 1 0 0 0 0 0 0 200 200 0 0] );    
            else
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1.5 1.5 1 100 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000],'lower',[ 0.5 0.5 0 1 0 0 0 0 0 0 200 200 0 0] );                    
            end     
        case 13
            function_string=['pressure_broadening_VdW_and_coloumb(2*pi*299792485./x,c,c2,a_p,2*pi*299792485/lambda,gamma,sigma,coloumb_width)'];
            fitted_peaks=zeros(length(index_peaks),7);
            if(b_guesses(1)>0)
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1 1000 1000 1000 1000 1000,1000],'lower',[ 0 0 0 0 0 200 0] );    
            else
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 1000 1000 1000 1000 1000 1000],'lower',[ 1 0 0 0 0 200 0] );    
            end 
    end
    opts.MaxFunEval=1000;
   opts.TolX=1e-7;
    %Translate back to non-logarithmic intensity, if needed
    if(log_intensity)
        intensity=exp(intensity);
        log_intensity=0;
    end
    
    opts.Display = 'iter';
    %initialize loop counter
    i=1;
    %Make a fit for every peak
    while i<=length(index_peaks)
        [xData, yData] = prepareCurveData( lambda(index_peaks(i)-lambda_intervall_index_length_m(i):index_peaks(i)+lambda_intervall_index_length_p(i)), intensity(index_peaks(i)-lambda_intervall_index_length_m(i):index_peaks(i)+lambda_intervall_index_length_p(i)) );
        yData=yData/max(yData);
         %opts.Weights=(yData/max(yData))+0.5;
        % Set up fittype and options.
        ft = fittype(function_string, 'independent', 'x', 'dependent', 'y' );
        
        %Set starting guesses, based on previous data treatment. The
        %parameters for the linearized respons function and extra baseline
        %compensation is assumed to be small and therefore we set zero as
        %starting guess
switch model
    case 1
        opts.StartPoint = [intensity(index_peaks(i)) gamma_guesses(i) lambda_peaks(i) sigma_guesses(i)];
    case 2
        opts.StartPoint = [intensity(index_peaks(i)) intensity(index_peaks(i+1)) gamma_guesses(i) gamma_guesses(i+1) lambda_peaks(i) lambda_peaks(i+1) sigma_guesses(i) sigma_guesses(i+1)];
    case 3
        opts.StartPoint = [intensity(index_peaks(i)) intensity(index_peaks(i+1)) intensity(index_peaks(i+2)) gamma_guesses(i) gamma_guesses(i+1) gamma_guesses(i+2) lambda_peaks(i) lambda_peaks(i+1) lambda_peaks(i+2) sigma_guesses(i) sigma_guesses(i+1) sigma_guesses(i+2)];            
    case 4
        opts.StartPoint = [intensity(index_peaks(i)) gamma_guesses(i) lambda_peaks(i) sigma_guesses(i) 0 0];
    case 5
        opts.StartPoint = [intensity(index_peaks(i)) intensity(index_peaks(i+1)) gamma_guesses(i) gamma_guesses(i+1) lambda_peaks(i) lambda_peaks(i+1) sigma_guesses(i) sigma_guesses(i+1) 0 0];
    case 6
        opts.StartPoint = [intensity(index_peaks(i)) gamma_guesses(i) gamma_guesses(i) lambda_peaks(i) sigma_guesses(i) sigma_guesses(i)];   
    case 7
        if(b_guesses(i)<0)
            function_string=['pressure_broadened_profile(2*pi*299792485./x,c,3,2*pi*299792485/lambda,gamma,sigma,-1)'];
        else
            function_string=['pressure_broadened_profile(2*pi*299792485./x,c,3,2*pi*299792485/lambda,gamma,sigma,1)'];
        end
        opts.StartPoint = [c_guesses(i) gamma_guesses(i) lambda_peaks(i) sigma_guesses(i)];
    case 8
        opts.StartPoint = [intensity(index_peaks(i+1)) c_guesses(i) c_guesses(i+1) gamma_guesses(i) gamma_guesses(i+1) lambda_peaks(i) lambda_peaks(i+1) sigma_guesses(i) sigma_guesses(i+1)];   
    case 9
        opts.StartPoint = [0.5 b_guesses(i) c_guesses(i) gamma_guesses(i) lambda_peaks(i) 0 sigma_guesses(i)];    
    case 10
        opts.StartPoint = [0.5 0.5 0.8 b_guesses(i) b_guesses(i+1) c_guesses(i) c_guesses(i+1) gamma_guesses(i) gamma_guesses(i+1) lambda_peaks(i) lambda_peaks(i+1) 0 0 sigma_guesses(i) sigma_guesses(i+1)];   
    case 11
            if(b_guesses(i)>0)
                opts.upper=[1 1000 1000 1000 1200 10000];
                opts.lower=[ 0 0 0 0 200 0];
               % opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1 1000 1000 1000 1000 1000],'lower',[ 0 0 0 200 0 0] );    
                opts.StartPoint = [1/PLRratio c_guesses(i) c_guesses(i) gamma_guesses(i) lambda_peaks(i) sigma_guesses(i)]; 
            else
               opts.upper=[100 1000 1000 1000 1200 10000];
               opts.lower=[ 1 0 0 0 200 0];                
               opts.StartPoint = [PLRratio c_guesses(i) c_guesses(i) gamma_guesses(i) lambda_peaks(i) sigma_guesses(i)]; 
               %opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 1000 1000 1000 1000 1000],'lower',[ 10 0 0 200 0 0] );    
            end
    case 12
        if(b_guesses(i)>0&&b_guesses(i+1)>0)
            opts.upper=[1.5 1.5 1 1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
            opts.lower=[0.5 0.5 0 0 0 0 0 0 0 0 200 200 0 0];
            opts.StartPoint = [1 1 1/PLRratio 1/PLRratio c_guesses(i) c_guesses(i) c_guesses(i+1) c_guesses(i+1) gamma_guesses(i) gamma_guesses(i+1) lambda_peaks(i) lambda_peaks(i+1) sigma_guesses(i) sigma_guesses(i+1)];             
        elseif(b_guesses(i)<0&&b_guesses(i+1)>0)
            opts.upper=[1.5 1.5 100 1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
            opts.lower=[0.5 0.5 1 0 0 0 0 0 0 0 200 200 0 0];
            opts.StartPoint = [1 1 PLRratio 1/PLRratio c_guesses(i) c_guesses(i) c_guesses(i+1) c_guesses(i+1) gamma_guesses(i) gamma_guesses(i+1) lambda_peaks(i) lambda_peaks(i+1) sigma_guesses(i) sigma_guesses(i+1)];             
        elseif(b_guesses(i)<0&&b_guesses(i+1)<0)
            opts.upper=[1.5 1.5 1 1 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
            opts.lower=[ 0.5 0.5 0 0 0 0 0 0 0 0 200 200 0 0];
            opts.StartPoint = [1 1 PLRratio PLRratio c_guesses(i) c_guesses(i) c_guesses(i+1) c_guesses(i+1) gamma_guesses(i) gamma_guesses(i+1) lambda_peaks(i) lambda_peaks(i+1) sigma_guesses(i) sigma_guesses(i+1)];                         
        else
            opts.upper=[1.5 1.5 100 100 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
            opts.lower=[0.5 0.5 1 1 0 0 0 0 0 0 200 200 0 0];
            opts.StartPoint = [1 1 1/PLRratio PLRratio c_guesses(i) c_guesses(i) c_guesses(i+1) c_guesses(i+1) gamma_guesses(i) gamma_guesses(i+1) lambda_peaks(i) lambda_peaks(i+1) sigma_guesses(i) sigma_guesses(i+1)];             
        end
    case 13
            if(b_guesses(i)>0)
                opts.upper=[1 1000 1000 1000 1000 1000 10000];
                opts.lower=[ 0 0 0 0 200 0 0];
               % opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[1 1000 1000 1000 1000 1000],'lower',[ 0 0 0 200 0 0] );    
                opts.StartPoint = [1/PLRratio c_guesses(i) c_guesses(i) coloumb_width_guesses(i) gamma_guesses(i) lambda_peaks(i) sigma_guesses(i)]; 
            else
               opts.upper=[100 1000 1000 1000 1000 1000 10000];
               opts.lower=[ 1 0 0 0 200 0 0];                
               opts.StartPoint = [PLRratio c_guesses(i) c_guesses(i) coloumb_width_guesses(i) gamma_guesses(i) lambda_peaks(i) sigma_guesses(i)]; 
               %opts = fitoptions( 'Method', 'NonlinearLeastSquares','upper',[100 1000 1000 1000 1000 1000],'lower',[ 10 0 0 200 0 0] );    
            end         
end
        
        
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        
        % Plot fit with data.
        figure( 'Name', 'untitled fit 1' );
        h = plot( fitresult, xData, yData );
        legend( h, 'intensity vs. lambda', 'untitled fit 1', 'Location', 'NorthEast' );
        % Label axes
        xlabel lambda
        ylabel intensity
        grid on
       switch model
           case 1
                fitted_peaks(i,:)=[fitresult.A fitresult.lambda fitresult.gamma fitresult.sigma];
                i=i+1;
           case 2     
               fitted_peaks(i,:)=[fitresult.A fitresult.lambda fitresult.gamma fitresult.sigma fitresult.A2 fitresult.lambda2 fitresult.gamma2 fitresult.sigma2];
               i=i+2;
           case 3
               fitted_peaks(i,:)=[fitresult.A fitresult.lambda fitresult.gamma fitresult.sigma fitresult.A2 fitresult.lambda2 fitresult.gamma2 fitresult.sigma2 fitresult.A3 fitresult.lambda3 fitresult.gamma3 fitresult.sigma3 ];
               i=i+3;
           case 4
               fitted_peaks(i,:)=[fitresult.A fitresult.lambda fitresult.gamma fitresult.sigma fitresult.t fitresult.z];
                i=i+1;
           case 5 
               fitted_peaks(i,:)=[fitresult.A fitresult.lambda fitresult.gamma fitresult.sigma fitresult.A2 fitresult.lambda2 fitresult.gamma2 fitresult.sigma2 fitresult.t fitresult.z];
               i=i+2;
           case 6
               fitted_peaks(i,:)=[fitresult.A fitresult.lambda fitresult.gamma1 fitresult.gamma2 fitresult.sigma2 fitresult.sigma2];
               i=i+1;
           case 7
               fitted_peaks(i,:)=[fitresult.c fitresult.lambda fitresult.gamma fitresult.sigma];
               i=i+1;
           case 8
               fitted_peaks(i,:)=[fitresult.A2 fitresult.c fitresult.c2 fitresult.lambda fitresult.lambda2 fitresult.gamma fitresult.gamma2 fitresult.sigma fitresult.sigma2];
               i=i+2;
           case 9
               fitted_peaks(i,:)=[fitresult.omega_shift fitresult.lambda fitresult.omega_shift fitresult.a fitresult.b fitresult.c fitresult.gamma fitresult.sigma];
               i=i+1;
           case 10
               fitted_peaks(i,:)=[fitresult.lambda fitresult.lambda2 fitresult.omega_shift fitresult.omega_shift2 fitresult.A2 fitresult.a fitresult.a2 fitresult.b fitresult.b2 fitresult.c fitresult.c2 fitresult.gamma fitresult.gamma2 fitresult.sigma fitresult.sigma2];
               i=i+2;
           case 11
               fitted_peaks(i,:)=[fitresult.lambda fitresult.c fitresult.c2 fitresult.a_p fitresult.gamma fitresult.sigma];
               i=i+1; 
           case 12
               fitted_peaks(i,:)=[fitresult.lambda fitresult.lambda2 fitresult.c fitresult.c21 fitresult.c2 fitresult.c22 fitresult.a_p fitresult.a_p2 fitresult.gamma fitresult.gamma2 fitresult.N1 fitresult.N2 fitresult.sigma fitresult.sigma2];
               i=1+2;
           case 13
               fitted_peaks(i,:)=[fitresult.lambda fitresult.c fitresult.c2 fitresult.coloumb_width fitresult.a_p fitresult.gamma fitresult.sigma];
               i=i+1;
       end
    end