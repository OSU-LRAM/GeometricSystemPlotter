function [pinch_wave_matched]= pinched_sin_scale

% Generate a scale factor so that the maximum angle reached by the sin and
% pinched sin basis functions is the same
n = logspace(-1,1.5,100);

x = linspace(0,pi/2);

base_area = trapz(x,cos(x));

for idx = 1:numel(n)
    
    scale_factor(idx) = base_area/trapz(x,cos(x).^n(idx));
    
end

P = polyfit(n,scale_factor,3);

% figure(18)
% 
% plot(n,scale_factor,n,polyval(P,n))


% Now generate the pinched sin function

syms s n A phi a1 a2 

% Raw wave function
raw_wave = cos(2*pi*s+phi);

% Pinched in wave
pinched_wave = A * (abs(raw_wave))^n * sign(raw_wave);

% Substitute in the a1 and a2 values
amp = (a1^2+a2^2)^.5;
phase = atan2(a2,a1);
pinched_wave_restored = subs(pinched_wave,[A phi],[amp phase]);


% % Wave with sin and cosine amplitudes
% raw_wave = (a1*cos(2*pi*s)+a2*sin(2*pi*s));
% 
% % Amplitude of the wave
% amp_wave = (a1^2+a2^2)^(.5);
% 
% % Normalize the wave
% norm_wave = raw_wave/amp_wave;
% 
% % Pinch the wave
% pinched_wave = (abs(norm_wave)^n) * sign(norm_wave);
% 
% % Restore the scale of the wave
% pinched_wave_restored = pinched_wave*amp_wave;




% Scale to match max angle with original wave
pinch_wave_matched = pinched_wave_restored * (P*[n^3;n^2;n;1]);





end
