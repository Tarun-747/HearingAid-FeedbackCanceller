% Feedback cancellation in Single Micropone Single loudspeaker based Digital Hearing Aid
% September 2017
% Somanath Pradhan
% Use of this code with out prior permission is a crime
%% *************************************************************************
clc;
%close all;
clear all;
disp('Simulation Started');
aa_1 = [zeros(1,60) 0 0 -0.00125 -0.065 0.15 0.12 -0.035 -0.13 -0.07 0.02 0.07  0.02 -0.05 -0.005 0.025 0.03 -0.005 -0.035 0.04 0.035 -0.045 -0.025 0.025 0.02 -0.02 -0.025 0.0175 0.022 0.0175 -0.02 -0.015 0.015 0.005 -0.01 -0.005 zeros(1,65)]';

% a_1=load('Normal_Impulse.mat');% Impulse response
% aa_1=(a_1.Normal_Impulse)';
 
V1=load('V.mat');% Impulse response
IN=50000*(V1.V);%Impulsive Noise

%[speech, fs] = audioread('C:\Users\mtaru\Desktop\Media1.wav');
%speech = speech(:, 1); % Use only the first channel if stereo
%if fs ~= 8000
%    speech = resample(speech, 8000, fs);
%    fs = 8000;
%end
%Input = speech';
%N = length(Input); % Input sample length


% s1=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp01.wav');
% s2=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp02.wav');
% s3=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp03.wav');
% s4=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp04.wav');
% s5=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp05.wav');
% s6=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp06.wav');
% s7=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp07.wav');
% s8=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp08.wav');
% s9=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp09.wav');
% s10=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp10.wav');
% s11=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp11.wav');
% s12=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp12.wav');
% s13=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp13.wav');
% s14=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp14.wav');
% s15=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp15.wav');
% s16=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp16.wav');
% s17=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp17.wav');
% s18=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp18.wav');
% s19=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp19.wav');
% s20=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp20.wav');
% s21=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp21.wav');
% s22=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp22.wav');
% s23=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp23.wav');
% s24=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp24.wav');
% s25=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp25.wav');
% s26=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp26.wav');
% s27=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp27.wav');
% s28=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp28.wav');
% s29=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp29.wav');
% s30=audioread('C:\Users\Dell\Desktop\HA_Code\Wideband_Speech\sp30.wav');

% Input Signal
Input=filter(1,[1 -0.8], randn(1,200000));
N=length(Input);%Input sample length
IN=IN(1:N); %Impulsive Noise

M=length(aa_1);%Length of Impulse response of acoustic paths

F_dB=40;%forward path gain in dB scale
F_amp =8;%10.^(F_dB./20);%forward path gain in normal scale
delay_forw=48;
delta = 1e-5;
mu_h=0.01;

b1=0.5;
b2=0.5;%beamformer coefficient
Nfreq=M;
H_1=fft(aa_1,Nfreq);


for itr=1:1
e_bar_delay = zeros(N+delay_forw,1); 
feedback_tap_1= zeros(M,1);
feedback_tap_2= zeros(M,1);
feedback_cancel_tap_1= zeros(M,1);
feedback_cancel_tap_2= zeros(M,1);
h_1_cap= zeros(M,1);
h_2_cap= zeros(M,1);

for n=1:N
    
    h_1=aa_1;

    x_1(n)=Input(n); %input to microphone-1=speech+impulsive noise
    u(n)=F_amp*e_bar_delay(n);%loudspeaker output

    feedback_tap_1=[u(n); feedback_tap_1(1:end-1)];%tap delayed vector
    
    feedback_cancel_tap_1=[u(n); feedback_cancel_tap_1(1:end-1)];%tap delayed vector
    
    f_1(n)=feedback_tap_1'*h_1;%output of feedback path-1
    
    m_1(n) = x_1(n) + f_1(n)% + IN ; %output of microphone-1 with impulsive noise
    
    f_1_cap(n)=feedback_cancel_tap_1'*h_1_cap;%output of feedback canceller-1
    
    e(n)=m_1(n)-f_1_cap(n); % subtract feedback canceller-1 output from microphone-1 output
        
    e_bar_delay(n+delay_forw)=e(n); %insert delay in the forward path
    
    h_1_cap = h_1_cap + (mu_h  / (norm(feedback_cancel_tap_1)^2 + delta)) * feedback_cancel_tap_1 .* e(n);%update canceller-1
    
    H_1_cap = fft( h_1_cap,Nfreq);
    
    diff = H_1(1:Nfreq/2+1) - H_1_cap(1:Nfreq/2+1);% fft difference-1
    
    MSG_without(n)=20*log10 ( min (1./  (F_amp*abs(diff) )  )  );% maximum stable gain without AFC
    
    MSG(n)=20*log10 ( min (1./  abs(diff)   )  );% maximum stable gain

    ASG(n)=MSG(n)-( 20*log10( min( 1./ (  abs(H_1(1:Nfreq/2+1)))  )  )  );% Added stable gain
     
    MIS(n)=10*log10( ( (norm(diff)^2)  ) / ( (norm(H_1(1:Nfreq/2))^2)  ));% Misalignment

    
   
end
clc
iteration=itr

end
%% Performance Measures
fs=8000;  
t = linspace(0,(N-1)/fs,N);


figure;
plot(t,MIS,'b','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('Misalignment (dB)');

figure;
plot(t,MSG,'b','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('MSG (dB)');

figure;
plot(t,ASG,'b','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('ASG (dB)');

figure;
plot(h_1);hold on;
plot(h_1_cap,'r');


% figure;
% plot(t,MSG_without,'b','linewidth',2);grid on;
% xlabel('Time (sec)');
% ylabel('MSG (dB)');
% 
% figure;
% plot(autocorr(Input,'NumLags',100))

display('Simulation Completed')



% [hz1,w1]=freqz(aa_1,1,512,16000);
% phi1=180*unwrap(angle(hz1))/pi;
% [hz2,w2]=freqz(aa_2,1,512,16000);
% phi2=180*unwrap(angle(hz2))/pi;
% 
% [hz3,w3]=freqz(aa_3,1,512,16000);
% phi3=180*unwrap(angle(hz3))/pi;
% [hz4,w4]=freqz(aa_4,1,512,16000);
% phi4=180*unwrap(angle(hz4))/pi;
% 
% figure;
% plot(w1,20*log10(abs(hz1)),'b','linewidth',2);hold on;
% plot(w2,20*log10(abs(hz2)),'r','linewidth',2);grid on;
% plot(w3,20*log10(abs(hz3)),'g','linewidth',2);hold on;
% plot(w4,20*log10(abs(hz4)),'black','linewidth',2);grid on;
% xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
% legend('|G_1(\omega)|','|G_2(\omega)|','|G_1(\omega)|(Close object)','|G_2(\omega)|(Close object)');

