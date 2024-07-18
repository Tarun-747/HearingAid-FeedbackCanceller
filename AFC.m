%% Feedback Cancellation in Hearing Aids using Adaptive Decorrelation NLMS
% 04 April 2017
% Based on prediction error method, which used Levinson Durbin Algorithm
clc;
clear all;
%close all;
%% Initialization
speech_1=load('INPUT.mat');
speech=speech_1.INPUT;
fp_gain=5; % Forward path gain
feedback_path=[zeros(1,60) 0 0 -0.00125 -0.065 ...
    0.15 0.12 -0.035 -0.13 -0.07 0.02 0.07  0.02 ...
    -0.05 -0.005 0.025 0.03 -0.005 -0.035 0.04 ...
    0.035 -0.045 -0.025 0.025 0.02 -0.02 -0.025 ...
    0.0175 0.022 0.0175 -0.02 -0.015 0.015 0.005 ...
    -0.01 -0.005 zeros(1,25)];
length_input=length(speech);
N=length_input;
delay_length=60;
canceller_length=length(feedback_path);
error=zeros(1, length_input);
error_delayed=zeros(1, length_input+delay_length);
feedback_tap=zeros(1,length(feedback_path));
canceller_tap=zeros(1,canceller_length);
canceller_copy_tap=zeros(1,canceller_length);
feedback_canceller=zeros(1,canceller_length);
feedback_canceller_copy=zeros(1,canceller_length);
weight_update_tap=zeros(1,canceller_length);
q_filter=zeros(1,51);
q_tap=zeros(1,length(q_filter));
l_filter=zeros(1, length(q_filter)+1); % Adaptive decorrelation filter
l1_tap=zeros(1,length(l_filter));
l2_tap=zeros(1,length(l_filter));
frame_size=160;
pem_tap=zeros(1,frame_size);
deltaW=zeros(1,length_input);
h1=zeros(1,length_input);

V1=load('V.mat');% Impulse response
%variance =15;
SNR=5;
n_std=10^(-SNR/20)*sqrt(var(speech));
noise=n_std*randn(length(speech),1);
noise=1*noise;

outlierr=zeros(1,length_input);
outlierr(1e5)=0;

% Adjust the outlier generation
% n_std = 10^(-SNR/20) * sqrt(var(speech));
% outlierr = zeros(length(speech), 1); % Initialize the outlier vector
% 
% % Probability and magnitude for outliers
% outlier_prob = 0.01; % Probability of an outlier (adjust to control sparsity)
% outlier_mean = 0; % Mean of the Gaussian distribution
% outlier_std_dev = n_std; % Standard deviation of the Gaussian distribution
% 
% % Generate sparse Gaussian outliers
% outlier_indices = rand(length(speech), 1) < outlier_prob; % Logical indices for outliers
% outlierr(outlier_indices) = outlier_mean + outlier_std_dev * randn(sum(outlier_indices), 1); % Assign Gaussian outliers
%outlierr=0*outlierr;
% Algorithm name
% Parameters
%%
%%
%NLMS
mu_NLMS=0.002;

%ZA-NLMS
rho_ZA=6e-9;
mu_ZANLMS=2e-3;

%RZA-NLMS
rho_RZA=1e-8;
eps_RZA=800;
mu_RZANLMS=2e-3;
%l0-NLMS
mu_l0NLMS=2e-3;
eps_l0=1e-2;
rho_l0=1e-20;
delta=1e-8;
fft_length=120;
 %% Feedback Cancellation
 for n=1:length_input
      n
     forward_opt(n)=fp_gain*error_delayed(n);
     feedback_tap=[forward_opt(n) feedback_tap(1:end-1)];
     feedback_signal(n)=feedback_tap*feedback_path';
     mic_signal(n) = speech(n) + feedback_signal(n) + noise(n) + outlierr(n);
     canceller_copy_tap=[forward_opt(n) canceller_copy_tap(1:end-1)];
     canceller_copy_opt(n)=canceller_copy_tap*feedback_canceller_copy';
     error(n)= mic_signal(n)-canceller_copy_opt(n);
     error_delayed(n+delay_length)=error(n);
     
     canceller_tap=[forward_opt(n) canceller_tap(1:end-1)];
     canceller_opt(n)=canceller_tap*feedback_canceller';
     error_pem(n)=mic_signal(n)-canceller_opt(n);
     pem_tap=[error_pem(n) pem_tap(1:end-1)];
     [r,lg]=xcorr(pem_tap','biased'); 
     r(lg<0)=[];
     pem=levinson(r,length(q_filter));
     
     l1_tap=[error_pem(n) l1_tap(1:end-1)];
     l1_opt(n)=l1_tap*pem';
     l2_tap=[forward_opt(n) l2_tap(1:end-1)];
     l2_opt(n)=l2_tap*pem';
     weight_update_tap=[l2_opt(n) weight_update_tap(1:end-1)];
     
     %% Different Algorithms
     ez=feedback_canceller;
     %NLMS
     feedback_canceller=ez + (mu_NLMS/(weight_update_tap*weight_update_tap'+delta))*weight_update_tap.*l1_opt(n); % NLMS Weight Update
     deltaW(n)=norm(feedback_canceller-ez);
     h1(n)=norm((mu_NLMS/(weight_update_tap*weight_update_tap'+delta))*weight_update_tap);
     %ZA-NLMS
     %feedback_canceller=ez -rho_ZA*sign(ez) + (mu_ZANLMS/(weight_update_tap*weight_update_tap'+delta))*weight_update_tap.*l1_opt(n); % ZA_NLMS Weight Update
     
     %RZA-NLMS
     %feedback_canceller=ez -rho_RZA*sign(ez)./(1+eps_RZA*abs(ez)) + (mu_RZANLMS/(weight_update_tap*weight_update_tap'+delta))*weight_update_tap.*l1_opt(n); % RZA_NLMS Weight Update
     
     %New l0 norm -NLMS
     %feedback_canceller=ez -rho_l0*((ez)./(ez.^2+eps_l0^2)-(ez.^3)./(ez.^2+eps_l0^2).^2) + (mu_l0NLMS/(weight_update_tap*weight_update_tap'+delta))*weight_update_tap.*l1_opt(n); % RZA_NLMS Weight Update
     %%%%%%%%%%%%%%%%%%%%%%%%
     feedback_canceller_copy=feedback_canceller;
%      figure(1)
%      plot(feedback_canceller_copy);
     %% Performance Measure
     canceller_fft=fft(feedback_canceller, fft_length);
     feedback_path_fft=fft(feedback_path, fft_length);
     fft_difference=feedback_path_fft(1:fft_length/2+1)...
         -canceller_fft(1:fft_length/2+1);
     MSG(n)=20*log10(min(1./abs(fft_difference))); % The maximum difference between the two paths in the frequency domain for every iteration
     MIS(n)=10*log10((norm(fft_difference)^2)/(norm(feedback_path_fft(1:fft_length/2+1)))^2);
     ASG(n)=MSG(n)-20*log10(min(1./abs(feedback_path_fft(1:fft_length/2+1))));
 end
 % figure;
 % plot(feedback_path);
 % hold on;
 % plot(feedback_canceller,'r');
 % figure
 % subplot(3,1,1);
 % plot(MSG);
 % xlabel('Samples');
 % ylabel('MSG (dB)');
 % subplot(3,1,2);
 % plot(MIS);
 % xlabel('Samples');
 % ylabel('MIS (dB)');
 % subplot(3,1,3);
 % plot(ASG);
 % xlabel('Samples');
 % ylabel('ASG (dB)');
%figure
hold on

plot(MIS);
xlabel('Samples');
ylabel('MIS (dB)');
 %% Testing
%  sound(speech(1:100000),16000);   % Original speech
%  sound(mic_signal(1:100000),16000); % Speech + Feedback
%  sound(forward_opt(1:100000),16000); % As heard by the user
 audiowrite('Processed_Basic.wav',forward_opt(end-50000:end),8000);
 audiowrite('Clean_Speech.wav',speech(end-50000:end)+noise(end-50000:end),8000);
 pesq_score=pesq('Clean_Speech.wav','Processed_Basic.wav')