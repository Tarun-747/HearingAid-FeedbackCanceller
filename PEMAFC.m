
clc;
close all;
clear all;
disp('Simulation Started');
aa_1 = [zeros(1,60) 0 0 -0.00125 -0.065 0.15 0.12 -0.035 -0.13 -0.07 0.02 0.07  0.02 -0.05 -0.005 0.025 0.03 -0.005 -0.035 0.04 0.035 -0.045 -0.025 0.025 0.02 -0.02 -0.025 0.0175 0.022 0.0175 -0.02 -0.015 0.015 0.005 -0.01 -0.005 zeros(1,65)]';

Input=filter(1,[1 -0.8], randn(1,200000));
V1=load('V.mat');% Impulse response
IN=50000*(V1.V);%Impulsive Noise
N=length(Input);%Input sample length
IN=IN(1:N); %Impulsive Noise


x=Input;

N=length(Input(1:200000));%Input sample length
M=length(aa_1);%Length of Impulse response of acoustic paths

F_dB=40;%forward path gain in dB scale
F_amp = 8;%10.^(F_dB./20);%forward path gain in normal scale
delay_forw=48;
delta = 1e-5;
mu_h=0.01;

Nfreq=M;
H_1=fft(aa_1,Nfreq);

for itr=1:1
e_delay = zeros(N+delay_forw,1); 
feedback_tap= zeros(M,1);
feedback_cancel_tap= zeros(M,1);
h_cap= zeros(M,1);

Ma=21;% Order of AR model
Frame_size = 160;
tap=zeros(Frame_size,1); % frame size for Levinson-Durbin recurtion-1


error_tap=zeros(Ma,1);
loudspeaker_tap=zeros(Ma,1);
cancel_tap=zeros(M,1);
cancel_update_tap=zeros(M,1);

h_hat= zeros(M,1);
H_hilbert=[4.54486433487e-05, -0.03191360456117,-3.31277815639e-05, -0.02604783588745,-1.594271701847e-06, -0.03732840213574,3.469425875731e-06, -0.05311074599375,-1.033785626931e-06, -0.07669953750135,-5.966834480336e-06,  -0.1168538267669,1.666246290771e-05,  -0.2058058100572,-1.893028029711e-05,  -0.6344726293921, 0,   0.6344726293921,1.893028029711e-05,   0.2058058100572,-1.666246290771e-05,   0.1168538267669,5.966834480336e-06,  0.07669953750135, 1.033785626931e-06,  0.05311074599375,-3.469425875731e-06,  0.03732840213574, 1.594271701847e-06,  0.02604783588745, 3.31277815639e-05,  0.03191360456117, -4.54486433487e-05]';
Hilbert_length=length(H_hilbert);
Hilbert_Delay=16;
u_tap=zeros(Hilbert_length,1);

for n=1:N
    
    h=aa_1;
    u(n)=F_amp*e_delay(n);%loudspeaker output

    feedback_tap=[u(n); feedback_tap(1:end-1)];
    feedback_cancel_tap=[u(n); feedback_cancel_tap(1:end-1)];
    f(n)=feedback_tap'*h;%output of feedback path-1
    
    m(n) = x(n) + f(n)+ IN(n); %output of microphone-1
    f_cap(n)=feedback_cancel_tap'*h_cap;%output of feedback canceller-1
    e(n)=m(n)-f_cap(n); % subtract feedback canceller-1 output from microphone-1 output
    e_delay(n+delay_forw)=e(n); %insert delay in the forward path
    
    cancel_tap=[u(n); cancel_tap(1:end-1)]; %Input vector to  adaptive shadow filter-1
    f_hat(n)=cancel_tap'*h_hat; % output of adaptive shadow filter
    error(n)=m(n)-f_hat(n);%error used to update shadow filter-1
    
    tap=[e(n);tap(1:end-1,1)];
    [r,lg] = xcorr(tap,'biased');
    r(lg<0) = [];
    L1 = levinson(r,Ma-1); %Levinson-Durbin algorithm-1
    L=L1';
    
    error_tap=[error(n); error_tap(1:end-1)];
    error_prefilter(n)=error_tap'*L;%error signal prefiltered
    
    loudspeaker_tap=[u(n); loudspeaker_tap(1:end-1)];
    loudspeaker_prefiltered(n)=loudspeaker_tap'*L;%prefiltering of loudspeaker signal
    cancel_update_tap=[loudspeaker_prefiltered(n); cancel_update_tap(1:end-1)]; %Input vector to  adaptive shadow filter-1
    %=====================USE Different Algorithm==============   
    h_hat = h_hat + (mu_h  / (norm(cancel_update_tap)^2 + delta)) * cancel_update_tap .* error_prefilter(n);%update shadow filter-1
    
    h_cap=h_hat;
    %==================================================
    
    HH_1_cap = fft( h_cap,Nfreq);
    
    diff = H_1(1:Nfreq/2+1) - HH_1_cap(1:Nfreq/2+1);% fft difference-1
    
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
plot(t,MIS,'r','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('Misalignment (dB)');

figure;
plot(t,MSG,'r','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('MSG (dB)');

figure;
plot(t,ASG,'r','linewidth',2);grid on;
xlabel('Time (sec)');
ylabel('ASG (dB)'); 


figure;
plot(h);
hold on;
plot(h_cap,'r');

disp('Simulation Completed');




