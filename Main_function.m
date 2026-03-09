%% Clean the memory and worksapace 
close all ;
clear     ;
clc       ;

%% Configure the system simulation condition 
fs  =   16000    ; % The system sampling rate.
T   =   3        ; % The duration of the simulation.
t   =   0:1/fs:(T-1/fs) ; 
N   =   length(t); % The number of the data.

Len_N = 512      ; % Seting the length of the control filter.
snrr  = 50;
noise_type = 'white';  %white noise and pulse noise

%% Build the broad band noise for training set
% Loading path 
load('path\P1.mat')    ;
load('path\S11.mat')   ;
Pri_path = conv(P1,S11);


Track_num = 3         ; % Seting the number of the track for the trainning noise. 
if exist('Primary_noise.mat', 'file') == 2
    disp('Primary_noise exists in the current path.\n');
    % Loading the primary noise 
    load('Primary_noise.mat');
    load('Disturbance.mat')  ;
    load('Reference.mat')    ;
else
Noise     = randn(N,1);
% filter 
filter_1 = fir1(512,[400 1000]*2/fs);
filter_2 = fir1(512,[800 2200]*2/fs) ;
filter_3 = fir1(512,[2000,3000]*2/fs) ;

% Primary noise 
Pri_1 = filter(filter_1,1,Noise) ;
Pri_2 = filter(filter_2,1,Noise) ;
Pri_3 = filter(filter_3,1,Noise) ;
aircrat = load('.\Real_recorded_noise\707_Sound_for_Simulation.mat'); 
Pri_4 = aircrat.PilingNoise(fs*4+1:fs*7); % aircraft noise
helicpt = audioread('.\Real_recorded_noise\Helicopter.wav');
Pri_5 = helicpt(1:3:fs*T*3);  %downsample 48k to 16k 

t = 0:1/fs:(T-1/fs);
f0 = [500,1200,3000];
PolyNoise = (cos(2*pi*f0(1)*t) + cos(2*pi*f0(2)*t) + cos(2*pi*f0(3)*t))/3;
Pri_6 = PolyNoise'; %multi-tone noise

% Drawing fiture 
data = [Pri_1,Pri_2,Pri_3,Pri_4,Pri_5,Pri_6];
figure ;
len_fft = length(Pri_1)   ;
len_hal = round(len_fft/2);
title('Frequency spectrum of primary noises')
for ii = 1:Track_num
    freq = 20*log(abs(fft(data(:,ii))));
    subplot(Track_num,1,ii);
    plot(0:(fs/len_fft):(len_hal-1)*(fs/len_fft), freq(1:len_hal));
    grid on   ;
    title("The "+num2str(ii)+"th primary noise")
    xlabel('Frequency (Hz)')
end
% Save primary noise into workspace 
save('Primary_noise.mat','Pri_1','Pri_2','Pri_3','Pri_4','Pri_5');
% Generating Distrubance 
Dis_1 = filter(Pri_path,1,Pri_1);
Dis_2 = filter(Pri_path,1,Pri_2);
Dis_3 = filter(Pri_path,1,Pri_3);
Dis_4 = filter(Pri_path,1,Pri_4);
Dis_5 = filter(Pri_path,1,Pri_5);
Dis_6 = filter(Pri_path,1,Pri_6);
% Save distrubancec into workspace 
save('Disturbance.mat','Dis_1','Dis_2','Dis_3','Dis_4','Dis_5','Dis_6');
% Genrating Filtered reference signal 
Rf_1 = filter(S11,1,Pri_1);
Rf_2 = filter(S11,1,Pri_2);
Rf_3 = filter(S11,1,Pri_3);
Rf_4 = filter(S11,1,Pri_4);
Rf_5 = filter(S11,1,Pri_5);
Rf_6 = filter(S11,1,Pri_6);
% Save filter reference signal into workspace 
save('Reference.mat','Rf_1','Rf_2','Rf_3','Rf_4','Rf_5','Rf_6');
end

%% Radomly sampling the noise tracks to build dataset for the MAML algorithm

if exist('Sampe_data_N_set.mat', 'file') == 2
    disp('Sampe_data_N_set in the current path.\n');
    load('Sampe_data_N_set.mat');
else
N_epcho  = 3e5                         ; % Setting the number of the epcho 
Trac     = randi(Track_num,[N_epcho,1]); % Randomly choosing the different tracks. 
Len_data = length(Dis_1)               ;
% Seting the N steps 
len   = 2*Len_N -1 ;
Fx_data = zeros(Len_N,N_epcho);
Di_data = zeros(Len_N,N_epcho);

Ref_data = [Rf_1,Rf_2,Rf_3,Rf_4,Rf_5,Rf_6]   ;
Dis_data = [Dis_1,Dis_2,Dis_3,Dis_4,Dis_5,Dis_6];
% add noise 
Dis_data = add_noise(Dis_data, snrr, noise_type);
Ref_data = add_noise(Ref_data, snrr, noise_type);

% Randomly Sampling
for jj = 1:N_epcho
    End = randi([len,Len_data]);
    Di_data(:,jj) = Dis_data(End-511:End,Trac(jj));
    Fx_data(:,jj) = Ref_data(End-511:End,Trac(jj));
end
save('Sampe_data_N_set.mat','Di_data','Fx_data');
end

%% Using Modified MAML algorithm to get the best initial control filter

if exist('Step_initiate_Nstep_forget.mat', 'file') == 2
    disp('Step_initiate_Nstep_forget in the current path.\n');
    load('Step_initiate_Nstep_forget.mat');
else
% Create a MAML algorithm
a  = MAML_train_test_forget(0,Len_N);
N  = size(Di_data,2)         ; % The number of the sample in training set.
Er = zeros(N,1)              ; % Residual error vector
miu = zeros(N,1)             ; % the initial step size

% Seting the initial value of the control filter, 
W    = zeros(Len_N,1)        ;
% Seting the forget factor 
lamda = 0.9                ;
% Seting the learning for MAML 
epslon = 1e-7 ;
% Runing the MAML algorithm 
for jj = 1:N
    [a, Er(jj)] = a.MAML_initial(Fx_data(:,jj),Di_data(:,jj),lamda,epslon);
    miu(jj) = a.Phi;
end

% Save the MAML step size
Mux = miu(end);
save('Step_initiate_Nstep_forget.mat','Mux','Mux');

% Drawing the residual error of the Modified MAML algorihtm 
cc = orderedcolors('gem12'); % Set the color of lines
figure;subplot(2,1,1)
plot(Er)  ;
grid on   ;
title('(a) Learning curve','Interpreter','latex');
xlabel('Epoch','Interpreter','latex');
ylabel('Residual error','Interpreter','latex');
xlim([-1e2,1e5])
ylim([-10,10])

subplot(2,1,2)
plot(miu,'Color',cc(2,:));
grid on   ;
title('(b) Step size iteration process','Interpreter','latex');
xlabel('Epoch','Interpreter','latex');
ylabel('Step size','Interpreter','latex');
xlim([-1e2,1e5])
ylim([6e-4,1.2e-3])

end

%% Testing broadband noise cancellation by using MAML 
% generate broadband noise
Noise     = randn(fs*10,1); 
filter_1 = fir1(512,[200 800]*2/fs);
Pri_1 = filter(filter_1,1,Noise) ;
% Generating the disturbacne 
Dis_1   = filter(Pri_path,1,Pri_1);
% Generating the filter reference             
Rf_1    = filter(S11,1,Pri_1);
% add noise
Rf_1 = add_noise(Rf_1, snrr, noise_type);
Dis_1 = add_noise(Dis_1, snrr, noise_type);


% Runing the FxLMS with the empirical step size value    
Wc_initial = zeros(Len_N,1);
muw1        = 1/(norm(Rf_1))^2/2; % empirical step size value
Er1        = FxLMS(Len_N, Wc_initial, Dis_1, Rf_1, muw1);

% Runing the FxNLMS step size      
Wc_initial = zeros(Len_N,1);
Beta       = 0.05 ; % the parameter for adjusting the step size
[Er2,muw2] = FxNLMS(Len_N, Wc_initial, Dis_1, Rf_1, Beta);

% Runing the FxLMS with the VSS    
Wc_initial = zeros(Len_N,1);
Beta       = 0.005 ; % the parameter for adjusting the step size
[Er3,muw3] = FxLMS_VSS2(Len_N, Wc_initial, Dis_1, Rf_1, Beta);

% Runing the FxLMS with the CSS   
Wc_initial = zeros(Len_N,1);
Beta       = 0.75 ; % the parameter for adjusting the step size
mu1        = 0.001;
mu2        = 0.0001;
[Er4,muw4] = FxLMS_CSS(Len_N, Wc_initial, Dis_1, Rf_1, Beta, mu1, mu2);

% Runing the FxLMS with the MAML step size    
Wc_initial = zeros(Len_N,1);
muw5        = Mux;
[Er5] = FxLMS(Len_N, Wc_initial, Dis_1, Rf_1, muw5);


%% Drawing results
cc = orderedcolors('gem12'); % Set the color of lines
figure; set(groot,'defaultAxesTickLabelInterpreter','latex');
subplot(2,1,1);
plot((0:length(Dis_1)-1)*(1/fs),Dis_1,'Color',cc(1,:));
hold on
plot((0:length(Er1)-1)*(1/fs),Er1,'Color',cc(2,:));
plot((0:length(Er2)-1)*(1/fs),Er2,'Color',cc(3,:));
plot((0:length(Er3)-1)*(1/fs),Er3,'Color',cc(9,:));
plot((0:length(Er4)-1)*(1/fs),Er4,'Color',cc(5,:));
plot((0:length(Er4)-1)*(1/fs),Er5,'Color',cc(7,:));
xlim([0.1,5]);ylim([min(Dis_1),max(Dis_1)])
title('(a) Error signal','Interpreter','latex')
xlabel('Time (seconds)','Interpreter','latex')
ylabel('Magnitude','Interpreter','latex')
legend({'ANC off','Empirical step size','Normalized step size', ...
    'Variable step size','Combined step size','MAML-based step size'},'Interpreter','latex')
grid on;
hold off

smooth_num = 2000;
subplot(2,1,2);
plot((0:length(Er1)-1)*(1/fs),smooth(20*log10(abs(Er1(1:end))),smooth_num),'Color',cc(2,:))
hold on
plot((0:length(Er2)-1)*(1/fs),smooth(20*log10(abs(Er2(1:end))),smooth_num),'LineStyle','--','Color',cc(3,:))
plot((0:length(Er3)-1)*(1/fs),smooth(20*log10(abs(Er3(1:end))),smooth_num),'Color',cc(9,:))
plot((0:length(Er4)-1)*(1/fs),smooth(20*log10(abs(Er4(1:end))),smooth_num),'LineStyle','--','Color',cc(5,:))
plot((0:length(Er5)-1)*(1/fs),smooth(20*log10(abs(Er5(1:end))),smooth_num),'Color',cc(7,:))
ylim([-45,15]);xlim([0.1,5])
title('(b) Mean square error ','Interpreter','latex')
xlabel('Time (seconds)','Interpreter','latex')
ylabel('MSE(dB)','Interpreter','latex')
grid on;
hold off


