%% example.m
% Author: Yongjie Zhuang
% Example demonstrating preconditioning filter design
% based on the simple example described in:
% "Causal preconditioning filters design for real-time multichannel ANC"

clear; close all; clc;
% --- Add src to MATLAB path ---
addpath(fullfile(fileparts(mfilename('fullpath')),'..','src'));

%% Parameters
% the following parameters are defined as the same in the paper
Nr = 2; % the number of reference microphones
Ne = 2; % the number of error microphones
Ns = 2; % the number of control speakers
N_Finv   = 10;  % length of inverse whitening filter
N_mininv = 10;  % length of inverse minimum-phase filter
plot_flag = 1;  % enable plots

fs = 1000;  % sampling frequency is 1 kHz
c0 = 340;   % sound speed is 340 m/s
M = [0.75, 0.3; 0.3, 1];    % the mixing matrix, which is the Fxx for precondition filter design function
l1 = 2.0;   % l1 is defined as 2 meter
l2 = 3.0;   % l2 is defined as 2 meter
N1 = round(l1*fs/c0);   % the delay N1, which should be 6
N2 = round(l2*fs/c0);   % the delay N2, which should be 9

%% construct system matrics based on parameters
% Fxx is such that reference signal = Fxx * v where v is white noise vector
% in this simple example, Fxx is the mixing matrix M
% but we need to reformat it to the correct tensor form
% the filter length for Fxx is set to be 16 (we only need 1 in this case, but just use a
% longer array to be more general)
Fxx = zeros(16,Nr,Nr);
Fxx(1,:,:) = M;

% for h_Ge, it is a 2 by 2 matrix
% for the filter length, we only need to be larger than max(N1,N2)
% let's pick either 16 or max(N1,N2)
N_Ge = max(16,max(N1,N2));
h_Ge = zeros(N_Ge,Ne,Ns);
% those are defined in the paper
h_Ge(N1,1,1) = 1;
h_Ge(N2,1,2) = l1/l2;
h_Ge(N2,2,1) = l1/l2;
h_Ge(N1,2,2) = 1;

%% --- Compute preconditioning filters ---
[Fxx_inv, Ge_min_inv, Ge_all] = precond_obtain_filter(Fxx,h_Ge,N_Finv,N_mininv,plot_flag);

disp(Fxx_inv)
disp(Ge_min_inv)
disp(Ge_all)
disp('Preconditioning filter computation complete.');
disp('Check plots for magnitude/phase and impulse responses.');