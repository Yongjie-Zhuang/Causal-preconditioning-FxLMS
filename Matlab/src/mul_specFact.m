%% mul_specFact: Multichannel spectral factorization
% Author: Yongjie Zhuang
% Sxx: Power spectrum density Sxx(z) expressed as a polynomial matrix of
% order q, dimension is Nx*Nx*(q+1) (indexed from z^0 to z^-q)
% L and Re: result of spectral factorization, such that Sxx = L*Re*L';

function [F,L,Re]= mul_specFact(Sxx)
%% obtain dimension and check input data
[Nq,Nx,Nx2] = size(Sxx);
Nq = Nq-1;
comp1 = [Nx];
comp2 = [Nx2];
if any(comp1~=comp2)
    error('The dimension of input data does not match');
end
%% Initialization
F = kron(diag(ones(1,Nq-1),-1),eye(Nx));
h = kron([zeros(1,Nq-1),1],eye(Nx));
N_bar = zeros(Nx*Nq,Nx);
for ii = 1:Nq
    N_bar((ii-1)*Nx+1:ii*Nx,:) = squeeze(Sxx(Nq+2-ii,:,:));
end
%% solve the Discrete algebraic Riccati equation
R0 = squeeze(Sxx(1,:,:));
[Sigma,~,~] = idare(F',h',zeros(Nx*Nq),-R0,-N_bar,eye(Nx*Nq));
Re = R0 - h*Sigma*h';
g = (N_bar - F*Sigma*h')/Re;
%% reshape to get L
L = zeros(Nx,Nx,Nq);
for ii = 1:Nq
    L(:,:,Nq+1-ii) = g((ii-1)*Nx+1:ii*Nx,:);
end
% get filter F
chol_Re = chol(Re);
chol_Re = chol_Re';
F = zeros(Nq+1,Nx,Nx);
F(1,:,:) = chol_Re;
for ii = 2:Nq+1
    F(ii,:,:) = L(:,:,ii-1)*chol_Re;
end
end