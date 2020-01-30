%%% Main File for pre-processing and calling the piMF function to recover
%%% reaction networks and the extents of reaction.
%%% An example is used to demonstrate the algorithm. Data_Example.mat
%%% contains time series data and atomic matrix for the different chemicals
%%% involved in the system as follows:
%%% time: 50 x 1 dimentional matrix of time stamps
%%% C   : 50 x 7 dimentional matrix of concentration for seven species
%%% A   : 4  x 7 dimentional matrix of atomic matrix
clear all
clc
close all
%%% Load data
load Data_Example
%%% size of data (observations (m_C) x number of species (n_C))
[m_C,n_C]=size(C);
%%% Computing Reaction variant form of C
C_RV = (C-ones(m_C,1)*C(1,:));
%%% Maximum number of reactions
R_max=n_C-rank(A);
%%% Observed number of reactions through singular value decomposition
s_obs=svd(C_RV);
%%% Number of significant eigenvalues is computed using 97% thersholding
%%% (s_thres)
s_thres = (sum(s_obs(1:3))/sum(s_obs))*100;
%%% s_thers is 97% for three singular values. 
%%%Note that the maximum number singular values that can be chosen is 3=R_max.

%%% Value that stoichiometric coefficient can take
S_set=[-1, 0, 1];
%%% Number of reactions in the system
R= R_max;
% Sparsity: Maximum number of species that can participate
K=4;
%%% Search diameter for sphere decoding algorithm
dia_sphere_decoding = 0.5;

%%%% Initialization of algorithm
W_initial = sort(abs(randn(50,3)), 'ascend')
%%%% Calling piMF 
[X_out,N_out,fvalopt]=piMF(C_RV,S_set,R,A,K,dia_sphere_decoding,W_initial);

plot(time,X_out,'o')
N_out
