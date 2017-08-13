function [out] = RunMoxModel5(icond, init)
%RunMoxModel5
% Methane Oxidation Model, closed system
%
% Pools
%   1 = CH4, outside of cell
%   2 = "CH3OH" including all downstream
%
% Reactions
%   (1) 1 -> 2
%
% Species
% 61, 62, 63, 64, 65 = 12CH4, 13CH4, 12CH3D, 13CH3D, 12CH2D2
%
% All reactions are first-order in reactant.

% last modified 21 Mar 2017 // MATLAB R2012b (8.0.0.783) 64-bit Windows 7
% David T. Wang (dtw@mit.edu)

% clear all; close all

R_PDB = 0.0111802;
R_SMOW = 0.00015576; 

%% Load Conditions
% Rows = diffusion, AeOM, AOM, OH, Cl
% Columns = a13	aD	g64	g65
conditions = csvread('./oxidation/conditions.csv',1,2);    

in = conditions(icond,:)

%% SETUP

dCi = init(1)/1000;
dDi = init(2)/1000;
D64i = init(3)/1000;
D65i = init(4)/1000;

aC = in(1);
aD = in(2);
g64 = in(3);
g65 = in(4);

k = 0.1;

N61i = 1e3;

MassMatinit = [N61i, N61i*R_PDB*(dCi+1), N61i*R_SMOW*4*(dDi+1), N61i*R_PDB*R_SMOW*4*(dCi+1)*(dDi+1)*(D64i+1), N61i*6*(R_SMOW*(dDi+1))^2*(D65i+1); 
    0, 0, 0, 0, 0];
MassMat = MassMatinit;
% MassMat, npool x nspecies
%           SPECIES
%  POOL  61  62  63  64  65
%   1
%   2

kmat = [k, k*aC, k*aD, k*aC*aD*g64, k*aD^2*g65];
% kmat, nrxn x nspecies
%           SPECIES
% RXN   61  62  63  64  65
%  1

rxtsSPECmat = [1 0];
prodSPECmat = [0 1];
% ****SPECmat, nrxn x npool
%         POOL
% RXN   1   2 
%  1

%% Do reaction

tsteps = 1e3;

MMout = zeros(size(MassMat,1), size(MassMat,2), tsteps);
% MMout, npool x nspecies x tsteps

for i = 1:tsteps,
    % (1) "REACT", calculate masses reacted
    % (2) Calculate MassMat after "REACT" step
    % (3) Speciate products
    % (4) Calculate MassMat after speciating products
    MMrxtd = bsxfun(@times, MassMat, ...        
        bsxfun(@times, permute(kmat,[3 2 1]), ...
            permute(rxtsSPECmat,[2 3 1])));     % MMrxtd, npool x nspecies x nrxn
    MassMat = MassMat - sum(MMrxtd,3);          % MassMat, npool x nspecies
    MMprod = bsxfun(@times,sum(MMrxtd,1), ...   
        permute(prodSPECmat,[2 3 1]));          % MMprod, npool x nspecies x nrxn
    MassMat = MassMat + sum(MMprod,3);
    MMout(:,:,i) = MassMat;
end

MMout = cat(3,MassMatinit,MMout);

%% Reduce output

N61mat = squeeze(permute(MMout(:,1,:),[3 1 2]));    % tsteps x npool
R13mat = squeeze(permute(MMout(:,2,:)./MMout(:,1,:),[3 1 2]));
RDmat  = squeeze(permute(MMout(:,3,:)./MMout(:,1,:),[3 1 2]));
R13Dmat= squeeze(permute(MMout(:,4,:)./MMout(:,1,:),[3 1 2]));
RDDmat = squeeze(permute(MMout(:,5,:)./MMout(:,1,:),[3 1 2]));

f = N61mat(:,1)./N61i;
dCmat = 1000*(R13mat./R_PDB - 1);
dDmat  = 1000*(RDmat./R_SMOW/4 - 1);
D64mat = 1000*(R13Dmat./R13mat./RDmat - 1);
D65mat = 1000*(RDDmat./(6/16*RDmat.^2) - 1);

resid = [dCmat(:,1), dDmat(:,1), D64mat(:,1), D65mat(:,1)];
fi = [100:-1:1]';

for ii = 1:size(resid,2)
    residi = interp1(f*100,resid,fi);
end
    
out = [fi residi];
