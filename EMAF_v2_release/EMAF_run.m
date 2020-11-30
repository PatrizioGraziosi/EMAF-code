% Copyright 2020 Patrizio Graziosi, Neophytos Neophytou                   %
% Written and developed by                                                %
% Patrizio Graziosi, patrizio.graziosi@cnr.it, during his                 %
% Marie Curie - Individual Fellowships  GENESIS - project ID 788465       %
% Generic transport simulator for new generation thermoelectric materials %
% under the supervision of Neophytos Neophytou                            %
% ----------------------------------------------------------------------- %
% This file is distributed under the terms of the GNU                     %
% General Public License. See the file `LICENSE' in  the root directory   %
% of the present distribution.                                            %
% ----------------------------------------------------------------------- %
%                                                                         %
% Please cite the code source and the paper [1] when publishing results   %
% obtained  using the present EMAF code                                   %                      
% [1] P. Graziosi, C. Kumarasinghe, N. Neophytou, ACS Appl. Energy Mater. %
% 3, 6, 5913 (2020) ; https://pubs.acs.org/doi/abs/10.1021/acsaem.0c00825 %
%                                                                         %
% ----------------------------------------------------------------------- %


% equivalent effective mass finder from the bands
% 1) EF far from the edge --> carrier dentiy + T --> Nv --> m*
% 2) thermal injection velocity iwth non-SC code

% it reqiroes as input the 4D E(k) for the equivalent DOS mass
% also the three k-arrays for the conductivity DOS

clearvars
% ---------------------- input instructions ------------------------------
% ------------------------------------------------------------------------
material_name = 'HfCoSb' ; %    'HfNiSn'; % 'Mg3Sb2' ; %
T = 300; % K

k_units_SI = 'yes'; % in 1/m
k_units_pi_a = 'no'; %  for numerically built bands 

Spin_Orbit_flag = 'off'; % if off, the two spins are equivalent, i.e. each band contains two electrons

plot_velocities = 'yes'; % to plot automatically the injection velocities
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

savefilename = ['code\EMAF_input.mat'];
save(savefilename)
runfilename = ['code\main_EMAF.m'];
run(runfilename)
