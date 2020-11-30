% Copyright 2020 Patrizio Graziosi, Neophytos Neophytou                   %
% Written and developed by                                                %
% Patrizio Graziosi, patrizio.graziosi@cnr.it, during the                 %  
% Marie Curie - Individual Fellowships  GENESIS - project ID 788465       %
% Generic transport simulator for new generation thermoelectric materials %
% ----------------------------------------------------------------------- %
% This file is distributed under the terms of the GNU                     %
% General Public License. See the file `LICENSE' in  the root directory   %
% of the present distribution.                                            %
% ----------------------------------------------------------------------- %
%                                                                         %
% Please cite the code source and the paper [1] when publishing results   %
% obtained  using the present EMAF code                                   %                             
% [1] P. Graziosi, C. Kumarasinghe, N. Neophytou, ACS Appl. Energy Mater. 
% 3, 6, 5913 (2020) ; https://pubs.acs.org/doi/abs/10.1021/acsaem.0c00825 %
%                                                                         %
% ----------------------------------------------------------------------- %

function m_equivalent_dos = EMAF_DOS_funct(T, Spin_Orbit_flag, dVk, Ek)

q0=1.609e-19;             % [col]
kB=1.38e-23;              % [J/K]

T_emaf = T;

EF_emaf_dos = -0.6; % for holes all energies have been reversed

N_imp_emaf_dos = (strcmp(Spin_Orbit_flag,'off')+1)/(2*pi)^3*dVk*sum(sum(sum(sum(1./(exp((Ek-EF_emaf_dos)./(kB*T/q0))+1))))); % in m^-3

 E_edge = min(min(min(min(Ek)))); % Ek is already separated into CB and VB in the main code
% the band edge has been put to zero
Nband = N_imp_emaf_dos*1e-6/exp(-abs(E_edge-EF_emaf_dos)/(kB*T_emaf/q0)); % Nv or Nc in parabolic band approx. 

m_equivalent_dos = (Nband/4.84e15)^(2/3)/T_emaf; % M. Guzzi, Principi di fisica dei semiconduttori or R. F. Pierret, Advanced Semiconductor Fundamentals

% in Si conduction bands gives 1.0597 rather than 1.0618, 
% in Si valence band gives 0.854
% in GaAs conduction band gives 0.0629 rather than 0.063
% I'm really happy with this