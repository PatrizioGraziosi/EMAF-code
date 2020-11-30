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
% [1] P. Graziosi, C. Kumarasinghe, N. Neophytou, ACS Appl. Energy Mater. 
% 3, 6, 5913 (2020) ; https://pubs.acs.org/doi/abs/10.1021/acsaem.0c00825 %
%                                                                         %
% ----------------------------------------------------------------------- %

% equivalent effective mass finder from the bands
% 1) EF far from the edge --> carrier dentiy + T --> Nv --> m*
% 2) thermal injection velocity with non-SC code

% it reqiroes as input the 4D E(k) for the equivalent DOS mass
% also the three k-matrixes, i.e. k points ccoordinates, for the
% conductivity mass

clearvars
load('code\EMAF_input.mat');

%---- % Physical Constants ---------------------------------------
%-----------------------------------------------------------------
m0=9.1094e-31;            % [kg]
hbar=(6.6261e-34)/(2*pi); % [J-sec]
q0=1.609e-19;             % [col]
% T=300;                    % [K]
kB=1.38e-23;              % [J/K]
kT_q = kB*T/q0;           % [eV]   
%------------------------------------------------------------------------

% ----------------------------- code -------------------------------------
% ------------------------------------------------------------------------

tic
% load the file with the E(k) and the k-arrays

Load_filename=['..\Bands\Ek_',material_name,'.mat'];   % this file contains the Ek and the k arrays
load(Load_filename);




Ek_input = Ek;


dimensionality = '3D';
fsize = 40;
lwbor = 4;
lwpl  = 4;
if strcmp(dimensionality,'3D')
 % we check if the k_matrixes have been provided, because especially for
    % numerically built bandstructures, only the k axes are given. In this
    % case, we assume them to be perpendicular
    if exist('kx_matrix','var') == 0 
        disp('composing k_matrixes as orthogonal axes');
      for id_z = size(kz_array,2):-1:1
        for id_y = size(ky_array,2):-1:1
            for id_x = size(kx_array,2):-1:1
              kx_matrix(id_x,id_y,id_z) = kx_array(id_x);
              ky_matrix(id_x,id_y,id_z) = ky_array(id_y);
              kz_matrix(id_x,id_y,id_z) = kz_array(id_z);
            end
        end
      end  
    end
     

    % units of the k-arrays into S.I.
    if strcmp(k_units_pi_a,'yes')
        kx_matrix=kx_matrix/1e-9; ky_matrix=ky_matrix/1e-9; kz_matrix=kz_matrix/1e-9;
    end


    dkx_v = [kx_matrix(2,1,1) ky_matrix(2,1,1) kz_matrix(2,1,1)] - [kx_matrix(1,1,1) ky_matrix(1,1,1) kz_matrix(1,1,1)];
    dky_v = [kx_matrix(1,2,1) ky_matrix(1,2,1) kz_matrix(1,2,1)] - [kx_matrix(1,1,1) ky_matrix(1,1,1) kz_matrix(1,1,1)];
    dkz_v = [kx_matrix(1,1,2) ky_matrix(1,1,2) kz_matrix(1,1,2)] - [kx_matrix(1,1,1) ky_matrix(1,1,1) kz_matrix(1,1,1)];
    dVk = dot(cross(dkx_v,dky_v),dkz_v);
    
elseif strcmp(dimensionality,'2D')    
elseif strcmp(dimensionality,'1D')
end


if max(max(max(max(Ek)))) > 0
 
    carriers = 'electrons';
    
    if strcmp(dimensionality,'3D')
        
        %this shifts the conduction bands edge to zero
        [pseudoconduction_bands, pseudovalence_bands, Ek] = shifting_bands_EMAF_funct(Ek) ;
        Ek = Ek(:,:,:,pseudoconduction_bands(1):pseudoconduction_bands(size(pseudoconduction_bands,2)));
        
        % run the code to calculate the DOS equivalent effective mass
        m_equivalent_dos = EMAF_DOS_funct(T, Spin_Orbit_flag, dVk, Ek) ;
      
        % run the code to calculate the conductivtiy equivalent effective
        % mass and eventually plot the results
        nonSC_charge_calc_EMAF_3D_script
        
        %
        [nv, k_dist] = valleys_parameter_EMAF_funct (Ek, kx_matrix, ky_matrix, kz_matrix, pseudoconduction_bands, T) ;
        
         nv_CB = nv ;
         k_dist_CB = k_dist ;
         vp_CB = nv * k_dist ;
    
    elseif strcmp(dimensionality,'2D')
    elseif strcmp(dimensionality,'1D')
    
    end

    equivalent_DOS_mass_electrons = m_equivalent_dos;
    equivalent_injection_mass_electrons = m_equivalent_injection;
    disp('DOS effective mass for the electrons') 
    disp(equivalent_DOS_mass_electrons)
    disp('conductivity effective mass for the electrons') 
    disp(equivalent_injection_mass_electrons)
    disp('weighted number of valleys for the conduction band') 
    disp(nv_CB)
    disp('valleys parameter for the conduction band') 
    disp(vp_CB)
    
    me_inj_x = m_inj_x;
    me_inj_y = m_inj_y;
    me_inj_z = m_inj_z;

end




if min(min(min(min(Ek_input)))) < 0
    
    carriers = 'holes';
    % reversing the bands to make with the holes as with the electrons
    clear Ek
    Ek = -Ek_input;
    if strcmp(dimensionality,'3D')
   
        %this shifts the conduction bands edge to zero
        clear pseudoconduction_bands
        [pseudoconduction_bands, pseudovalence_bands, Ek] = shifting_bands_EMAF_funct(Ek) ;
        Ek = Ek(:,:,:,pseudoconduction_bands(1):pseudoconduction_bands(size(pseudoconduction_bands,2)));
        % run the code to calculate the DOS equivalent effective mass
        m_equivalent_dos = EMAF_DOS_funct(T, Spin_Orbit_flag, dVk, Ek) ;
        
        % run the code to calculate the conductivtiy equivalent effective
        % mass and eventually plot the results
        nonSC_charge_calc_EMAF_3D_script
        
        %
        %
        [nv, k_dist] = valleys_parameter_EMAF_funct (Ek, kx_matrix, ky_matrix, kz_matrix, pseudoconduction_bands, T) ;
        
         nv_VB = nv ;
         k_dist_VB = k_dist ;
         vp_VB = nv * k_dist ;
         
         
    elseif strcmp(dimensionality,'2D')    
    elseif strcmp(dimensionality,'1D')    
    end

    equivalent_DOS_mass_holes = m_equivalent_dos;
    equivalent_injection_mass_holes = m_equivalent_injection;
    disp('DOS effective mass for the holes') 
    disp(equivalent_DOS_mass_holes)
    disp('conductivity effective mass for the holes') 
    disp(equivalent_injection_mass_holes)
    disp('weighted number of valleys for the valence band') 
    disp(nv_VB)
    disp('valleys parameter for the valence band') 
    disp(vp_VB)
    
    mh_inj_x = m_inj_x;
    mh_inj_y = m_inj_y;
    mh_inj_z = m_inj_z;
    
end
Ek = Ek_input;

clearvars -except Ek material_name equivalent_DOS_mass_electrons equivalent_DOS_mass_holes equivalent_injection_mass_electrons equivalent_injection_mass_holes nv_VB nv_CB vp_CB vp_VB k_dist_CB k_dist_VB T me_inj_x me_inj_y me_inj_z mh_inj_x mh_inj_y mh_inj_z
save_filename=['..\parameters_',material_name,'.mat'];   % this file contains the Ek and the k arrays
save(save_filename);
toc