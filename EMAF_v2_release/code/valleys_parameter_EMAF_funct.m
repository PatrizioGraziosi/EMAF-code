% Copyright 2020 Patrizio Graziosi                                        %
% A creation of Patrizio Graziosi, written and developed by               %
% Patrizio Graziosi, patrizio.graziosi@cnr.it                             %  
% ----------------------------------------------------------------------- %
% This file is distributed under the terms of the GNU                     %
% General Public License.                                                 %
% ----------------------------------------------------------------------- %
%                                                                         %
% Please cite the code source and the paper [1] when publishing results   %
% obtained  using the present  code                                       %
% [1] P. Graziosi, C. Kumarasinghe, N. Neophytou, ACS Appl. Energy Mater. %
% 3, 6, 5913 (2020) ; https://pubs.acs.org/doi/abs/10.1021/acsaem.0c00825 %
%                                                                         %
% ----------------------------------------------------------------------- %

% ---------------------- input instructions ------------------------------
% ------------------------------------------------------------------------

function [nv, k_dist] = valleys_parameter_EMAF_funct (Ek, kx_matrix, ky_matrix, kz_matrix, pseudoconduction_bands, T)

%--------------------- Physical Constants ------------------------
%-----------------------------------------------------------------
q0=1.609e-19;             % [col]
kB=1.38e-23;              % [J/K] 
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

Valley_threshold = 5 * kB*T ;

% Ek = Ek(:,:,:,pseudoconduction_bands(1):pseudoconduction_bands(size(pseudoconduction_bands,2)));


% -----------------------  valleys count --------------------------------
Ev_CB = struct();
for id_band = size(pseudoconduction_bands,2):-1:1
    Ev_CB(id_band).E = [];
end

for id_band = size(pseudoconduction_bands,2):-1:1

    Ek_temp = Ek(:,:,:,id_band);
    
    nv_temp = size( find( Ek_temp( imregionalmin(Ek_temp,26) ) < Valley_threshold ) );
    
    if not(nv_temp(1) == 0)
        E_valleys_temp =  ( Ek_temp( imregionalmin(Ek_temp,26) ) ) ;
        Ev_CB(id_band).E  = E_valleys_temp;
        valley_weight =  ( 1 ./ ( exp( ( E_valleys_temp - 0 ) ./ (kB*T/q0) ) + 1 ) );
%         weighted_nv_CB(id_band) = nv_temp(1) * mean(valley_weight);
        weighted_nv_CB(id_band) = sum(valley_weight);
    end

end
[~,n] = find( weighted_nv_CB > 0 ) ;
nv = mean(weighted_nv_CB(n) ) ;

% ------------------------------------------------------------------------



% -------------------- valleys distance ----------------------------------
for id_band = size(pseudoconduction_bands,2):-1:1

    Ek_temp = Ek(:,:,:,id_band) ;
    
    id_k = imregionalmin(Ek_temp,26) ;
    E_minima = Ek_temp(id_k) ;
    for id_E = size(E_minima,1) : -1 : 1 
        if E_minima(id_E) < Valley_threshold

            kx_min = [];
            ky_min = [];
            kz_min = [];
            for i = 1:size(find(id_k),1)
                kx_min(i) = kx_matrix(i);
                ky_min(i) = ky_matrix(i);
                kz_min(i) = kz_matrix(i);
            end

            for i = 1:size(find(id_k),1)
                delta_array = (kx_min(i)-kx_min).^2 + (ky_min(i)-ky_min).^2 + (kz_min(i)-kz_min).^2 ;
                [~,n] = find( delta_array == 0 ) ;
                delta_array(n) = [] ;
                delta(i) = mean( sqrt( delta_array ) ) ;
            end

            delta_band(id_band) = mean(delta) ;
            
        end
    end

end
[~,n] = find( delta_band > 0 ) ;

k_dist = mean( delta_band(n) ) *1e-9 ;

