% Copyright 2020 Patrizio Graziosi, Neophytos Neophytou                   %
% Written and developed by                                                %
% Patrizio Graziosi and Neophytos Neophytou                               %
% patrizio.graziosi@cnr.it and n.neophytou@warwick.ac.uk                  %
% under the                                                               %
% Marie Curie - Individual Fellowships -  GENESIS - project ID 788465     %
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


function [pseudoconduction_bands, pseudovalence_bands, Ek] = shifting_bands_EMAF_funct(Ek) %#codegen
% shifting of the bands to zero, that have to be already flipped in the
% case of holes, this is why I say pseudoconduction

    Emaximalextreme = zeros(1,size(Ek,4)) ;
    Eminimalextreme = zeros(1,size(Ek,4)) ;
    for id_n = 1:size(Ek,4)
        Emaximalextreme(id_n) = max(max(max(Ek(:,:,:,id_n))));
        Eminimalextreme(id_n) = min(min(min(Ek(:,:,:,id_n))));
    end

    [~,pseudoconduction_bands] = find( abs(Emaximalextreme) > abs(Eminimalextreme) );
    [~,pseudovalence_bands] = find( abs(Emaximalextreme) < abs(Eminimalextreme) );


    new_Fermi = min(min(min(min(Ek(:,:,:,pseudoconduction_bands(1):pseudoconduction_bands(size(pseudoconduction_bands,2)))))));
    Ek = Ek - new_Fermi; % this the conduction band, or the flipped valence band, start from zero, positive values of the EF_array are into the band and negative EF values are into the gap
end