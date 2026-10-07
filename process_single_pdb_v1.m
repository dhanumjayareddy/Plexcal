function diffusion_length = process_single_pdb(pdb_filename, microtubule_diameter_angstrom)
    % Process a single PDB file and return diffusion length
    % Inputs:
    %   pdb_filename - path to PDB file
    %   microtubule_diameter_angstrom - diameter in Angstroms
    % Output:
    %   diffusion_length - calculated diffusion length in Angstroms
    
    try
        % EXACT from original code - start
        varsbefore = who; %// get names of current variables (note 1)
        
        % create a cylinder with a defined height and radius.
        Microtubule_diameter_Angstrom = microtubule_diameter_angstrom ; 
        Microtubule_radius_Angstrom = Microtubule_diameter_Angstrom./2 ;
        
        % Read PDB file
        Tubulin_load = pdbread(pdb_filename);
        
        % Extract atom information
        Locations = Tubulin_load.Model.Atom;
        Length_Locations = length(Locations);
        
        Residues_Name = string({Locations.resName})';
        Residues_Number = double(string({Locations.resSeq})');
        
        Atom_Names    = string({Locations.AtomName})' ;   % <-- atom names ("CG","CD1",...) 
        
        Residues_X = double(string([Locations.X]))';
        Residues_Y = double(string([Locations.Y]))';
        Residues_Z = double(string([Locations.Z]))';
        
        % Initialize chromophore arrays
        X_tryptophan = zeros(Length_Locations, 1);
        Y_tryptophan = zeros(Length_Locations, 1);
        Z_tryptophan = zeros(Length_Locations, 1);
        Residue_Number_of_trp = zeros(Length_Locations, 1);
        Names_tryptophan = strings(Length_Locations, 1);
        
        X_tyrosine = zeros(Length_Locations, 1);
        Y_tyrosine = zeros(Length_Locations, 1);
        Z_tyrosine = zeros(Length_Locations, 1);
        Residue_Number_of_tyr = zeros(Length_Locations, 1);
        Names_tyrosine = strings(Length_Locations, 1);
        
        X_phenylalanine = zeros(Length_Locations, 1);
        Y_phenylalanine = zeros(Length_Locations, 1);
        Z_phenylalanine = zeros(Length_Locations, 1);
        Residue_Number_of_phen = zeros(Length_Locations, 1);
        Names_phenylalanine = strings(Length_Locations, 1);
        
        % Initialize logical arrays
        Doyouhave_tryptophan = false(Length_Locations, 1);
        Doyouhave_tyrosine = false(Length_Locations, 1);
        Doyouhave_phenylalanine = false(Length_Locations, 1);
        
        for k = 1:Length_Locations
        
            Doyouhave_tryptophan(k) = (contains(Residues_Name(k),'TRP','IgnoreCase',true)) ;
            Numberfiles_tryptophan = nnz(Doyouhave_tryptophan) ;
            Find_tryptophan = find(Doyouhave_tryptophan,1) ;
            
            if Doyouhave_tryptophan(k) == 1
                X_tryptophan(k,:) = Residues_X(k) ; 
                Y_tryptophan(k,:) = Residues_Y(k) ; 
                Z_tryptophan(k,:) = Residues_Z(k) ; 
                Residue_Number_of_trp(k,:) = Residues_Number(k) ;
                Names_tryptophan(k,:) = Atom_Names(k) ;   % not Residues_Name(k)
            end
            
            
            Doyouhave_tyrosine(k) = (contains(Residues_Name(k),'TYR','IgnoreCase',true)) ;
            Numberfiles_tyrosine = nnz(Doyouhave_tyrosine) ;
            Find_tyrosine = find(Doyouhave_tyrosine,1) ;
            
            if Doyouhave_tyrosine(k) == 1
                X_tyrosine(k,:) = Residues_X(k) ; 
                Y_tyrosine(k,:) = Residues_Y(k) ; 
                Z_tyrosine(k,:) = Residues_Z(k) ; 
                Residue_Number_of_tyr(k,:) = Residues_Number(k) ;
                Names_tyrosine(k,:) = Atom_Names(k) ;   % not Residues_Name(k)
            end
            
            
            Doyouhave_phenylalanine(k) = (contains(Residues_Name(k),'PHE','IgnoreCase',true)) ;
            Numberfiles_phenylalanine = nnz(Doyouhave_phenylalanine) ;
            Find_phenylalanine = find(Doyouhave_phenylalanine,1) ;
            
            if Doyouhave_phenylalanine(k) == 1
                X_phenylalanine(k,:) = Residues_X(k) ; 
                Y_phenylalanine(k,:) = Residues_Y(k) ; 
                Z_phenylalanine(k,:) = Residues_Z(k) ; 
                Residue_Number_of_phen(k,:) = Residues_Number(k) ;
                Names_phenylalanine(k,:) = Atom_Names(k) ;   % not Residues_Name(k)
            end
            
        end
        
        
        % --- Which rows are actual TRP/TYR/PHE atoms? ---
        mask_trp  = ismissing(Names_tryptophan)    | Names_tryptophan    == "" ;
        mask_tyr  = ismissing(Names_tyrosine)      | Names_tyrosine      == "" ;
        mask_phen = ismissing(Names_phenylalanine) | Names_phenylalanine == "" ;
        
        % --- Atom counts from the masks (rows that are NOT empty) ---
        Number_of_Atoms_of_tryptophan    = nnz(~mask_trp) ;
        Number_of_Atoms_of_tyrosine      = nnz(~mask_tyr) ;
        Number_of_Atoms_of_phenylalanine = nnz(~mask_phen) ;
        
        % --- Remove the non-atom rows ---
        X_tryptophan(mask_trp,:)    = [] ;
        Y_tryptophan(mask_trp,:)    = [] ;
        Z_tryptophan(mask_trp,:)    = [] ;
        Names_tryptophan(mask_trp,:) = [] ;
        Residue_Number_of_trp(mask_trp,:)  = [] ;
        
        X_tyrosine(mask_tyr,:)    = [] ;
        Y_tyrosine(mask_tyr,:)    = [] ;
        Z_tyrosine(mask_tyr,:)    = [] ;
        Names_tyrosine(mask_tyr,:) = [] ;
        Residue_Number_of_tyr(mask_tyr,:)  = [] ;
        
        X_phenylalanine(mask_phen,:)    = [] ;
        Y_phenylalanine(mask_phen,:)    = [] ;
        Z_phenylalanine(mask_phen,:)    = [] ;
        Names_phenylalanine(mask_phen,:) = [] ;
        Residue_Number_of_phen(mask_phen,:) = [] ;
        
        
        % --- Residue counts from the trimmed residue-number vectors ---
        Number_of_tryptophan    = numel(unique(Residue_Number_of_trp)) ;
        Number_of_tyrosine      = numel(unique(Residue_Number_of_tyr)) ;
        Number_of_phenylalanine = numel(unique(Residue_Number_of_phen)) ;
        Number_of_molecules     = Number_of_tryptophan + Number_of_tyrosine + Number_of_phenylalanine ;
        
        if Number_of_molecules < 2
           diffusion_length = 0;
           return;
        end
        
        
        % ---- Group TRP atoms by residue number ----
        
        if Number_of_tryptophan > 0
        
           unique_trp = unique(Residue_Number_of_trp) ;
           nTRP = numel(unique_trp) ;
           maxA_trp = 0 ;
           
           for r = 1:nTRP
               maxA_trp = max(maxA_trp, sum(Residue_Number_of_trp == unique_trp(r))) ;
           end
           
           x_coord_tryptophan    = nan(maxA_trp, nTRP) ;
           y_coord_tryptophan    = nan(maxA_trp, nTRP) ;
           z_coord_tryptophan    = nan(maxA_trp, nTRP) ;
           name_coord_tryptophan = strings(maxA_trp, nTRP) ;
           
           for r = 1:nTRP
               idx = find(Residue_Number_of_trp == unique_trp(r)) ;
               n = numel(idx) ;
               x_coord_tryptophan(1:n,r)    = X_tryptophan(idx) ;
               y_coord_tryptophan(1:n,r)    = Y_tryptophan(idx) ;
               z_coord_tryptophan(1:n,r)    = Z_tryptophan(idx) ;
               name_coord_tryptophan(1:n,r) = Names_tryptophan(idx) ;
           end
           
        else
           x_coord_tryptophan    = nan(0,0) ;
           y_coord_tryptophan    = nan(0,0) ;
           z_coord_tryptophan    = nan(0,0) ;
           name_coord_tryptophan = strings(0,0) ;
        end
        
        
        % ---- Group TYR atoms by residue number ----
        
        if Number_of_tyrosine > 0
        
           unique_tyr = unique(Residue_Number_of_tyr) ;
           nTYR = numel(unique_tyr) ;
           maxA_tyr = 0 ;
           
           for r = 1:nTYR
               maxA_tyr = max(maxA_tyr, sum(Residue_Number_of_tyr == unique_tyr(r))) ;
           end
           
           x_coord_tyrosine    = nan(maxA_tyr, nTYR) ;
           y_coord_tyrosine    = nan(maxA_tyr, nTYR) ;
           z_coord_tyrosine    = nan(maxA_tyr, nTYR) ;
           name_coord_tyrosine = strings(maxA_tyr, nTYR) ;
           
           for r = 1:nTYR
               idx = find(Residue_Number_of_tyr == unique_tyr(r)) ;
               n = numel(idx) ;
               x_coord_tyrosine(1:n,r)    = X_tyrosine(idx) ;
               y_coord_tyrosine(1:n,r)    = Y_tyrosine(idx) ;
               z_coord_tyrosine(1:n,r)    = Z_tyrosine(idx) ;
               name_coord_tyrosine(1:n,r) = Names_tyrosine(idx) ;
           end
           
        else
           x_coord_tyrosine    = nan(0,0) ;
           y_coord_tyrosine    = nan(0,0) ;
           z_coord_tyrosine    = nan(0,0) ;
           name_coord_tyrosine = strings(0,0) ;
        end
        
        
        % ---- Group PHE atoms by residue number ----
        
        if Number_of_phenylalanine > 0
        
           unique_phe = unique(Residue_Number_of_phen) ;
           nPHE = numel(unique_phe) ;
           maxA_phe = 0 ;
           
           for r = 1:nPHE
               maxA_phe = max(maxA_phe, sum(Residue_Number_of_phen == unique_phe(r))) ;
           end
           
           x_coord_phenylalanine    = nan(maxA_phe, nPHE) ;
           y_coord_phenylalanine    = nan(maxA_phe, nPHE) ;
           z_coord_phenylalanine    = nan(maxA_phe, nPHE) ;
           name_coord_phenylalanine = strings(maxA_phe, nPHE) ;
           
           for r = 1:nPHE
               idx = find(Residue_Number_of_phen == unique_phe(r)) ;
               n = numel(idx) ;
               x_coord_phenylalanine(1:n,r)    = X_phenylalanine(idx) ;
               y_coord_phenylalanine(1:n,r)    = Y_phenylalanine(idx) ;
               z_coord_phenylalanine(1:n,r)    = Z_phenylalanine(idx) ;
               name_coord_phenylalanine(1:n,r) = Names_phenylalanine(idx) ;
           end
           
        else
           x_coord_phenylalanine    = nan(0,0) ;
           y_coord_phenylalanine    = nan(0,0) ;
           z_coord_phenylalanine    = nan(0,0) ;
           name_coord_phenylalanine = strings(0,0) ;
        end
        
        
        % --- Centroid of each residue, one column at a time ---
        Centroid_tryptophan    = nan(Number_of_tryptophan, 3) ;
        Centroid_tyrosine      = nan(Number_of_tyrosine, 3) ;
        Centroid_phenylalanine = nan(Number_of_phenylalanine, 3) ;
        
        for r = 1:Number_of_tryptophan
        
            Centroid_tryptophan(r,:) = [ mean(x_coord_tryptophan(:,r), 'omitnan'), mean(y_coord_tryptophan(:,r), 'omitnan'), mean(z_coord_tryptophan(:,r), 'omitnan') ] ;
            
        end
        
        
        for r = 1:Number_of_tyrosine
            
            Centroid_tyrosine(r,:) = [ mean(x_coord_tyrosine(:,r), 'omitnan'), mean(y_coord_tyrosine(:,r), 'omitnan'), mean(z_coord_tyrosine(:,r), 'omitnan') ] ;
        
        end
        
        
        for r = 1:Number_of_phenylalanine
        
            Centroid_phenylalanine(r,:) = [ mean(x_coord_phenylalanine(:,r), 'omitnan'), mean(y_coord_phenylalanine(:,r), 'omitnan'), mean(z_coord_phenylalanine(:,r), 'omitnan') ] ;
            
        end
        
        
        positions_xyz = vertcat(Centroid_tryptophan,Centroid_tyrosine,Centroid_phenylalanine);
        
        
        % Calculate distance matrix - EXACT from original ; PDB records coordinates in Angstroms only
        Distances_matrix = zeros(Number_of_molecules,Number_of_molecules);
        for First_molecule = 1:Number_of_molecules
            for Second_molecule = 1:Number_of_molecules
                Distances_matrix(First_molecule,Second_molecule) = pdist2(positions_xyz(First_molecule,:),positions_xyz(Second_molecule,:));
            end
        end
        
        % FRET parameters
        Debye_units = 3.33564*10^(-30);
        
        J_TRP_TRP = 5.5953e-06;
        J_TYR_TYR = 4.7933e-04;
        J_PHE_PHE = 6.1123e-04;
        J_TRP_TYR = 2.3247e-05;
        J_TYR_TRP = 9.2550e-04;
        J_PHE_TYR = 0.0043;
        J_PHE_TRP = 0.0039;
        J_TRP_PHE = 1.3034e-07;
        J_TYR_PHE = 1.8891e-05;
        
        mu_Trp = 2.074*Debye_units;
        mu_Tyr = 1.18*Debye_units;
        mu_Phen = 0.28*Debye_units;
        
        % Build J and mu matrices - EXACT from original code
        J_matrix = zeros(Number_of_molecules,Number_of_molecules);
        mu_squared_matrix = zeros(Number_of_molecules,Number_of_molecules);
        
        trp_locations_on_the_J_and_mu_matrix = 1:Number_of_tryptophan;
        tyr_locations_on_the_J_and_mu_matrix = (Number_of_tryptophan+1):(Number_of_tyrosine+Number_of_tryptophan);
        phen_locations_on_the_J_and_mu_matrix = (Number_of_tyrosine+Number_of_tryptophan+1):Number_of_molecules;
        
        % Cross-chromophore J values - EXACT from original
        % Same chromophore J values
        J_matrix(trp_locations_on_the_J_and_mu_matrix,trp_locations_on_the_J_and_mu_matrix) = J_TRP_TRP;
        J_matrix(tyr_locations_on_the_J_and_mu_matrix,tyr_locations_on_the_J_and_mu_matrix) = J_TYR_TYR;
        J_matrix(phen_locations_on_the_J_and_mu_matrix,phen_locations_on_the_J_and_mu_matrix) = J_PHE_PHE;
        
        % Cross-chromophore J values
        J_matrix(trp_locations_on_the_J_and_mu_matrix,tyr_locations_on_the_J_and_mu_matrix) = J_TRP_TYR;
        J_matrix(tyr_locations_on_the_J_and_mu_matrix,trp_locations_on_the_J_and_mu_matrix) = J_TYR_TRP;
        J_matrix(trp_locations_on_the_J_and_mu_matrix,phen_locations_on_the_J_and_mu_matrix) = J_TRP_PHE;
        J_matrix(phen_locations_on_the_J_and_mu_matrix,trp_locations_on_the_J_and_mu_matrix) = J_PHE_TRP;
        J_matrix(tyr_locations_on_the_J_and_mu_matrix,phen_locations_on_the_J_and_mu_matrix) = J_TYR_PHE;
        J_matrix(phen_locations_on_the_J_and_mu_matrix,tyr_locations_on_the_J_and_mu_matrix) = J_PHE_TYR;
        
        % Cross-chromophore mu values - EXACT from original
        % Same chromophore mu values
        mu_TRP_TRP = mu_Trp*mu_Trp;
        mu_TYR_TYR = mu_Tyr*mu_Tyr;
        mu_PHE_PHE = mu_Phen*mu_Phen;
        
        % Cross-chromophore mu values
        mu_TRP_TYR = mu_Trp*mu_Tyr;
        mu_TYR_PHE = mu_Tyr*mu_Phen;
        mu_TRP_PHE = mu_Trp*mu_Phen;
        
        % Fill mu matrix - EXACT from original
        mu_squared_matrix(trp_locations_on_the_J_and_mu_matrix,trp_locations_on_the_J_and_mu_matrix) = mu_TRP_TRP;
        mu_squared_matrix(tyr_locations_on_the_J_and_mu_matrix,tyr_locations_on_the_J_and_mu_matrix) = mu_TYR_TYR;
        mu_squared_matrix(phen_locations_on_the_J_and_mu_matrix,phen_locations_on_the_J_and_mu_matrix) = mu_PHE_PHE;
        
        % Cross-chromophore mu matrix entries
        mu_squared_matrix(trp_locations_on_the_J_and_mu_matrix,tyr_locations_on_the_J_and_mu_matrix) = mu_TRP_TYR;
        mu_squared_matrix(tyr_locations_on_the_J_and_mu_matrix,trp_locations_on_the_J_and_mu_matrix) = mu_TRP_TYR;
        mu_squared_matrix(trp_locations_on_the_J_and_mu_matrix,phen_locations_on_the_J_and_mu_matrix) = mu_TRP_PHE;
        mu_squared_matrix(phen_locations_on_the_J_and_mu_matrix,trp_locations_on_the_J_and_mu_matrix) = mu_TRP_PHE;
        mu_squared_matrix(tyr_locations_on_the_J_and_mu_matrix,phen_locations_on_the_J_and_mu_matrix) = mu_TYR_PHE;
        mu_squared_matrix(phen_locations_on_the_J_and_mu_matrix,tyr_locations_on_the_J_and_mu_matrix) = mu_TYR_PHE;
        
        % Calculate coupling constants - EXACT physics from original code
        refractive_index = 1.4 ; % protein refractive index.
        epsilon_0 = (8.854*10^-12) ; % units being in Farads per meter.
        constants = 1./(4*pi*(epsilon_0)) ;
        
        
        % this is the first step for concatenating the tryp x coordinates with those of tyr and phen. Here, we are making the number of molecules effectively the 'same' by adding NaN rows
        % EXACT kappa calculation from original code
        
        Number_of_rows_kappa = max([ size(x_coord_tryptophan,1), size(x_coord_tyrosine,1), size(x_coord_phenylalanine,1) ]) ;
        
        x_coord_tryptophan    = [ x_coord_tryptophan    ; nan(Number_of_rows_kappa - size(x_coord_tryptophan,1),    size(x_coord_tryptophan,2))    ] ;
        y_coord_tryptophan    = [ y_coord_tryptophan    ; nan(Number_of_rows_kappa - size(y_coord_tryptophan,1),    size(y_coord_tryptophan,2))    ] ;
        z_coord_tryptophan    = [ z_coord_tryptophan    ; nan(Number_of_rows_kappa - size(z_coord_tryptophan,1),    size(z_coord_tryptophan,2))    ] ;
        name_coord_tryptophan = [ name_coord_tryptophan ; strings(Number_of_rows_kappa - size(name_coord_tryptophan,1), size(name_coord_tryptophan,2)) ] ;
        
        x_coord_tyrosine      = [ x_coord_tyrosine      ; nan(Number_of_rows_kappa - size(x_coord_tyrosine,1),      size(x_coord_tyrosine,2))      ] ;
        y_coord_tyrosine      = [ y_coord_tyrosine      ; nan(Number_of_rows_kappa - size(y_coord_tyrosine,1),      size(y_coord_tyrosine,2))      ] ;
        z_coord_tyrosine      = [ z_coord_tyrosine      ; nan(Number_of_rows_kappa - size(z_coord_tyrosine,1),      size(z_coord_tyrosine,2))      ] ;
        name_coord_tyrosine   = [ name_coord_tyrosine   ; strings(Number_of_rows_kappa - size(name_coord_tyrosine,1), size(name_coord_tyrosine,2)) ] ;
        
        x_coord_phenylalanine    = [ x_coord_phenylalanine    ; nan(Number_of_rows_kappa - size(x_coord_phenylalanine,1),    size(x_coord_phenylalanine,2))    ] ;
        y_coord_phenylalanine    = [ y_coord_phenylalanine    ; nan(Number_of_rows_kappa - size(y_coord_phenylalanine,1),    size(y_coord_phenylalanine,2))    ] ;
        z_coord_phenylalanine    = [ z_coord_phenylalanine    ; nan(Number_of_rows_kappa - size(z_coord_phenylalanine,1),    size(z_coord_phenylalanine,2))    ] ;
        name_coord_phenylalanine = [ name_coord_phenylalanine ; strings(Number_of_rows_kappa - size(name_coord_phenylalanine,1), size(name_coord_phenylalanine,2)) ] ;

        
        X_coord_molecules = horzcat(x_coord_tryptophan,x_coord_tyrosine,x_coord_phenylalanine) ;
        Y_coord_molecules = horzcat(y_coord_tryptophan,y_coord_tyrosine,y_coord_phenylalanine) ;
        Z_coord_molecules = horzcat(z_coord_tryptophan,z_coord_tyrosine,z_coord_phenylalanine) ;
        
        Name_molecules = horzcat(name_coord_tryptophan, name_coord_tyrosine, name_coord_phenylalanine) ;
        
        
        ring_atoms.TRP = ["CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"] ;
        ring_atoms.TYR = ["CG","CD1","CD2","CE1","CE2","CZ"] ;
        ring_atoms.PHE = ["CG","CD1","CD2","CE1","CE2","CZ"] ;
        
        col_TRP = 1 : Number_of_tryptophan ;
        col_TYR = Number_of_tryptophan + (1 : Number_of_tyrosine) ;
        col_PHE = Number_of_tryptophan + Number_of_tyrosine + (1 : Number_of_phenylalanine) ;
        
        
        ring_normal = zeros(Number_of_molecules, 3) ;
        
        colgroups = {col_TRP, 'TRP'; col_TYR, 'TYR'; col_PHE, 'PHE'} ;
        
        
        for g = 1:size(colgroups,1)
            cols = colgroups{g,1} ;
            rtype = colgroups{g,2} ;
            wanted = ring_atoms.(rtype) ;
            
            
            for k = cols
            
                rows = find(ismember(Name_molecules(:,k), wanted)) ;   % <-- name lookup
                
                coords = [ X_coord_molecules(rows,k), ...
                           Y_coord_molecules(rows,k), ...
                           Z_coord_molecules(rows,k) ] ;
                           
                coords = coords(~any(isnan(coords),2), :) ;            % drop NaN rows
                
                c = mean(coords, 1) ;
                [~,~,V] = svd(coords - c, 'econ') ;
                
                % n = V(:,end)' ;
                % n = n / norm(n) ;
                % if n(3) < 0, n = -n ; end
                % ring_normal(k,:) = n ;
                
                n = V(:,end)' ;
                n = n / norm(n) ;
                
                % orientation convention: normal points from ring centroid
                % toward the overall residue centroid (i.e. toward the backbone side).
                
                c_overall = positions_xyz(k,:) ;      % already computed above
                
                v = c_overall - c ;           % c = ring centroid from this loop
                
                if dot(n, v) < 0
                   n = -n ;
                end
                
                ring_normal(k,:) = n ;
                
            end
        end
        
        %% Kappa matrix using centroid-joining vector (positions_xyz)
        kappa_matrix = zeros(Number_of_molecules, Number_of_molecules) ;
        
        for i = 1:Number_of_molecules
            for j = 1:Number_of_molecules
                if i == j
                    kappa_matrix(i,j) = 0 ;   % no self-hop
                    continue
                end
                
                r_ij  = positions_xyz(j,:) - positions_xyz(i,:) ;
                r_hat = r_ij / norm(r_ij) ;          % joining vector from actual centroid
                
                n_i = ring_normal(i,:) ;
                n_j = ring_normal(j,:) ;
                
                cos_ij = dot(n_i, n_j) ;
                cos_i  = dot(n_i, r_hat) ;
                cos_j  = dot(n_j, r_hat) ;
                
                kappa_matrix(i,j) = cos_ij - 3 * cos_i * cos_j ;
            end
        end
        
        Coupling_Constant_joules = zeros(Number_of_molecules,Number_of_molecules) ;
        Fermi_rate_constant_joules = zeros(Number_of_molecules,Number_of_molecules) ;
        
        Distances_matrix_m = Distances_matrix*10.^-10 ;
        
        for First_molecule_CC = 1:Number_of_molecules 
            for Second_molecule_CC = 1:Number_of_molecules 
                Coupling_Constant_joules(First_molecule_CC,Second_molecule_CC) = (constants)*(kappa_matrix(First_molecule_CC,Second_molecule_CC))*(mu_squared_matrix(First_molecule_CC,Second_molecule_CC))./((Distances_matrix_m(First_molecule_CC,Second_molecule_CC)).^3) ;
            end
        end
        
        Coupling_Constant_cm_minus_1 = Coupling_Constant_joules*(5.03*10^22) ; % to convert from Joules to cm-1, multiply by 5.03*10^22.
        hbar_cm_minus_1 = (5.29*10^-12) ; % units are cm-1.s
        c = 3*10^10 ; % units are  cm.s -1.
        Fermi_rate_constant_per_second = (1/(c.*(hbar_cm_minus_1)^2)).*(1./(refractive_index.^4)).*((J_matrix).*Coupling_Constant_cm_minus_1.^2) ;
        % surf(CC_cm_1)
        
        % select only those distances that are within a cutoff distance.
        Cutoff_distance_Angstrom = Inf ;
        Distances_matrix_cutoff = Distances_matrix ; 
        Distances_matrix_cutoff(Distances_matrix_cutoff>Cutoff_distance_Angstrom) = NaN ;
        
        % call the coupling constants that are within this distance.
        Fermi_rate_constant_per_sec_cutoff = Fermi_rate_constant_per_second ;
        Fermi_rate_constant_per_sec_cutoff(isnan(Distances_matrix_cutoff)) = NaN ;
        
        % % now propagate according to a random walk.
        
        % first, make the total prob that anything will happen equal 1.
        Total_hops = sum(sum(~isnan(Fermi_rate_constant_per_sec_cutoff))) ;
        
        Fermi_rate_constant_per_sec_cutoff_zero = Fermi_rate_constant_per_sec_cutoff ;
        Fermi_rate_constant_per_sec_cutoff_zero(isinf(Fermi_rate_constant_per_sec_cutoff_zero)|isnan(Fermi_rate_constant_per_sec_cutoff_zero)) = 0;
        
        % Lifetime parameters - EXACT from original code in nanoseconds
        lifetime_tryptophan_in_tubulin_nanoseconds = 3.6 ; % from Kalra et al., 2023 experiments.
        lifetime_tyrosine_in_tubulin_nanoseconds = 7.5 ; % tyrosine lifetime from citation Feitelson et al., 1964.
        lifetime_phenylalanine_in_tubulin_nanoseconds = 3.39 ; % phenylalanine lifetime from citation Guzow et al., 2004.
        
        % creating a matrix from chromophore lifetimes - EXACT from original
        % lifetime_chromophores_matrix = zeros(Number_of_molecules,Number_of_molecules) ;
        lifetime_chromophores_matrix = zeros(Number_of_molecules,1) ;
        lifetime_chromophores_matrix(trp_locations_on_the_J_and_mu_matrix,1) = lifetime_tryptophan_in_tubulin_nanoseconds ;
        lifetime_chromophores_matrix(tyr_locations_on_the_J_and_mu_matrix,1) = lifetime_tyrosine_in_tubulin_nanoseconds ;
        lifetime_chromophores_matrix(phen_locations_on_the_J_and_mu_matrix,1) = lifetime_phenylalanine_in_tubulin_nanoseconds ;
        
        % Calculate rates
        K_radiative_units_per_ns = 1./lifetime_chromophores_matrix;
        K_ET_units_per_ns = Fermi_rate_constant_per_sec_cutoff_zero*10^-9 ;
        
        
        K_all_per_ns = horzcat(K_ET_units_per_ns, K_radiative_units_per_ns) ;
        

        
        % Monte Carlo simulation - EXACT from original code
        Total_number_of_input_photoexcitations = 10^4 ; % Reduced from 10^5 for speed
        Total_time_ns = 1000 ; % Keep full physical time scale
        
        % setting up the matrices - EXACT from original
        Distances_travelled_as_function_of_time = zeros(Total_number_of_input_photoexcitations,Total_time_ns) ;
        Where_it_moved_all_values = zeros(Total_number_of_input_photoexcitations , Total_time_ns) ;
        
        for photoexcitation_number = 1:Total_number_of_input_photoexcitations 

                % pick an entry of this chormophore containing cylinder - EXACT from original
                
                Starting_molecule(photoexcitation_number) = randsample(Number_of_molecules,1) ;

                Where_it_moved_all_values(photoexcitation_number,1) =  Starting_molecule(photoexcitation_number) ;
                
                Where_it_moved_location = Starting_molecule(photoexcitation_number) ;

            for Time = 1:1:Total_time_ns
                
                acceptor_index = randsample(size(K_all_per_ns,2), 1, true, K_all_per_ns(Where_it_moved_location,:)) ;
                
                if acceptor_index == size(K_all_per_ns,2)
                        
                        break
                        
                else
                
                    % i.e. its non-radiative (or ET), and the hop actually takes place.
                    % Only in non-radiative processes does the 'hop' take place.
                    
                    Where_it_moved_all_values(photoexcitation_number,Time+1) = acceptor_index ;
                    Where_it_moved_location = Where_it_moved_all_values(photoexcitation_number,Time+1) ;
                    Distances_travelled_as_function_of_time(photoexcitation_number,Time+1) = Distances_matrix(Starting_molecule(photoexcitation_number),Where_it_moved_location) ;
                    
                end
            end 
        end

        %identify the last number in the row of values, where the excitation died - EXACT from original
        for all_rows = 1:Total_number_of_input_photoexcitations
                [row,col,Where_it_moved_all_values_last_value(all_rows)]= find(Where_it_moved_all_values(all_rows,1:end),1,'last') ;
        end 

        % find the distance travelled at the end of each excitation - EXACT from original
        Distances_travelled_by_each_excitation = zeros(Total_number_of_input_photoexcitations,1) ;

        for excitation_number_for_diff_length = 1:Total_number_of_input_photoexcitations 
                 Ending_molecule = Where_it_moved_all_values_last_value(excitation_number_for_diff_length) ;
                 Distances_travelled_by_each_excitation(excitation_number_for_diff_length,1) = Distances_matrix(Starting_molecule(excitation_number_for_diff_length),Ending_molecule) ;

                 Diffusion_length_Angstrom_record(excitation_number_for_diff_length) = mean(Distances_travelled_by_each_excitation(1:excitation_number_for_diff_length)) ;
        end

        % making the zeros in the Where_it_moved_all_values array equal to NaN - EXACT from original
        % This will make sure that the zeros are not included in the average.
        Where_it_moved_all_values_remove_zeros = Where_it_moved_all_values ;
        for excitation_number_remove_zeros = 1:1:Total_number_of_input_photoexcitations 
            for time_remove_zeros = 1:1:(Total_time_ns) 
                if time_remove_zeros+1<=Total_time_ns && Where_it_moved_all_values_remove_zeros(excitation_number_remove_zeros,time_remove_zeros) == Where_it_moved_all_values_remove_zeros(excitation_number_remove_zeros,time_remove_zeros+1) 
                   Where_it_moved_all_values_remove_zeros(excitation_number_remove_zeros,time_remove_zeros+1) = NaN ;
                   Where_it_moved_all_values_remove_zeros(excitation_number_remove_zeros,time_remove_zeros+1:end) = NaN ;
                end
            end
        end

        Diffusion_length_Angstrom = mean(Distances_travelled_by_each_excitation)
        
        % Assign to output variable
        diffusion_length = Diffusion_length_Angstrom;
        
    catch ME
        fprintf('Error processing %s: %s\n', pdb_filename, ME.message);
        diffusion_length = NaN;
    end
end
