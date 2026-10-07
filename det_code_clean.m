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
        
        
        positions_xyz = vertcat(Centroid_tryptophan,Centroid_tyrosine,Centroid_phenylalanine) ;
        
        % Calculate distance matrix - EXACT from original ; PDB records coordinates in Angstroms only
        Distances_matrix_Angstrom = zeros(Number_of_molecules,Number_of_molecules) ;
        
        
        for First_molecule = 1:Number_of_molecules 
            for Second_molecule = 1:Number_of_molecules 
                Distances_matrix_Angstrom(First_molecule,Second_molecule) = pdist2(positions_xyz(First_molecule,:),positions_xyz(Second_molecule,:)) ; 
                % magn_xy = norm(positions_xyz(First_molecule,1:2))*norm(positions_xyz(Second_molecule,1:2)) ;
                % z_12 = positions_xyz(Second_molecule,3) - positions_xyz(First_molecule,3) ;
                % theta = acos(dot(positions_xyz(First_molecule,1:2),positions_xyz(Second_molecule,1:2))/(magn_xy));
                % if abs(theta)<1e-6
                %    theta=0;
                % end
                % x_arc = Microtubule_radius_Angstrom*theta ;
                % Distances_matrix_Angstrom(First_molecule,Second_molecule) = real(sqrt(x_arc^2 + z_12^2)) ;
            end
        end
        
        
        Distances_matrix_m = Distances_matrix_Angstrom*10.^-10 ;
        
        % from here on, we calculate J values.
        % these are the J-values of each chromophore pair.
        J_TRP_TRP = 0.00058 ; % units are cm, from your calculations.
        J_TYR_TYR = 0.0012 ; % units are cm, from your calculations.
        J_PHE_PHE = 0.00098 ; % from your calculations.
        
        J_TRP_TYR = 0.000487 ; % units are cm, from your calculations.
        J_TYR_TRP = 0.0011 ; % units are cm, from your calculations.
        J_PHE_TYR = 0.0012 ; % units are cm, from your calculations.
        J_PHE_TRP = 0.0039 ;
        J_TRP_PHE = 0.00031 ;
        J_TYR_PHE = 0.00098 ;
        
        
        J_matrix = zeros(Number_of_molecules,Number_of_molecules) ;
        
        trp_locations_on_the_J_matrix = 1:Number_of_tryptophan ; 
        tyr_locations_on_the_J_matrix = (Number_of_tryptophan+1):(Number_of_tyrosine+Number_of_tryptophan) ; 
        phen_locations_on_the_J_matrix = (Number_of_tyrosine+Number_of_tryptophan+1):Number_of_molecules ;
        % this is to generate the matrix containing the J value relationships of each chromophore.
        
        % this is for energy transfer in the forward direction i.e from higher
        % wavelengths (of donor chromophores) to lower wavelengths (of accepter
        % chromophores), but those of the same chromophore. So trp-to-trp for
        % example.
        J_matrix(trp_locations_on_the_J_matrix,trp_locations_on_the_J_matrix) = J_TRP_TRP ;
        J_matrix(tyr_locations_on_the_J_matrix,tyr_locations_on_the_J_matrix) = J_TYR_TYR ;
        J_matrix(phen_locations_on_the_J_matrix,phen_locations_on_the_J_matrix) = J_PHE_PHE ;
        
        % this is to generate the matrix containing the J value relationships of each chromophore.
        % this is for energy transfer in the forward direction i.e from higher
        % wavelengths (of donor chromophores) to lower wavelengths (of accepter
        % chromophores). So trp-to-tyr for example.
        
        % tryptophan to tyrosine.
        J_matrix(trp_locations_on_the_J_matrix,tyr_locations_on_the_J_matrix) = J_TRP_TYR ;
        
        %tyrosine to tryptophan
        J_matrix(tyr_locations_on_the_J_matrix,trp_locations_on_the_J_matrix) = J_TYR_TRP ;
        
        % tryptophan to phenylalanine.
        J_matrix(trp_locations_on_the_J_matrix,phen_locations_on_the_J_matrix) = J_TRP_PHE ;
        
        % phen to tryptophan.
        J_matrix(phen_locations_on_the_J_matrix,trp_locations_on_the_J_matrix) = J_PHE_TRP ;
        
        % tyrosine to phenylalanine.
        J_matrix(tyr_locations_on_the_J_matrix,phen_locations_on_the_J_matrix) = J_TYR_PHE ;
        
        % phen to tyrosine.
        J_matrix(phen_locations_on_the_J_matrix,tyr_locations_on_the_J_matrix) = J_PHE_TYR ;
        
        % Now we will determine coupling constants here on out.
        
        % Given that this is the distance and coupling constant relationship, find the
        % coupling values for all other distances.
        
        
        DET_coupling_constant_extrapolated_eV = ((3.22e-2)/( 1.4 * 1.4)) * exp(-(1e10*2.5579)*(Distances_matrix_m - 4.15e-10));
        
        DET_coupling_constant_extrapolated_cm_minus_1 = DET_coupling_constant_extrapolated_eV*(1.602e-19)*(5.03*10^22) ; % to convert from eV to cm-1, multiply by elementary charge 1.602e-19, multiply by 5.03*10^22. spectroscopic wavelength number ;
        
        
        % make all the diagonal elements zero in the DET coupling constant matrix.
        DET_coupling_constant_extrapolated_cm_minus_1 = DET_coupling_constant_extrapolated_cm_minus_1 - diag(diag(DET_coupling_constant_extrapolated_cm_minus_1));
        
        hbar_cm_minus_1 = (5.29*10^-12) ; % units are cm-1.s
        c = 3*10^10 ; % units are  cm.s -1.
        DET_Fermi_rate_constant_extrapolated_per_second = (1/(c.*(hbar_cm_minus_1)^2)).*((J_matrix).*(DET_coupling_constant_extrapolated_cm_minus_1.^2)) ;
        
        % select only those distances that are within a cutoff distance.
        Cutoff_distance_m = Inf ;
        Distances_matrix_cutoff_m = Distances_matrix_m ; 
        Distances_matrix_cutoff_m(Distances_matrix_cutoff_m>Cutoff_distance_m) = NaN ;
        
        % call the rate constants for these -cutoff- distance.
        DET_Fermi_rate_constant_per_sec_cutoff = DET_Fermi_rate_constant_extrapolated_per_second ;
        DET_Fermi_rate_constant_per_sec_cutoff(isnan(Distances_matrix_cutoff_m)) = NaN ;
        
        % these three lines make the diagonal elements zero in the rate matrix.
        % DET_Fermi_rate_constant_per_sec_cutoff_size = size(DET_Fermi_rate_constant_per_sec_cutoff);
        % index_rate = 1:DET_Fermi_rate_constant_per_sec_cutoff_size(1)+1:DET_Fermi_rate_constant_per_sec_cutoff_size(1)*DET_Fermi_rate_constant_per_sec_cutoff_size(2);  % Indices of the main diagonal
        % DET_Fermi_rate_constant_per_sec_cutoff(index_rate) = DET_Fermi_rate_constant_per_sec_cutoff(index_rate) * 0 ; 
        
        % % now propagate according to a random walk.
        
        % first, make the total prob that anything will happen equal 1.
        Total_hops = sum(sum(~isnan(DET_Fermi_rate_constant_per_sec_cutoff))) ;
        
        
        % converted to per microsecond
        DET_rate_constant_per_microsecond_cutoff_zero = (10^-6).*DET_Fermi_rate_constant_per_sec_cutoff ;
        DET_rate_constant_per_microsecond_cutoff_zero(isinf(DET_rate_constant_per_microsecond_cutoff_zero)|isnan(DET_rate_constant_per_microsecond_cutoff_zero)) = 0;
        
        % lifetime of tryptophan, tyrosine and phenylalanine. in microsecond
        lifetime_tryptophan_in_tubulin_microseconds = 7.3 ; % from Bent and Hayons 1975 paper. Table 1 in the trp paper.
        lifetime_tyrosine_in_tubulin_microseconds = 5.6 ; % from Bent and Hayons 1975 paper. Table 1 in the tyr paper.
        lifetime_phenylalanine_in_tubulin_microseconds = 3.1 ; % from Bent and Hayons 1975 paper. Table 1 in the phen paper.
        
        % creating a matrix from chromophore lifetimes.
        lifetime_chromophores_matrix = zeros(Number_of_molecules,1) ;
        lifetime_chromophores_matrix(trp_locations_on_the_J_matrix,1) = lifetime_tryptophan_in_tubulin_microseconds ;
        lifetime_chromophores_matrix(tyr_locations_on_the_J_matrix,1) = lifetime_tyrosine_in_tubulin_microseconds ;
        lifetime_chromophores_matrix(phen_locations_on_the_J_matrix,1) = lifetime_phenylalanine_in_tubulin_microseconds ;
        
        
        
        % these are all the number of events in a second. This is where you
        % interconvert time units from ns to s or whatever ; NO NEED
        
        % these are all the number of events in a nanosecond.
        K_radiative_units_per_microsecond  = 1./lifetime_chromophores_matrix ; % units are per nanosecond. 
        
        
        K_DET_units_per_microsecond = DET_rate_constant_per_microsecond_cutoff_zero ; % units are per microsecond. 
        
        
        K_all_per_microsecond = horzcat(K_DET_units_per_microsecond, K_radiative_units_per_microsecond) ;
        
        
        % the number of times you want to input a photon onto the system.
        Total_number_of_input_photoexcitations = 10^4 ;
        
        % this is the total time duration of the random walk units are nanosecond ;
        Total_time_microsecond = 1000 ;
        
        % setting up the matrices.
        Distances_travelled_as_function_of_time = zeros(Total_number_of_input_photoexcitations,Total_time_microsecond) ;
        Where_it_moved_all_values = zeros(Total_number_of_input_photoexcitations , Total_time_microsecond) ;
        
        
        for photoexcitation_number = 1:Total_number_of_input_photoexcitations 
        
            % pick an entry from all available chromophores
            Starting_molecule(photoexcitation_number) = randsample(Number_of_molecules,1) ;
            
            Where_it_moved_all_values(photoexcitation_number,1) =  Starting_molecule(photoexcitation_number) ;
            
            Where_it_moved_location = Starting_molecule(photoexcitation_number) ;
            
            
            for Time = 1:1:Total_time_microsecond
                
                acceptor_index = randsample(size(K_all_per_microsecond,2), 1, true, K_all_per_microsecond(Where_it_moved_location,:)) ;
                
                if acceptor_index == size(K_all_per_microsecond,2)
                        
                        break
                        
                else
                
                    % i.e. its non-radiative (or ET), and the hop actually takes place.
                    % Only in non-radiative processes does the 'hop' take place.
                    
                    Where_it_moved_all_values(photoexcitation_number,Time+1) = acceptor_index ;
                    Where_it_moved_location = Where_it_moved_all_values(photoexcitation_number,Time+1) ;
                    Distances_travelled_as_function_of_time(photoexcitation_number,Time+1) = Distances_matrix_Angstrom(Starting_molecule(photoexcitation_number),Where_it_moved_location) ;
                    
                end
            end
        end
        
        
        %identify the last number in the row of values, where the excitation died.
        for all_rows = 1:Total_number_of_input_photoexcitations
            [row,col,Where_it_moved_all_values_last_value(all_rows)]= find(Where_it_moved_all_values(all_rows,1:end),1,'last') 
        end 
        
        % find the distance travelled at the end of each excitation.
        Distances_travelled_by_each_excitation = zeros(Total_number_of_input_photoexcitations,1) ;
        
        for excitation_number_for_diff_length = 1:Total_number_of_input_photoexcitations 
            Ending_molecule = Where_it_moved_all_values_last_value(excitation_number_for_diff_length) ;
            Distances_travelled_by_each_excitation(excitation_number_for_diff_length,1) = Distances_matrix_Angstrom(Starting_molecule(excitation_number_for_diff_length),Ending_molecule) ;
        end
        
        Diffusion_length_Angstrom = mean(Distances_travelled_by_each_excitation)
        
        diffusion_length = Diffusion_length_Angstrom;
        
    catch ME
        fprintf('Error processing %s: %s\n', pdb_filename, ME.message);
        diffusion_length = NaN;
    end
end

