% All codes were created by using MATLAB version 2019a Update 9 (9.6.0.1472908) 64-bit (maci64), 
% 2023b Update 1 (23.2.0.2380103) 64-bit (maca64), or 2024b Update 1 (24.2.0.2740171) 64-bit (maca64).
% All of them have been verified to work correctly with 2024b Update 1 (24.2.0.2740171) 64-bit (maca64).


% "Statistics and Machine Learning Toolbox" is required.


% Description of the variables in the attached MAT-file:

    % Raman data --- Source: this paper
    % - 'rawraman_14strains'
    %     14-by-2 cell array; 1st column: strain name, 2nd column: Raman spectra.
    %     In each cell of the second column, each row represents a single-cell Raman spectrum. 
    %     All the 14 strains are included. 
    % - 'rawraman_9strains'
    %     9-by-2 cell array; 1st column: strain name, 2nd column: Raman spectra. 
    %     In each cell of the second column, each row represents a single-cell Raman spectrum. 
    %     Out of the 14 strains in total, nine with transcriptome data are included.
    % - 'rawramanindex2Ramanshift'
    %     2-by-2 cell array; 1st row: indices of spectra in 'rawraman_14strains' and 'rawraman_9strains', 2nd row: corresponding Raman shifts (1/cm).
    %     (within the fingerprint region)

    % transcriptome data --- Source: NCBI GEO accession numbers 
    %                                   GSE43643 (Y. Katayama, et al., Antimicrobial Agents and Chemotherapy 60, 3730 (2016).), 
    %                                   GSE42218 (Y. Aiba et al., Antimicrobial Agents and Chemotherapy 57, 4861 (2013).), 
    %                                   GSE50826 (M. Saito, et al., Antimicrobial agents and chemotherapy 58, 5024 (2014).)
    % The following variables have the corresponding order of rows (9 strains) and columns (2628 genes).
    % - 'genesmean'
    %     9-by-2628 array; gene expression abundances.
    % - 'genesse'
    %     9-by-2628 array; gene expression abundance standard errors.
    % - 'genes_strainnames'
    %     9-by-1 array; strain names.
    % - 'genes_description'
    %     3-by-2628 array; 1st row: ORF ID, 2nd row: gene name, 3rd row: product.

    % vancomycin MIC data --- Source: this paper
    % The following variables have the corresponding order of rows (14 strains).
    % - 'MIC_VCM_14strains48h'
    %     14-by-1 array; MIC in μg/mL.
    % - 'MIC_VCM_14strains48h_strainnames'
    %     14-by-1 array; strain names.

    % The Raman data, the transcriptome data, and the MIC data have the same order of the strains.


% Codes:

% [Analysis of all the 14 strains]
    % Apply PC-LDA to Raman spectra
    % (Necessary to perform all analyses involving Raman data)
    raman_pclda_14strains

        % Plot PC-LDA Raman (Fig. 1D, 1E, and 1F)
        plot_ldaraman
        visualize_ldaraman_tsne

        % Correlation between Raman LDA1 and vancomycin MIC (Fig. 1G)
        scatterplot_lda1_mic

        % Vancomycin MIC prediction from Raman data and model selection (Fig. 1H)
        mic_prediction_modelselection

% [Analysis of 9 strains with transcriptome data out of all the 14 strains]
    % Apply PC-LDA to Raman spectra
    % (Necessary to perform all analyses involving Raman data)
    raman_pclda_9strains

        % Plot PC-LDA Raman (Fig. 2A and 2B)
        plot_ldaraman

        % Predict Raman data from transcriptome data and perform permutation test (Fig. 2C and 2D)
        raman_prediction_loosv

        % Predict proteome data from Raman data and perform permutation test (Fig. 3A, 3B and 3C)
        omics_prediction_loocv

        % Low-dimensional map of transcriptome profiles on LDA Raman space (Fig. 4)
        normalized_coefficient_space

        % heterogeneity of single-cell Raman spectra (Fig. 5A and 5B)
        spectral_heterogeneity

