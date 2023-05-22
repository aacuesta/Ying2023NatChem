# Ying2023NatChem

**Pre-requisites**
1. Install blast plus from NCBI.
	1. FTP download: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
	1. NCBI Website: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
1. Install Perl
	1. Download ActiveState Perl installer here: https://www.perl.org/get.html#win32
1. Use the command line to download blast databases
	1. Blast and update_blastdb.pl are usually installed to "C:\Program Files\NCBI\blast-2.10.0+\bin"
	1. From command prompt, run the following commands (you may need to run command prompt as administrator):
		1. cd C:\Program Files\NCBI\blast-2.10.0+\bin\
		1. perl update_blastdb.pl --decompress pdbaa
		1. perl update_blastdb.pl --decompress swissprot
	1. Add a user environment variable called "BLASTDB" and set it to point to the directory where the blast databases are installed
		1. Instructions to edit environment variables can be found here: https://docs.telerik.com/teststudio/features/test-runners/add-path-environment-variables
1. Install the Microsoft Visual C ++ Redistributable package from here, if you do not already have it on your system: https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads
1. A windows binary file is available as listed above. To run from the command line, python3 is required. It is possible to run the python scripts on MacOS after setting the BLASTDB environment variable. 
  1. python pf_lnr_gui.py to start 


**User Instructions**

This program will take an excel file that contains a list of uniprot IDs and associated residue position numbers. With this, it will run a blast search against the pdbaa database, and get a list of PDB files that contain the residue of interest. The resulting PDB IDs will be added to the input file.

Next, it will automatically download each PDB file, and search for whether there is a ligand nearby. It will add the identified protein name, ligand and distance to the residue of interest to the initial input file. It will also determine whether there is an adjancent Lys/Arg and whether the residue of interest is at the interface of two distinct chains in the PDB file.

Note: This will identify the structures of homologous proteins and identify the position in the protein that matches the residue of interest. If this amino acid is nucleophilic, it will find ligands near the nucleophilic atoms. Otherwise, it will look for the nearest atom in the amino acid to the neighboring ligand.  


**Main Window**

![alt text](https://github.com/aacuesta/Ying2023NatChem/blob/main/Main%20Panel.png)
1. The input file must be an .xlsx file.
1. It must also have a column named "Protein" that contains the Uniprot ID, and a column named "Position" that contains the residue number of interest.
1. The output of the program will be an xlsx file.

**Options window**

![alt text](https://github.com/aacuesta/Ying2023NatChem/blob/main/Options%20Panel.PNG)
1. In the options panel, you can specify the parameters you would like to use to run the search.
	1. **Blastp query length** defines the maximum length of the amino acid sequence to blast against the pdb database, pdbaa.
	1. **Max e-value to report**: Blast e-value is usually reported in scientific notation. Set the size of the exponent using this field. Higher values will limit search results to higher quality alignments.
	1. **Max number of PDB files to report**: Will report at most the indicated number of alignments, ordered by best e-value in the resulting Excel file.
	1. **Min Percent Identity** defines the minimum percent identity between the blastp query sequence and resulting PDB file sequence 
		1. ALL results will be used for the ligand finding steps, regardless of the value chosen here. Limit ligand finding to high quality alignments by setting an appropriate e-value.
	1. **Max distance to ligand**: Maximum distance in angstroms between the residue of interest and any atom in a neighboring ligand.
	1. **PDB Directory**: Place to find and store PDB files
		1. This program will attempt to download the PDB files of interest. Be sure to have enough free space on the drive containing the PDB files--they may take **50-100 GB** of space!!

1. You can save the settings you choose as a ".dat" file and load them again later. You can also save you favorite settings as defaults for convenience.

Please contact me with questions and comments:

Adolfo Cuesta

adolfo.cuesta@ucsf.edu

05/21/2023
