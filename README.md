The Sequence Calculator 

When developing novel peptide compounds in drug discovery, scientists have to manually calculate the concentration, volume, and mass of each individual amino acid in the sequence. For short peptides this is not a problem, but for very long peptides (100 mer or more) this ia time consuming task which if mistakes are made, can cost scientists weeks. 

The Sequence Calculator is a simple tool that allows scientists to put in their peptide sequence and get out an excel file that gives them the required information to synthesise their peptide of choice on an automated peptide synthesiser.

The scientists need to know how many vials they will need, how many racks they will need, and how much of each amino acid they will need. Further, the scientists are limited by the volume of the vial so they need to know how many times a particular amino acid appears in a sequence before they need a new vial. In this tool, the reactor loop is 2.5 mL and the vials hold 
16 mL of solvent. So the maximum number of times a vial can be sampled is 6 before a new vial is needed (16 ml vial /2.5 ml reactor coil is 6.4 samplings). The calculator first calculates the occurrence of a particular amino acid (the number of times it appears in the peptide) and then calculates the number of vials required. If the amino acid T appears 8 times, the excel sheet will output T, T2 and split the vials so the correct number of samples is taken. You cannot take 8 samplings from a vial that has a 6 sample limit!

Then, in a separate sheet in the same excel sheet, the calculator will output the order the vials are samples for a particular peptide by giving the vials a number. 

The amino_acid.csv file is kept separate from the .exe but in the same directory (folder) so that scientists can add novel un-natural amino acids to the csv file so that the novel sequences can still be run without causing errors. The csv contains the amino acid code (eg, T or W) and their molecular weight (MW) which is needed for calculations.

Due to un-natural amino acids not having a single letter designation, they are put in as Pra for example. This means that the sequence must be put in with spaces between each amino acid letter: T T V Q I Pra P R A instead of TTVQIPraPRA because Pra and P R A are seen as the same thing and so the code has been written to strip spaces between the letters to capture
the three letter un-natural amino acids. (Note, this can be adapted for three letter codons as well).

Inner workings

def resource_path(relative_path) checks that the csv file is located in the same directory as the .exe file

def validate_sequence(sequence, df) uses a panda dataframe validate the amino acids in a sequence and then splits them into tokens

def calculate_sequence_mass(tokens, df) uses tokenised amino acids to calculate mass. The amino acids in the sequence are made into tokens because single letter amino acids and three letter amino acids must have spaces to avoid conflicts that result in inproper interpretation of the peptide sequence.

def find_occurrence(tokens) finds the number of times a single amino acid appears in the peptide sequence so that later, this data can be used to calculate the number of vials required for a particular amino acid

def calculate_aa_mass_volume(tokens, df, conc=0.4, max_per_vial=6, max_volume=16) does the meat of the calculations. First, creates a dictionary and zips the datafram AA and MW together, gets the number of occurrences, and uses max per vial to ensure calculate how to split vials. A vial map using a dictionary is created because the racks are limited to 30 positions. Each amino acid (and split vials) is given a specific numerical position in the rack that the scientist will need to know. This data is output to the excel file for the scientist. Once the for loop hits position 30, the rack number is increased and the new rack is limited to 30 again. This will loop depending on how many vials and racks are required. The for loops determine how many vials are required using the logic that if count (occurrence) is greater than the max_per_vial (set at 6 currently), then the amino acid is split into two vials. The vial map is output as a panda data frame.

def build_synthesis_plan(tokens, vial_map, max_per_vial=6) takes the above funtions outputs to check for related vials and ensure they are named correctly in the output file. If T has an occurrence of 8, this function ensures the vials in the map are named T, T2 for the sequence. This is returned as a dataframe. This function creates the table that tells the machine which order to sample the vials for the synthesis of the input sequence. 

def process_sequence() is all of the tkinter outputs and inputs

All of this codes takes a peptide sequence and gives the scientists the amount of each individual amino acid they need, the occurence of each in a sequence so that a synthesis plan can be devised, and then outputs the vial order required to synthesise the peptide on the automatic peptide synthesiser. 
