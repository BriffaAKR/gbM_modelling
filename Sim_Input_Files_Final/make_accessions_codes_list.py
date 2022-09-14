# load in some standard python libraries

import numpy as np
import pandas as pd
import os
import csv


input_file = "m1001_Col0like_samples_non-redundant_mCG_density_coverage_50.txt"
Col0like_accessions_df = pd.read_csv(input_file, sep="\t{1}",engine='python')
Col0like_accessions_list = Col0like_accessions_df['SRA_Accession'].tolist()
print(len(Col0like_accessions_list))

# write out numbered list of accessions (two columns) 

output_file_name = 'Accession_codes_Col0like_NonRed_Cov50.txt'
open(output_file_name,'w').close

temp_output_file = output_file_name
with open(temp_output_file, 'a', newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')

    temp_output_list = ['Accession_number','Accession_code']
    line_writer.writerow(temp_output_list)
    for i_ in range(len(Col0like_accessions_list)):
        temp_output_list = [i_+1,Col0like_accessions_list[i_]]
        line_writer.writerow(temp_output_list)

