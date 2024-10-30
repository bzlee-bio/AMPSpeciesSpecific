# Species-specific prediction of antimicrobial activity
The tool predicts antimicrobial activity of peptides for five bacterial species, <br><i>Bacillus subtilis, Escherichia coli, Pseudomonas aeruginosa, Staphylococcus aureus, Staphylococcus epidermidis</i>


## Dependencies
To install python dependencies, run: `pip install -r requirements.txt`

## Input file 
Input file type is fasta format in which amino acids are represented using single-letter codes.

Detailed information of fasta links: https://en.wikipedia.org/wiki/FASTA_format

## Output file
Output file contains peptide list with prediction probabilities as AMPs for five bacterial species.

## Usage
To run prediction tool, use the following command:
`python AMP_species_prediction.py --fasta <input_fasta_file.fasta> --out <output_file_name> --windowsize <window size> --step <step size>` <br><be>
#### Parameters:
-	--fasta: (Required) Path to the input FASTA file containing peptide sequences.
-	--out: (Required) Desired name for the output file.
-	--windowsize: (Optional) Size of the window to truncate input sequences.
-	--step: (Optional) Step size for truncating sequences.
### Notes:
- If the --windowsize is provided, input sequences will be truncated into segments of the specified window size with the defined step size.
### Example:
`python AMP_species_prediction.py --fasta peptides.fasta --out predictions.csv --windowsize 20 --step 5`
<br><br>

## Prediction results of <i>Pardosa astrigera</i>
File: PA_pred_results.csv
Leveraging transcriptomic data from <i>Pardosa astrigera</i>, a deep learning model was employed to predict antimicrobial peptides. The results demonstrate high prediction probabilities for peptides effective against the five specified bacterial species.


## Citation
If you use this tool or prediction results of <i>Pardosa astrigera</i> in your research, please cite the following paper:

### Prediction model
Lee B, Shin MK, Yoo JS, Jang W and Sung J-S (2022) Identifying novel antimicrobial peptides from venom gland of spider <i>Pardosa astrigera</i> by deep multi-task learning. Front. Microbiol. 13:971503. doi: 10.3389/fmicb.2022.971503
http://doi.org/10.3389/fmicb.2022.971503

### Prediction results for <i>Pardosa astrigera</i>
To be added.

