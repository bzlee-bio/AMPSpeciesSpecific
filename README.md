# Species-specific prediction of antimicrobial activity
The tool predicts antimicrobial activity of peptides for five bacterial species, <br><i>Bacillus subtilis, Escherichia coli, Pseudomonas aeruginosa, Staphylococcus aureus, Staphylococcus epidermidis</i>


## Dependencies
To install python dependencies, run: `pip install -r requirements.txt`

## Input file 
Input file type is fasta format in which amino acids are represented using single-letter codes.

Detailed information of fasta links: https://en.wikipedia.org/wiki/FASTA_format

## Output file
Output file contains peptide list with prediction probabilities as AMPs for five bacterial species.

## Running

`python AMP_species_prediction.py --fasta <input_fasta_file.fasta> --out <output_file_name> --windowsize <window size> --step <step size>` <br><br>
A window size and step size is optional.<br>
If the window size is provoded, input sequences will be truncated into the size of the window with step size.
<br><br>




## Citation
Lee B, Shin MK, Yoo JS, Jang W and Sung J-S (2022) Identifying novel antimicrobial peptides from venom gland of spider Pardosa astrigera by deep multi-task learning. Front. Microbiol. 13:971503. doi: 10.3389/fmicb.2022.971503
http://doi.org/10.3389/fmicb.2022.971503

