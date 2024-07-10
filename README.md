# contact_binding

## Overview

This code pipeline extracts the contact information of protein, RNA, and DNA complexes, comparing it to known binding sites.

## Setup and Installation

1. Install the following dependencies: [conda.3x](https://docs.anaconda.com/miniconda/)
2. Clone the repository.
```
git clone https://github.com/ChangLab/pli-analyzer.git
```
3. Navigate to the directory
```
cd pli-analyzer
```

4. Install dependencies

```
conda env create -f environment.yml
```

5. Activate the environment:
```
conda activate pli-analyzer
```

## Usage

## Input Parameters

| Parameter    | Description                              |
|--------------|------------------------------------------|
| `input_dir` | Path to the input folder with the pdb/cif files.                   |
| `output_dir`  | Path to where the output files should be created       |
| `binding_file_name`    | Name of CSV file containing known binding sequences |
| `distance_threshold`    | Angstrom threshold for residue distances: optional, default 10 |
| `padding`    | # of residues added as "padding" on either side of the residues: optional, default 3 |
| `hide_csv`    | doesn't create intermediate CSVs: default = False |

### Example Usage

```
python contact_binding.py --input_dir /Users/fralian/Desktop/proteins --output_dir /Users/fralian/Desktop/proteins --binding_file_name /Users/fralian/Desktop/proteins/bindings.csv --padding 2 
```

### Output

Intermediate CSV: One is produced for each input pdb/cif file

<img width="1077" alt="Screen Shot 2024-07-10 at 3 46 11 PM" src="https://github.com/frances-liang/contact_binding/assets/114785097/2447612b-0263-42ff-a8c7-81aa05ee7ef0">  


<br>

Final CSV with Binding Ranges:

<img width="670" alt="Screen Shot 2024-07-10 at 3 53 44 PM" src="https://github.com/frances-liang/contact_binding/assets/114785097/5aa74bb9-f345-48f2-82a3-1320226830f6">




Final File





