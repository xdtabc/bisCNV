# BisCNV - A pipeline which was used to generate and verify the simulated WGBS or WGS data.

BisCNV is a tool for simulating whole genome sequencing (WGS) and whole genome bisulfite sequencing (WGBS) data.

## Features

1. Simulate FASTA sequences

2. Simulate WGS sequencing data (FASTQ format)

3. Simulate WGBS sequencing data (FASTQ format)

## Installation

```shell
git clone https://github.com/your_username/BisCNV.git
cd BisCNV
```



## Dependencies

- Python 3.x
- Biopython
- NumPy
- Pandas
- Samtools

Please ensure that you have installed the above dependencies.

## Usage

BisCNV provides four main functionalities, which can be selected using different parameters:

### 1. Simulate FASTA Sequences

```bash
python BisCNV.py --SIMU_FASTA --input_fa <reference.fa> --input_bed <variants.bed> --output_dir <output_directory>
```

#### Required parameters:

- --SIMU_FASTA: Flag to activate FASTA simulation
- --input_fa: Input reference genome FASTA file
- --input_bed: Input variant BED file

#### Optional parameters:

- --output_dir: Output directory (default: current directory)

### 2. Simulate WGS Sequencing Data

```bash
python BisCNV.py --SIMU_WGS_FASTQ --reference <reference.fa> --seq_length <read_length> --number_of_sequences <num_reads> --output_dir <output_directory>
```

#### Required parameters:

- --SIMU_WGS_FASTQ: Flag to activate WGS simulation
- --reference: Reference genome FASTA file
- --seq_length: Simulated read length

#### Optional parameters:

- --number_of_sequences: Number of simulated reads (default: 1000000)
- --error_rate: Sequencing error rate (default: 0.01)
- --paired_end: Generate paired-end sequencing data
- --output_dir: Output directory (default: current directory)

### 3. Simulate WGBS Sequencing Data

```bash
python BisCNV.py --SIMU_WGBS_FASTQ --reference <reference.fa> --seq_length <read_length> --number_of_sequences <num_reads> --conversion_rate <bisulfite_conversion_rate> --output_dir <output_directory>
```

#### Required parameters:

- --SIMU_WGBS_FASTQ: Flag to activate WGBS simulation
- --reference: Reference genome FASTA file
- --seq_length: Simulated read length

#### Optional parameters:

- --conversion_rate: Uniform bisulfite conversion rate
- --CG_conversion: Conversion rate in CG context
- --CH_conversion: Conversion rate in CH context
- --number_of_sequences: Number of simulated reads (default: 1000000)
- --error_rate: Sequencing error rate (default: 0.01)
- --paired_end: Generate paired-end sequencing data
- --output_dir: Output directory (default: current directory)

Note: Either --conversion_rate or both --CG_conversion and --CH_conversion must be specified.

## Additional Parameters

- --quality: Quality score for simulated reads (default: 40)
- --non_directional: Generate non-directional BS-Seq libraries
- --minfrag: Minimum fragment length for paired-end reads (default: 70)
- --maxfrag: Maximum fragment length for paired-end reads (default: 400)
- --fixed_length_adapter: Length of fixed adapter contamination
- --variable_length_adapter: Length of variable-length adapter contamination

For a complete list of parameters, run:

```bash
python BisCNV.py --help
```



## Output

- Simulated FASTA files (for FASTA simulation)
- Simulated FASTQ files (for WGS and WGBS simulation)

## Notes

- Ensure you have sufficient disk space for output files.
- For large simulations, the process may take some time. Please be patient.

## Examples

Simulate FASTA with variants:

```bash
python BisCNV.py --SIMU_FASTA --input_fa reference.fa --input_bed variants.bed --output_dir ./output
```

Simulate WGS data:

```bash
python BisCNV.py --SIMU_WGS_FASTQ --reference reference.fa --seq_length 150 --number_of_sequences 1000000 --paired_end --output_dir ./output
```

Simulate WGBS data:

```bash
python BisCNV.py --SIMU_WGBS_FASTQ --reference reference.fa --seq_length 150 --conversion_rate 0.98 --number_of_sequences 1000000 --paired_end --output_dir ./output
```



## License

[Add license information]

## Contact

[Add contact information]
