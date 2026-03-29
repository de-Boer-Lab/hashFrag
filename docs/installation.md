# Installation

It is recommended to execute hashFrag in a `conda` or `virtualenv` environment with Python version 3.10.

## pip install
We recommend installing hashFrag using the following pip command:
```
pip install hashFrag
```

## Clone repository
Clone the repository using the following command:
```
git clone https://github.com/de-Boer-Lab/hashFrag.git
```

Export the source directory to your `PATH`:
```
export PATH="$PATH:./hashFrag/src"
```

(**Optional**) To avoid running the above command every time you open a terminal, add it to your shell configuration file (e.g., `~/.bashrc`) with the following command:
```
echo 'export PATH="$PATH:./hashFrag/src"' >> ~/.bashrc
```

## Installing dependencies

If you are managing your virtual environment with Anaconda or Miniconda, you can directly install dependencies upon creation of the conda environment using the command:
```
conda env create -n hashFrag -f environment.yml
```
* This creates a `conda` environment named "hashFrag"

Alternatively, you can install dependencies located in the `requirements.txt` file with the folowing `pip` command:
```
pip install -r requirements.txt
```

## BLAST+ download

To install the suite of BLAST applications, follow the instructions found at the NCBI BLAST Command Line Applications User Manual [here](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

Directly download the BLAST+ package executables for different operating systems at the following FTP page:
> https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Follow the instructions to extract the downloaded file and export binaries to your `PATH`. 

Verify that `blastn` and `makeblastdb` commands were installed succesfully:
```
blastn -version
makeblastdb -version
```

### A minor note on the use of BLAST alignment scores in hashFrag

The `lightning` (default) mode of hashFrag uses pairwise local alignment scores derived from the BLAST algorithm. Rather than directly using the provided alignment scores, however, a *corrected* version of the alignment score is calculated.

Gap scoring is designed to reflect the biological occurrence of insertions and deletions in sequences. Typically, opening a gap incurs a larger penalty (`gapopen`), while subsequent extension of the same gap incurs smaller penalties (`gapextend`). Upon encountering a gap *opening* event, the BLAST algorithm applies the `gapextend` penalty in addition to a `gapopen` penalty. To conform to exact local alignment scoring conventions, we adjust the BLAST scores such that only the `gapopen` penalty is applied to gap opening events.  