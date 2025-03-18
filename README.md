# TopLib 

TopLib is a Python package for building and searching spectral libraries for
top-down mass spectrometry (MS).  

### Requirements
* Python version 3.8 or higher

If you use Toplib in your work, please cite the following publication:

* Kun Li, Haixu Tang, and Xiaowen Liu. TopLib: Building and searching top-down mass spectral libraries for proteoform identification (2024) bioRxiv preprint.
 

## 1. Building a spectral library using one top-down MS data file 

### 1.1 Top-down MS data preprocessing 

MsConvert and TopPIC Suite are needed for preprocessing top-down MS data files.

### Software tools  

* MsConvert: https://proteowizard.sourceforge.io/download.html
* TopPIC Suite version 1.7.8 or higher: https://www.toppic.org/software/toppic/index.html

In data preprocessing, raw MS files are converted into centroided mzML files using msConvert, 
mzML files are deconvoluted using TopFD, deconvoluted mass spectra are searched against 
a protein sequence database for spectral identification using TopPIC, and
finally proteoform-spectrum-matches (PrSMs) reported by TopPIC are filtered using a Python script. Suppose
the raw MS file is lib_spectra.raw and the protein sequence database is
proteins.fasta. Below is an example for data preprocessing.    

* Use msConvert to convert the raw MS data file lib_spectra.raw to a centroided mzML file lib_spectra.mzML 
* Use TopFD (https://www.toppic.org/software/toppic/) to deconvolute the tandem mass spectrometry (MS/MS) spectra in the mzML file lib_spectra.mzML to generate an msalign file ```lib_spectra_ms2.msalign```
* Use TopPIC (https://www.toppic.org/software/toppic/) to search the msalign file lib_spectra_ms2.msalign against the protein sequence database proteins.fasta for spectral identification. The resulting files ```lib_spectra_ms2_toppic_prsm_single.tsv``` and ```lib_spectra_ms2_toppic_proteoform_single.tsv``` are used in the next step.  
* Filter out inconsistent PrSMs reported by TopPIC 

  Run the command: 
  ```
  python3 tsv_file_processing.py lib_spectra_ms2_toppic_prsm_single.tsv lib_spectra_ms2_toppic_proteoform_single.tsv
  ```
The resulting file is ```lib_spectra_ms2_toppic_prsm_single_filtered.tsv```.
#### Standalone executable:
Command to filter out inconsistent PrSMs without requiring a Python environment:

```
tsv_file_processing.exe lib_spectra_ms2_toppic_prsm_single.tsv lib_spectra_ms2_toppic_proteoform_single.tsv
```

### 1.2 Building a top-down mass spectral library  
The Python script ms_library_building.py builds a top-down spectral library using a preprocessed top-down MS data set. 

* Input files: 
  * An msalign file: ```lib_spectra_ms2.msalign```
  * A tsv file containing filtered PrSM identifications: ```lib_spectra_ms2_toppic_prsm_single_filtered.tsv```

* Argument: 
  * The type of representative spectra: ```average``` or ```single```
  * The PPIR cutoff value (0-1): ```-c <0.5>```. Default value: 0
 
* Output: 
  * A top-down spectral library built using the MS/MS spectra and identifications in the input files, which is stored in a sqlite file lib_spectra_ms2.db     

Run the command to generate a library with average representative spectra: 

```
python3 ms_library_building.py lib_spectra_ms2.msalign lib_spectra_ms2_toppic_prsm_single_filtered.tsv average
```
#### Standalone executable:
Run ms_library_building.exe to create the spectral library without requiring a Python environment. 

```
ms_library_building.exe lib_spectra_ms2.msalign lib_spectra_ms2_toppic_prsm_single_filtered.tsv average
```
### 1.3 Convert a SQL-based library to a text-based spectrum library
The Python script db_msalign_tsv.py converts a SQL-based library (.db) to a text-based spectrum library. A text-based spectrum library is stored in two text files: one is a msalign file contained the meta information and the precursor and fragment masses of the spectra, the other is a tab-delimited text file containing the information of proteoform identifications of the spectra. 

* Input files: 
  * A db file: ```lib_spectra_ms2.db```
    
* Output: 
  * An msalign file: ```lib_spectra_m2_representative.msalign```
  * A tsv file: ```lib_spectra_m2_identification.tsv```       

Run the command to convert a SQL-based library to a text-based library:

```
python3 db_to_msalign_tsv.py lib_spectra_ms2.db
```

## 2. Top-down mass spectral identification by library search 

In top-down MS data preprocessing, suppose the raw query MS file is query_spectra.raw, use msConvert and TopFD to generate an msalign file ```query_spectra_ms2.msalign``` containing deconvoluted MS/MS spectra (see 1.1).      

The Python script ms_library_query.py searches top-down MS/MS spectra against a top-down spectral library for spectral identification.  

* Input files: 
  * A top-down mass spectral library: ```lib_spectra_ms2.db```
  * An msalign file containing deconvoluted MS/MS spectra: ```query_spectra_ms2.msalign```

* Arguments:
  * Precursor mass error tolerance type: ```-T <type>``` (type is ppm or Da). Default value: ppm
  * Precursor mass error tolerance: ```-E <a positive number>```. Default value: 10 ppm or 2.2 Da 
  * Fragment mass error tolerance (in ppm): ```-e <a positive number>```. Default value: 10 ppm  
  * Charge matching is required: ```-c <yes or no>```. Default value: no (not required). 

* Output: 
  * An TSV file containing the query results: ```toplib_output/query_res.tsv```

Run the command to search query_spectra_ms2.msalign against lib_spectra_ms2.db
with an error tolerance of 10 ppm for precursor and fragment masses. Charge matching is not required. 
```
python3 ms_library_query.py lib_spectra_ms2.db query_spectra_ms2.msalign 
```

## 3. Building a spectral library using multiple top-down MS data files
### 3.1. Creating a database for storing mass spectra

#### 3.1.1. Creating a sqlite database  

The command below creates a sqlite database toplib.db.

```
python3 db_gen.py 
```

#### 3.1.2. Adding a project to the database 
```
python3 db_add_project.py 
```

The script will ask the user to input the project name and project description.
An example is given below. 

  * Project name: Top-down MS study of human colorectal cells 
  * Project description: Top-down MS study of human colorectal metastatic (SW620) and nonmetastatic (SW480) cells


#### 3.1.3. Adding a sample associated with the project
```
python3 db_add_sample.py 
```
The script will ask the user to input the sample information. An example is
given below. 

  * Species name (choose one of the options by entering the corresponding number): ```1 = human; 2= mouse.```
  * Sample name: ```SW480 cells```
  * Sample description: ```SW480 cells```
  * Project id: ```1```

#### 3.1.4. Adding an MS method
```
python3 db_add_method.py   
```
The script will ask the user to input the method information. An example is
given below. 

  * Instrument (choose one of the options by entering the corresponding number): ```1 = Thermo Q Exactive HF; 2 = Thermo Orbitrap Fusion Lumos; 3 = Thermo Orbitrap Eclipse.```
  * Dissociation method (choose one of the options by entering the corresponding number): ```1 = HCD; 2 = CID; 3 = ECD; 4 = ETD.```
  * Collision energy: ```20%```
  * Resolution: ```120000```

#### 3.1.5. Adding an experiment associated with the project
```
python3 db_add_experiment.py 
```
The script will ask the user to input the experiment information. An example is
given below. 

  * Experiment name: ```SW480 2D replicate 1```
  * Experiment description: ```RPLC-CZE 2D seperation is used to analyze
    proteoforms in SW480 cells. The first replicate in a triplicate experiment.```
  * Sample id: ```1```
  * Method id: ```1```

The user can add projects, samples, methods, and experiments as they needed.

### 3.2 Adding MS data files
Use the methods in 1.1 to preprocess top-down MS data files. 

Use Python script db_add_file.py Add preprocessed MS data files one by one to the mass spectral database.  
Suppose the experiment ID is 1 for the preprocessed MS data files obtained in 1.1. The command below is used to add the mass spectra to the database. Average and single representative spectra are also generated and stored in the database.    

```
python3 db_add_file.py lib_spectra_ms2.msalign lib_spectra_ms2_toppic_prsm_single_filtered.tsv 1
```

### 3.3 Extracting representative spectra 
After all MS data files are added to the database, all representative spectra of the files can be extracted from the database using Python script db_query.py. 

* Input parameter:
  * Project ID: e.g., ```1```
  * Type of representative spectra: ```average``` or ```single```

* Output:
  * An msalign file containing all representaive spectra, which is stored as: ```lib_spectra_ms2.msalign```
  * A TSV file containing proteoform identifications of the representative spectra, which is stored as: ```lib_spectra_ms2_toppic_prsm_single_filtered.tsv```
    
Run the command to extract all single representative spectra of project 1:
```
python3 db_query.py 1 single
```
After running this command, an msalign file and a TSV file will be generated.

### 3.4 Building a spectral library for library search  
Follow the method in Section 1.2 to build a spectral library using the msalign
file and TSV file reported in Section 3.3. 

### 4. Standalone executable
You can also run TopLib as a standalone executable without needing to install a Python environment or additional packages. All required dependencies are packaged into .exe files. Run the command lines:

```
tsv_file_processing.exe lib_spectra_ms2_toppic_prsm_single.tsv lib_spectra_ms2_toppic_proteoform_single.tsv
ms_library_building.exe lib_spectra_ms2.msalign lib_spectra_ms2_toppic_prsm_single_filtered.tsv average
ms_library_query.exe lib_spectra_ms2.db query_spectra_ms2.msalign 
toplib_comprehensive_building.exe
db_add_file.exe lib_spectra_ms2.msalign lib_spectra_ms2_toppic_prsm_single_filtered.tsv 1
db_query.exe 1 single
```

## 5. Contact
If you have any questions, please send an email to Dr. Kun Li: kil7@tulane.edu.
