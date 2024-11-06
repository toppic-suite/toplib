# TopLib 

TopLib is a Python package for building and searching spectral libraries for
top-down mass spectrometry (MS).  

### Requirements
* Python version 3.8 or higher

If you use Toplib in your work, please cite the following publication:

* Kun Li, Haixu Tang, and Xiaowen Liu. TopLib: Building and searching top-down mass spectral libraries for proteoform identification (2024) bioRxiv preprint.
 

## 1. Building a spectral library using a top-down MS data set 

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


### 1.2 Building a top-down mass spectral library  
The Python script ms_library_building.py builds a top-down spectral library using a preprocessed top-down MS data set. 

* Input files: 
  * An msalign file: ```lib_spectra_ms2.msalign```
  * A tsv file containing filtered PrSM identifications: ```lib_spectra_ms2_toppic_prsm_single_filtered.tsv```

* Argument: 
  * The type of representative spectra: ```average``` or ```single``` 

 
* Output: 
  * A top-down spectral library built using the MS/MS spectra and identifications in the input files, which is stored in a sqlite file lib_spectra_ms2.db     

Run the command to generate a library with average representative spectra: 
```
python3 ms_library_building.py lib_spectra_ms2.msalign lib_spectra_ms2_toppic_prsm_single_filtered.tsv average
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
  * Charge matching is required: ```-c```. Default value: not required. 

* Output: 
  * An TSV file containing the query results: ```toplib_output/query_res.tsv```

Run the command to search query_spectra_ms2.msalign against lib_spectra_ms2.db
with an error tolerance of 10 ppm for precursor and fragment masses. Charge matching is not required. 
```
python3 ms_library_query.py lib_spectra_ms2.db query_spectra_ms2.msalign 
```

## 3. Building a spectral library using multiple top-down MS data sets
### 3.1 TopLib library creation
This creates TopLib database tables.

Run the command: 
```
python3 db_gen.py 
```

### 3.2 Adding MS Data sets 
Run the following command to add datasets to TopLib based on your experiment setup. This allows you to create a comprehensive library entries for your project.

* Input parameter:
  * project name: e.g., ```TDP study of human CRC cell lines```
  * project description: e.g., ```metastatic (SW620) and nonmetastatic (SW480)```

Run the command:   
```
python3 db_add_project.py 
```

And add sample information associated with the project.

* Input parameter:
  * species name (choose one of the options by entering the corresponding number): ```1 = human; 2= mouse.```
  * sample name: e.g.,```SW480 2d replicate 1```
  * sample description: ```SW480 cell```
  * project id: e.g., ```1```

Run the command: 
```
python3 db_add_sample.py 
```

Add method information for your experiment setting. 

* Input parameter: 
  * instrument (choose one of the options by entering the corresponding number): ```1 = Q Exactive HF; 2 = Orbitrap Fusion Lumos; 3 = LTQ Orbitrap Elite.```
  * dissociation method (choose one of the options by entering the corresponding number): ```1 = HCD; 2 = CID; 3 = ECD; 4 = EID; 5 = ETD.```
  * energy value (normalized collision energy): e.g.,```0.2```
  * resolution: e.g.,```120000```

Run the command:
```
python3 db_add_method.py   
```

Add file information.

* Input parameter:
  * msalign file (with extension): e.g., ```sw480_2d_combined_ms2.msalign```
  * identification file (with extension): e.g., ```sw480_2d_combined_ms2_toppic_prsm_single_filtered.tsv```
  * sample ID: e.g., ```1```
  * method ID: e.g., ```1```

Run the command:
```
python3 db_add_file.py  
```

After adding data to TopLib, spectral representatives table will be generated by inputing parameters:
  * file ID: e.g., ```1```

Run the command:
```
python3 db_rep_gen.py  
```

## 3.3 MS spectral library query 
After building TopLib, a user can execute a basic library query. This feature enables users to extract a subset of spectral library for targeted retrieval of specified project.

* Input parameter:
  * project ID: e.g., ```1```
  * representative method: e.g., ```Single```
    
Run the command:
```
python3 db_query.py  
```
After running this command, an msalign file and a TSV file will be generated.

## 3.4 Building a spectral library  


## Contact
For more information you can visit our code website or send an email to kil7@tulane.edu.
