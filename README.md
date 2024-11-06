# TopLib 

TopLib is a Python package for building and searching spectral libraries for
top-down mass spectrometry (MS).  

### Requirements
* Python v>=3.8

If you use Toplib in your work, please cite the following publication:

* Kun Li, Haixu Tang, and Xiaowen Liu. TopLib: Building and searching top-down mass spectral libraries for proteoform identification (2024) bioRxiv preprint.
 

### 1. Preprocessing top-down MS data files 

MsConvert and TopPIC suite are needed for preprocessing top-down MS data files.

### 1.1. Download software  

* msConvert ([download](https://proteowizard.sourceforge.io/download.html))
* TopPIC Suite [>=v1.7.8](https://www.toppic.org/software/toppic/index.html)

### 1.2. Top-down mass spectral preprocessing 

In data preprocessing, raw MS files are converted into centroided mzML files using msConvert, 
mzML files are deconvoluted using TopFD, deconvoluted mass spectra are searched against 
a protein sequence database for proteoform identification using TopPIC, and
finally proteoform identifications are filtered using a Python script. Suppose
the raw MS file is spectra.raw and the protein sequence database is
proteins.fasta. Below is an example for data preprocessing.    

* Use msConvert to convert the raw MS data file spectra.raw to a centroided mzML file spectra.mzML 
* Use TopFD (https://www.toppic.org/software/toppic/) to deconvolute the mzML file spectra.mzML to generate an msalign file spectra_ms2.msalign
* Use TopPIC (https://www.toppic.org/software/toppic/) to search the msalign file spectra_ms2.msalign to the protein sequence database proteins.fasta for proteoform identification. The resulting files spectra_ms2_toppic_prsm_single.tsv and spectra_ms2_toppic_proteoform_single.tsv are used in the next step.  
* Filter out inconsistent proteoform identifications reported by TopPIC 

  Run the command: 
  ```
  python3 tsv_file_processing.py spectra_ms2_toppic_prsm_single.tsv spectra_ms2_toppic_toppic_proteoform_single.tsv
  ```
The resulting file is spectra_ms2_toppic_prsm_single_filtered.tsv


## 2. Building a spectral library from a top-down MS data set 
This enables users to build spectral representatives from their own data using a msalign file and a TSV file. This feature allows users to create a small, customized library based on their spectral data for quick review and analysis. 

* Input parameter:
  * msalign file: e.g., ```sw480_combined_ms2.msalign```
  * tsv file: e.g, ```sw480_combined_ms2_toppic_prsm_single_filtered.tsv```
  * representative method (choose one of the options by entering the corresponding number): ```1 = Single; 2 = Average.```
 
output: 
  * spectral representatives from input msalign file and tsv file, stored as a .db file with the name as ```sw480_combined_ms2.db```.    

Run the command: 
```
python3 ms_library_building.py
```

In addition to build your own library, a well-built TopLib MS spectral library based on SW480 and SW620 datasets has been established in TopLib. This library is integrated into the TopLib package to generate a representative mass spectrometry database, facilitating query result analysis for top-down proteomics studies. 
This spectral library was written in txt and db format and can be downloaded [download](http://127.0.0.1:5500/index.html) for use in TDP studies. The TopLib library includes six library spectra datasets: SW480-2D-1, SW480-2D-2, SW480-2D-3, SW620-2D-1, SW620-2D-2 and SW620-2D-3. 

## 3. Query extraction
This allows users to query spectra for a generated spectral representative file based on their input parameters. 

* Input parameter:
  * library name: e.g.,```sw480_combined_ms2.db```
  * msalign file for querying: e.g.,```sw620_combined_ms.msalign```
  * precurosr mass error tolerance setting (choose one of the options by entering the corresponding number): ```1 = ppm; 2 = Da.```
  * charge usage: whether to consider charge state when querying (choose one of the options by entering the corresponding number): ```1 = Yes; 2 = No.```

output: 
  * An TSV file containing the query results saved in the ```toplib_output```folder, e.g., ```query_res.tsv```

Run the command: 
```
python3 ms_library_query.py 
```

## 4. Comprehensive MS spectral library building
### 4.1 TOPLib library creation
This creates TopLib database tables.

Run the command: 
```
python3 db_gen.py 
```

### 4.2 Datasets addition
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

## 5 Specified MS spectral library query and extraction
After building TopLib, a user can execute a basic library query. This feature enables users to extract a subset of spectral library for targeted retrieval of specified project.

* Input parameter:
  * project ID: e.g., ```1```
  * representative method: e.g., ```Single```
    
Run the command:
```
python3 db_query.py  
```
After running this command, an msalign file and a TSV file will be generated.


## Contact
For more information you can visit our code website or send an email to kil7@tulane.edu.
