# Welcome to Dalotia Coriaria DB Parser
<p align="center">
<img src="https://github.com/klcolon/DalotiaCoriaria_DB/blob/main/icon/beetle.png" alt="beetle icon" width="200" height="200">
</p>

## Create mRNA FASTA file from Dalotia coriaria genome to be fed into Probe Wizard 2.0
This app is used to generate mRNA fasta files that will be used as input for probe wizard 2.0 in designing fluorescent in situ hybridization (FISH) probes from Dalotia coriaria genome. Probes targeting exonic, intronic, or whole mRNA can be generated using this app. Probe wizard 2.0 is a propriety app of Cai Lab at Caltech.

## General Dependencies 
- python >= 3.7.0
- pandas >= 1.1.3
- biopython >= 1.79

## Make an executable app
1) Install pyinstaller `pip install pyinstaller`.
2) Go to DC_DB_Parser directory and run the following ```pyinstaller --icon=../icon/beetle.png --onefile DalotiaCoriariaDB.py```. This will create the app which should be compatible with most OS.
3) If you have a linux system, the --icon flag will not work. Instead, edit the paths in the provided dalotia_db.sh and dalotia.desktop files. The .desktop file should be placed into `~/.local/share/applications/` directory.

<p align="center">
<img src="https://github.com/klcolon/DalotiaCoriaria_DB/blob/main/icon/gui.png" alt="gui" width="400" height="200">
</p>

## Run DalotiaCoriariaDB.py script 
1) Go to DC_DB_Parser directory and install the neccessary dependencies if you wish to run DalotiaCoriariaDB.py `pip install -r requirements.txt`.
2) Run the following to open the app `python3 DalotiaCoriariaDB.py`

## How to use this app?
1) Load FASTA file, GFF file, and your desired gene list as a CSV file.
2) Select desired Biotype ('mRNA', 'exons', 'introns'). If you choose 'exons' or 'introns' then a BED file must be provided. If you choose 'mRNA' then it will use the entire strand.
3) Click the Generate Fasta button. A window will open to define where you wish to save the outputs.
4) Once the FASTA file is created, feed it into Probe Wizard 2.0 using your desired settings.
5) After probe genereration, open the app again load the probes list.
6) Click the BLAST button and you will get your final probes list.

Link: https://github.com/klcolon/DalotiaCoriaria_DB
