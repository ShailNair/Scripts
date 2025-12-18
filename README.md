# Scripts
# **foldseek_anno.py**    
Takes a Foldseek created m8 file, extracts unique Alphafold, PDB or ESMAtlas matched IDs, Rretrieves annotations via API (AlphaFold, PDB) or web scraping (MGnify) and and then produces a new tsv file with existing columns and added protein potential functional description.  


Requirements:  
Python ≥ 3.7  
aiohttp ```pip install aiohttp```  
tqdm	```pip install tqdm```  
Internet connection  

run as:  
python foldseek_anno.py \
-i path/to/foldseek.m8 \
-o path/to/foldseek.m8.annotation.tsv \
-d alphafold/mgnify/pdb
For full arguments:

python foldseek_anno.py -h


