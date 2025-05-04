# Scripts
# **foldseek_anno.py**    
Takes a Foldseek created m8 file, extracts all unique AlphaFold UniProt IDs, queries the AlphaFold API to fetch protein descriptions, and then produces a new tsv file with existing columns with added protein potential functional description.  

Requirements:  
Python ≥ 3.7  
aiohttp ```pip install aiohttp```  
tqdm	```pip install tqdm```  
Internet connection  

run as:  
python foldseek_anno.py \  
-i path/to/foldseek.m8 \  
-o path/to/foldseek.m8.annotation.tsv  
