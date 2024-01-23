# Projekt z Genomiki obliczeniowej


## Skrypty
 - `download.py` zawiera klasę służącą do pobierania i przechowywania proteomów.
 - `pipeline.py` uruchamia pipeline
 - `util.py` klasy służące do obsługi różnych narzędzi które są wykorzystywane w pipeline'ie
 - `calculate_RF.py` Służy do wyliczania odległości RF pomiędzy dwoma drzewami, jedno z których zawiera się w drugim.

### Sposoby użycia
#### `pipeline.py`
Potrzebuje zainstalowanych:
    
    - MAFFT
    - mmseqs2
    - fasturec
    - rapidNJ
    - dendropy
    - Biopython
Przyjmuje, że MAFFT i mmseqs2 są w PATH. Ścieżki do pozostałych narzędzi należy uzupełnic pod importami.

Włączenie z opcją `-h` powoduje wyświetlenie poniższej wiadomości:
``` 
usage: Supertrees project pipeline [-h] [--output_root OUTPUT_ROOT] [-d] [-v] [-n N] [--min_seq MIN_SEQ] [--min_org MIN_ORG] [--min_seq_id MIN_SEQ_ID] [--bootstrap BOOTSTRAP] [--logfile LOGFILE] input_file

positional arguments:
  input_file            Input file containing names of proteomes to run analysis on as well as urls to ftp containing proteomes.

options:
  -h, --help            show this help message and exit
  --output_root OUTPUT_ROOT
                        Output root directory.
  -d, --debug           Print lots of debugging statements
  -v, --verbose         Be verbose
  -n N                  Number of proteomes to analyse
  --min_seq MIN_SEQ     Minimum number of sequences in cluster
  --min_org MIN_ORG     Minimum number of organisms in a cluster
  --min_seq_id MIN_SEQ_ID
                        Minimum identity of sequence in clusters
  --bootstrap BOOTSTRAP, -b BOOTSTRAP
                        Number of bootstrap replicates to generate the trees.
  --logfile LOGFILE     Logfile in which to write

```

Do wyliczenia danych w prezentacji użyłem:
` python pipeline.py input.csv --min_seq 5 --output_root ./data/ --min_seq_id 0.5`

W wersji bez bootstrapu, oraz:

` python pipeline.py input.csv --min_seq 5 --output_root ./data/ --min_seq_id 0.5 -b 100`

Do wyliczenia drzew z bootstrapem.

#### `calculate_RF.py`
Wylicza odległość Robinsona - Fouldsa dla dwóch drzew, z których zbiór liści jednego jest podzbiorem zbioru liści drugiego.

Użycie z flagą `-h` daje:

```
usage: calculate_RF.py [-h] input_1 input_2

positional arguments:
  input_1     path to file containing one of the trees
  input_2     path to file containing tree that is a subtree of input1

options:
  -h, --help  show this help message and exit
```


## Wyniki
Zob. zawartość `./results/` oraz `prezentacja.pdf`