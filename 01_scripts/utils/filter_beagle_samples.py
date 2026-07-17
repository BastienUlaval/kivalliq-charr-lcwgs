# -*- coding: utf-8 -*-
#!/usr/bin/env python
import sys
import gzip

def get_sample_indices(bam_list_file, subset_bam_list_file):
    """Retourne les indices des echantillons a conserver (0-based)"""
    # Lire la liste complete
    with open(bam_list_file, 'r') as f:
        all_samples = [line.strip() for line in f]
    
    # Lire la liste des echantillons a conserver
    with open(subset_bam_list_file, 'r') as f:
        subset_samples = [line.strip() for line in f]
    
    # Trouver les indices
    indices = []
    for subset_sample in subset_samples:
        try:
            idx = all_samples.index(subset_sample)
            indices.append(idx)
        except ValueError:
            print(f"Attention: echantillon non trouve: {subset_sample}", file=sys.stderr)
    
    return sorted(indices)

def filter_beagle(input_file, output_file, sample_indices):
    """Filtre le fichier beagle pour ne garder que les echantillons specifies"""
    
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        for line_num, line in enumerate(infile):
            fields = line.strip().split('\t')
            
            if line_num == 0:  # Header
                # Garder les 3 premieres colonnes (marker, allele1, allele2)
                new_fields = fields[:3]
                # Ajouter les colonnes des echantillons selectionnes (3 colonnes par echantillon)
                for idx in sample_indices:
                    start_col = 3 + (idx * 3)
                    new_fields.extend(fields[start_col:start_col+3])
            else:
                # Meme traitement pour les donnees
                new_fields = fields[:3]
                for idx in sample_indices:
                    start_col = 3 + (idx * 3)
                    new_fields.extend(fields[start_col:start_col+3])
            
            outfile.write('\t'.join(new_fields) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_beagle_samples.py <bam_list> <subset_bam_list> <input_beagle>")
        sys.exit(1)
    
    bam_list = sys.argv[1]
    subset_bam_list = sys.argv[2] 
    input_beagle = sys.argv[3]
    output_beagle = input_beagle.replace('.beagle.gz', '.subset.beagle.gz')
    
    print(f"Filtrage de {input_beagle} vers {output_beagle}")
    
    # Obtenir les indices des echantillons
    sample_indices = get_sample_indices(bam_list, subset_bam_list)
    print(f"Echantillons a conserver (indices): {sample_indices}")
    print(f"Nombre d'echantillons: {len(sample_indices)}")
    
    # Filtrer le fichier beagle
    filter_beagle(input_beagle, output_beagle, sample_indices)
    
    print(f"Fichier filtre cree: {output_beagle}")
