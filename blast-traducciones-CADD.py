# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 17:16:38 2021

@author: android-2d
"""

from Bio.Blast import NCBIWWW,NCBIXML
from Bio import SeqIO
from io import StringIO

#leo un archivo fasta
f_record = next(SeqIO.parse("m_cold.fasta", "fasta"))

#busco en la BD
print("corriendo blast recuperando datos...")
result_handle = NCBIWWW.qblast("blastn", "nr", f_record.format("fasta"))

#guardo resultados
with open("m_cold_blast.out", "w") as save_file:
    blast_results = result_handle.read()
    save_file.write(blast_results)
result_handle.close()

print("Analizando los resultados y extrayendo información...")
# opción 1 - abre el archivo guardado para analizarlo
# opción 2 - crear un identificador a partir de la cadena y analizarlo

string_result_handle = StringIO(blast_results)
b_record = NCBIXML.read(string_result_handle)

# ahora obtengo la información de alineación 
# para todos los valores e mayores que algún umbral
E_VALUE_THRESH = 0.1
for alignment in b_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****Alineacion****")
            print("secuencia: %s" % alignment.title)
            print("longitud: %i" % alignment.length)
            print("e valor: %f" % hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")

"""
Requeriento 2 - Aliniar 2 secuenciAS
"""
seq1 = 'LIFAGKQLEDGRTLS'
seq2 = 'QLIFAAPKQLPGRT'

aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
alignments = aligner.align(seq1, seq2)
alignments = list(alignments)
print("Numero de aliniamientos: %d" % len(alignments))
alignment = alignments[0]
print("Score = %.1f" % alignment.score)
print("aliniacion:")
print(alignment)

"""
Requerimiento 3:Traducir secuencia de ARN mensajero a proteína
"""
#traducir arn mensajero a proteina
messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
proteina=messenger_rna.translate()
print ("3- messenger_rna:", messenger_rna)
print ("3- proteina:", proteina,"\n")

"""
Requerimiento 4:Traducir secuencia de ADN mensajero a proteína
"""
#traducir secuencia dna a proteina
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
proteina=coding_dna.translate()
messenger_rna=coding_dna.transcribe()
print ("4- coding_dna:", coding_dna)
print ("4- proteina:", proteina,"\n")
