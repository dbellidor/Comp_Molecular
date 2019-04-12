cad_dna="ATGGAAGTATTTAAAGCGCCACCTATTGGGATATAAG"
lista=list()
ini=0
fin=3
while (ini and fin) < len(cad_dna):
	temp=cad_dna[ini:fin]
	lista.append(temp)
	ini=fin
	fin=ini+3
	
print("Secuencia de ADN \n" + cad_dna )
#print(lista)
dna_tabla={"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
	   "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" }
lista_code=list()
for cad in lista:
	if dna_tabla.has_key(cad):
		temp=dna_tabla.get(cad)
		if temp!="STOP":
			lista_code.append(temp)
print ("Secuencia de Proteina: ")
print ''.join(lista_code)


