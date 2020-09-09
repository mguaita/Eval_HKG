# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 13:16:44 2020

@author: Maria
"""

"""
Generamos una tabla relacional entre los estudios y las plataformas
1 estudio - 1 plataforma
La clave primaria es combinatoria de GSE-GPL
No se repite la combinación GSE-GPL
"""
# Listas de almacenamiento, se correlacionan por posición
col_GSEIDs = []
col_GPLs = []
excluded_GSEs = []

i=0
with open('adipose_GSEs.txt') as f:
    for line in f:
        print(i)
        print(line)
        line = line.split('\t') #[0]-GSEID, [1]-GPLs, [2]-Study_tipes
        GSE_ID = line[0]
        GPLs = line[1] # si hay >1 plastaforma, habran ','
        print(GSE_ID)
        print(GPLs)
        if ',' in GPLs: # estudio multiplataforma
            GPLs_aux = GPLs.split(',') # convertimos la cadena de GPLs en lista
            print(GPLs_aux)
            for gpl in GPLs_aux: 
                col_GSEIDs.append(GSE_ID) # Añadimos tantas veces el estudio como plataformas tenga
                col_GPLs.append(gpl)
        elif 'related' in GPLs:
            excluded_GSEs.append(GSE_ID)
        else: # estudios monoplataforma - se añade una sola vez
            col_GSEIDs.append(GSE_ID) # Añadimos tantas veces el estudio como plataformas tenga
            col_GPLs.append(GPLs)
        input('ENTER')
        i+=1

f_out = open('gse_gpl_relational.txt', 'w')
for j in range(len(col_GSEIDs)):
    gseid = col_GSEIDs[j]
    gpl = col_GPLs[j]
    f_out.write('%s\t%s\n'%(gseid,gpl))
