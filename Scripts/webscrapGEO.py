from bs4 import BeatifulSoup
import codecs
import os

os.chdir('/home/biouser/Desktop/MariaGC/entornoR/filterGEO/adipose')

html_file = "adipose-source.html"
f = codecs.open(html_file, 'r', 'utf-8')

document = BeatifulSoup(f.read())

GSEID_tags = document.find_all('div', class_ ='resc')
GSEType_tags = document.find_all('dl', class_='details lefty')
GSEPlatform_tags = document.find_all('dl', class_='details')

GSEID_output = []
GSEType_output = []

for tag in GSEID_tags:
    accesion = tag.find_all('dd')
    GSEID_output.append(accesion[0].text)

for i in range(len(GSEType_tags)): # len(GSEType_tags) = 798 (399*2),tenemos 2 tags details lefty, organismo y tipo de estudio, y nos quedamos con las dd, dd1 = hsa, dd2= study type
    tag = GSEType_tags[i]
    tag_content = tag.find_all('dd')
    j= i+1
    if (j%2) == 0:
        studytype = tag_content[0].text
        GSEType_output.append(studytype)

GSEGPL_output = []
GSEsamples_output = []
contador = 0
n_gse_w_more_gpls = 0
indicador = False
for i in range(len(GSEPlatform_tags)):
    contador += 1
    if contador == 3:
        tag_gpl = GSEPlatform_tags[i]
        gpl_content = tag_gpl.find_all('dd')
        GPL = gpl_content[0].text
        if ('GDS' in GPL):
            GPL = gpl_content[1].text
        if ('GPL' not in GPL): # then is related platform or dataset
            gpl_content = tag_gpl.find_all('dt')
            GPL = gpl_content[0].text
            if ('Datasets' in GPL):
                GPL = gpl_content[1].text
            n_gse_w_more_gpls +=1 # related platforms > 4
            indicador = True # Cambia la etiqueta html
        #print(GPL)
        pos_final = GPL.rfind(' ')
        GPL = GPL[:pos_final]
        if 'GPL' in GPL: # else -> related platforms
            GPL = GPL.replace(' ', ',')
        # ltrim, rtrim, separar por comas
        contador = 0
        GSEGPL_output.append(GPL)
        if indicador  == False:
            gsm_content = tag_gpl.find_all('dt')
            n_gsm = gsm_content[0].text
        if indicador == True: # cambia la etiqueta:
            n_gsm = gpl_content[1].text
            indicador = False # restore default value
        print(n_gsm)
        GSEsamples_output.append(n_gsm)



print(len(GSEID_output))
print(len(GSEType_output))

# Generate output:
out_txt = 'adipose_GSEs.txt'
with open(out_txt, 'w') as f_out:
    f_out.write('GSEID\tGPL\tStudyType\n')
    for k in range(len(GSEID_output)):
        #GSE = (GSEID_output[k])
        #GSEt = (GSEType_output[k])
        f_out.write('%s\t%s\t%s\n'%(GSEID_output[k], GSEGPL_output[k], GSEType_output[k]))
        GSEt = GSEType_output[k]
        print(GSEt)

