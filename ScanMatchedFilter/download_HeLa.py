# se importa wget para los archivos
import wget

# se nombran de acuerdo a la table que hicimos (izquierda) y a la derecha el link de descarga
ENCFF388WMD = "https://www.encodeproject.org/files/ENCFF388WMD/@@download/ENCFF388WMD.bigWig"
ENCFF255NQJ = "https://www.encodeproject.org/files/ENCFF255NQJ/@@download/ENCFF255NQJ.bigWig"
ENCFF391MRG = "https://www.encodeproject.org/files/ENCFF391MRG/@@download/ENCFF391MRG.bigWig"
ENCFF045NNJ = "https://www.encodeproject.org/files/ENCFF045NNJ/@@download/ENCFF045NNJ.bigWig"
ENCFF573VUK = "https://www.encodeproject.org/files/ENCFF573VUK/@@download/ENCFF573VUK.bigWig"
ENCFF741QNA = "https://www.encodeproject.org/files/ENCFF741QNA/@@download/ENCFF741QNA.bigWig"
ENCFF010PHG = "https://www.encodeproject.org/files/ENCFF010PHG/@@download/ENCFF010PHG.bigWig"
ENCFF444SGK = "https://www.encodeproject.org/files/ENCFF444SGK/@@download/ENCFF444SGK.bigWig"
ENCFF118MMT = "https://www.encodeproject.org/files/ENCFF118MMT/@@download/ENCFF118MMT.bigWig"
ENCFF689TMV = "https://www.encodeproject.org/files/ENCFF689TMV/@@download/ENCFF689TMV.bigWig"
ENCFF866KTJ = "https://www.encodeproject.org/files/ENCFF866KTJ/@@download/ENCFF866KTJ.bigWig"
ENCFF526LYS = "https://www.encodeproject.org/files/ENCFF526LYS/@@download/ENCFF526LYS.bigWig"


# se nombra la herramienta junto con el nombre del archivo y el lugar en donde se guarda
wget.download(ENCFF388WMD, '/home/daniel/TFG/ENCFF388WMD.bigWig')
wget.download(ENCFF255NQJ, '/home/daniel/TFG/ENCFF255NQJ.bigWig')
wget.download(ENCFF391MRG, '/home/daniel/TFG/ENCFF391MRG.bigWig')
wget.download(ENCFF045NNJ, '/home/daniel/TFG/ENCFF045NNJ.bigWig')
wget.download(ENCFF573VUK, '/home/daniel/TFG/ENCFF573VUK.bigWig')
wget.download(ENCFF741QNA, '/home/daniel/TFG/ENCFF741QNA.bigWig')
wget.download(ENCFF010PHG, '/home/daniel/TFG/ENCFF010PHG.bigWig')
wget.download(ENCFF444SGK, '/home/daniel/TFG/ENCFF444SGK.bigWig')
wget.download(ENCFF118MMT, '/home/daniel/TFG/ENCFF118MMT.bigWig')
wget.download(ENCFF689TMV, '/home/daniel/TFG/ENCFF689TMV.bigWig')
wget.download(ENCFF866KTJ, '/home/daniel/TFG/ENCFF866KTJ.bigWig')
wget.download(ENCFF526LYS, '/home/daniel/TFG/ENCFF526LYS.bigWig')