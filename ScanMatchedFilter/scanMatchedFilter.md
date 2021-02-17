# scanMatchedFilter

Introducción de librerías necesarias

•	Matplotlib ---> es una biblioteca completa para crear visualizaciones estáticas, animadas e interactivas en Python.

•	Matplotlib.use('Agg') ---> solo procesa archivos PNG.

•	Sys ---> El módulo sys proporciona información sobre constantes, funciones y métodos del intérprete de Python.

•	Os ---> proporciona funciones sencillas que nos permiten interactuar y obtener información del sistema operativo e incluso controlar los procesos hasta un límite.

•	Un orderedDict es una subclase de diccionario que recuerda el orden en que se agregan sus contenidos.

•	El conjunto de programas BEDTools se utiliza ampliamente para la manipulación del intervalo genómico o "álgebra del genoma". Pybedtools envuelve y amplía BEDTools 
  y ofrece manipulaciones a nivel de función desde Python.
  
•	El objetivo de metaseq es unir una gran cantidad de software existente en un marco para explorar datos genómicos.

•	NumPy se utiliza para trabajar con matrices.

•	Scipy es el paquete científico, y de él importamos interpolar.

•	SciPy.optimize proporciona funciones para minimizar (o maximizar) funciones objetivas, posiblemente sujetas a restricciones.

•	Un backend de PDF matplotlib.

•	pylab proporciona una interfaz orientada a objetos a la biblioteca de trazado Matplotlib cercano a Matlab™.

•	testData ---> Un paquete simple que genera datos para pruebas.


Función principal del programa.   def main

Es la función principal del programa, y es la que se ejecuta nada más comenzar. Tiene como entradas:
-	bigWigFileList
-	metaProfileList
-	chrNameList
-	peakFileList
-	trainingPositives
-	trainingNegatives
-	opPrefix


Lo primero que hace es leer el nombre de los cromosomas y su longitud, para ello abre chrNameFile en modo lectura, y en la variable chromosomes acumula lo que devuelve 
la función readChromosomesAndSize (en concreto devuelve chrSize). Para finalizar se cierra lo que estaba abierto en modo lectura.

readChromosomesAndSize lee toda la lista del cromosoma línea por línea almacena en name el valor de la primera columna y en size el valor de la segunda columna. 
Posteriormente forma la variable chrSize de cada nombre con el tamaño de cada cromosoma y la devuelve.

El segundo paso es leer los meta-perfiles. Para ello se abre metaProfileList en modo lectura, se lee línea a línea y se genera la variable metaprofiles, formada por lo que 
devuelve la función readMetaprofiles (en concreto devuelve metaprofiles). Por último, se cierra la lectura.

readMetaprofiles lee toda la lista de metaprofiles línea por línea y almacena en sygnalType el valor de la primera columna y en filename el valor de la segunda columna. 
Para finalizar crea una nueva variable metaprofiles de cada señal con su correspondiente nombre y la devuelve.

Lo tercero que realiza el programa principal es leer los archivos bigWig. Para ello se abre en modo lectura bigWigFileList. Se lee línea a línea y se acumula en bigWigList y 
en bigWigFiles lo que devuelve la función readBigWigList (en concreto devuelve bigWigList y bigWigFiles). Para finalizar se cierra la lectura.

readBigWigList lee línea por línea y acumula en sygnalType el valor de la primera columna y en filename el valor de la segunda columna. Posteriormente crea una variable 
bigWigList para cada tipo de señal (metaseq.genomic_signal(filename, "bigWig") esta función es propia de las librerías y supongo que hace algo relacionado con el genoma 
que desconozco) y crea otra variable bigWigFiles para cada tipo de señal con el correspondiente nombre, y devuelve las dos variables.

Lo siguiente que se realiza es una comprobación de cordura de la señal de marcas y meta-perfiles. Para ello se utiliza la función sameMarks. Esta función tiene como entrada 
bigWigList y metaprofiles. Comprueba si todas las marcas con señales tienen meta-perfiles asociados. Con un bucle for recorre una marca en bigWigList y lo compara con su 
correspondiente metaprofiles, si no coincide lo elimina, así solamente se queda con las marcas en común.

El siguiente paso es leer la lista de archivos de pico, para ello se abre peakFileList en modo lectura, lee línea a línea y almacena en peakFiles el valor de la segunda columna. 
Posteriormente se vuelve a utilizar la función sameMarks para escoger las marcas en común de peakFiles y de metaprofiles.  

A continuación, se crea la variable pp con la función PdfPages, que la entra opPrefix.

Pasamos a obtener señales y regiones relevantes para verificar el meta-perfil. 

En toda esta parte me pierdo un poco, ya que tiene gran parte de biología, y voy a poner lo que creo que hacen, aunque no estoy seguro del todo de que sea lo correcto.

Lo primero que hacen es buscar “H3K27ac” en todos los bigWigFiles. Si lo encuentra pasa a obtener las señales. Para ello recorre una señal en bigWigFiles y obtiene signalList y 
regionList de cada señal con la función getRelevantSignals (que tiene como entradas: bigWigFiles, bigWigList, currRegions y currSignal), y devuelve signalList y currRegions.    

getRelevantSignals recorre currRegions haciendo cálculos (supongo que los cálculos tienen su explicación biológica, la cual yo no entiendo, solamente veo que se realizan con 
el fin de conseguir algo) y va generando una señal con bigWigList de currSignal y un vector de currFeature. También añade al final de signalList el valor de la primera columna 
de signal. Cuando se termina de recorrer el bucle se elimina bigWigList[currSignal] para crear una nueva bigWigList[currSignal] con la función 
metaseq.genomic_signal(bigWigFiles[currSignal], "bigWig"). Para finalizar devuelve signalList y currRegions.
A continuación, se calculan las puntuaciones de filtro coincidentes para diferentes señales. Para ello se recorre currSignal en bigWigFiles. Se crea currWidth y se inicializa
a 350.  Mientras que currWidth es menor o igual de 1100 se ejecuta scoreWithMatchedFilter y a currWidth se le suma 25, hasta que deja de cumplirse la condición y salimos del
bucle. Al finalizar se borran regionList y signalList.
La función scoreWithMatchedFilter tiene como entradas signalList, metaprofile, currWidth, currSignal y opPrefix. Crea variables con tamaños definidos que se utilizarán para 
hacer cálculos, y crea una variable op (opPrefix + currSignal + currWidth pasado a string + la terminación .bed en modo escritura, para crear un archivo .bed). en esta 
función se obtienen valores de x e y, que se interpolan y que supongo que son para la representación gráfica. Al finalizar se cierra el archivo que estaba en modo escritura.

Tras obtener esto, se calculan las puntuaciones normalizadas con la función calculateNormalizedScores, que tiene como entradas peakFiles, currWidth, opPrefix y pp. Esta 
función aun no la comprendo completamente, pero supongo que es para obtener valores y representar, y creo que puede ser una función problemática a la hora de ejecutar el 
código y que tenga algo que ver con el fallo que dio.

La parte final es la correspondiente al entrenamiento y al SVM model, y para finalizar la muestra de los datos obtenidos de todo el programa. 




