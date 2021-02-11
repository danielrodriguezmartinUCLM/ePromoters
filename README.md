

scanMatchedFilter.py escanea todo el genoma con un filtro emparejado y aplica el modelo SVM utilizando los datos de entrenamiento proporcionados con el código para predecir potenciadores y promotores de una manera específica del tipo de célula.
python scanMatchedFilter.py

dónde:

	<fileList> es un archivo con la lista de señales de cromatina en el formato (pestaña de 2 columnas delimitada con el nombre del conjunto de datos experimental en la columna 1 y el nombre del archivo en la columna 2). Las señales de cromatina están en el enriquecimiento de la señal de la señal de cambio logarítmico sobre el control.
	<metaProfileList> es la lista con perfiles de entrenamiento. Este es un archivo de 2 columnas delimitado por tabuladores con la primera columna que contiene el nombre del conjunto de datos experimental (por ejemplo, H3K4me1) y la segunda columna que contiene el nombre del archivo con metaperfil (proporcionado en los metaperfiles del subdirectorio).
	<chrNameList> es la lista con los nombres de los cromosomas. La primera columna contiene el nombre del cromosoma y la segunda columna contiene la longitud del cromosoma.
	<peakFileList> es el archivo con picos de cromatina. Estas regiones se eliminan durante el ajuste del modelo de fondo.
	<positiveScores> es el archivo que contiene las puntuaciones de todos los aspectos positivos del entrenamiento, que se proporciona en el directorio de datos de entrenamiento.
	<negativeScores> es el archivo que contiene las puntuaciones de todos los negativos de entrenamiento, que se proporciona en el directorio de datos de entrenamiento.
	<opPrefix> es el prefijo para todos los archivos de salida.
	
El archivo de salida final test_SVMpredScores.dat contiene las puntuaciones y predicciones de SVM.
