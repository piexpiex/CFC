python CFC_configuration/python_scripts/setup.py

fichero=$1

for imagen in $fichero/*.fits; do

	python CFC_configuration/python_scripts/size.py $imagen

	sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor1.sex


	psfex CFC_configuration/sextractor_result_files/test_psf.cat -c CFC_configuration/sextractor_configuration_files/psfex_config.psfex

	sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor2.sex

	python CFC_configuration/python_scripts/sex_analisis.py $imagen

done
