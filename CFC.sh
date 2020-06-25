 #!/bin/bash

python CFC_configuration/python_scripts/setup.py

arg=$1

case  $arg in
	*.csv)
		format="csv"
		line=0
		while IFS=',' read -r f1 f2 f3 f4
			
		do 
			if [ $line -gt 0 ];
				then 
					file=${f4#* } 
					folder=${f3#* } 
					imagen=science-imaging/$folder/$file

						python CFC_configuration/python_scripts/size.py $imagen

						sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor1.sex


						psfex CFC_configuration/sextractor_result_files/test_psf.cat -c CFC_configuration/sextractor_configuration_files/psfex_config.psfex

						sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor2.sex

						python CFC_configuration/python_scripts/sex_analisis.py $imagen $arg




			fi
			line=$((line + 1))
		done < "$arg"
		;;
	*)		format="folder"

			for imagen in $arg/*.fits; do

				python CFC_configuration/python_scripts/size.py $imagen

				sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor1.sex


				psfex CFC_configuration/sextractor_result_files/test_psf.cat -c CFC_configuration/sextractor_configuration_files/psfex_config.psfex

				sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor2.sex

				python CFC_configuration/python_scripts/sex_analisis.py $imagen 

			done

		;;
esac


