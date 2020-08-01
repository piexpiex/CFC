 #!/bin/bash

python CFC_configuration/python_scripts/setup.py $verbosity

arg=$1

if [ $overwrite = "no" ];
then
	case  $arg in
		*.csv)
			format="csv"
			line=0
			while IFS=',' read -r f1 f2 f3 f4
			
			do 
				if [ $line -gt 0 ];
					then 
						bh=${f4#* } 
						fh=${f3#* } 
						imagen=science-imaging/$fh/$bh

							python CFC_configuration/python_scripts/sex_analisis_saved.py $imagen $arg 


				fi
				line=$((line + 1))
			done < "$arg"
			;;
		"merge")
			python CFC_configuration/python_scripts/MERGER.py
			;;
		"MERGE")
			python CFC_configuration/python_scripts/MERGER.py
			;;
		*)		
			for imagen in $arg/*.fits; do

				python CFC_configuration/python_scripts/sex_analisis_saved.py $imagen 

			done
		;;
	esac
else
	case  $arg in
		*.csv)
			format="csv"
			line=0
			while IFS=',' read -r f1 f2 f3 f4
			
			do 
				if [ $line -gt 0 ];
					then 
						bh=${f4#* } 
						fh=${f3#* } 
						imagen=science-imaging/$fh/$bh

							python CFC_configuration/python_scripts/size.py $imagen

							sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor1.sex

							
							psfex CFC_configuration/sextractor_result_files/test_psf.cat -c CFC_configuration/sextractor_configuration_files/psfex_config.psfex

							sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor2.sex

							python CFC_configuration/python_scripts/sex_analisis.py $imagen $arg




				fi
				line=$((line + 1))
			done < "$arg"
			;;
		"merge")
			python CFC_configuration/python_scripts/MERGER.py
			;;
		"MERGE")
			python CFC_configuration/python_scripts/MERGER.py
			;;
		*)		
			format="folder"

			for imagen in $arg/*.fits; do

				python CFC_configuration/python_scripts/size.py $imagen

				sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor1.sex


				psfex CFC_configuration/sextractor_result_files/test_psf.cat -c CFC_configuration/sextractor_configuration_files/psfex_config.psfex

				sex $imagen -c CFC_configuration/sextractor_configuration_files/sextractor2.sex

				python CFC_configuration/python_scripts/sex_analisis.py $imagen 

			done

		;;
	esac
fi


