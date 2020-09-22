 #!/bin/bash

echo --------------------------------
echo --- Starting the CFC anaysis ---
echo --------------------------------



arg=$1

case  $arg in
	"merge")
		merge_statur="merge"
		;;
	"MERGE")
		merge_statur="merge"
		;;	
	*)
		python CFC_configuration/python_scripts/setup.py $verbosity
		;;
esac

case $overwrite in
	no)
		overwrite="no"
		;;
	NO)
		overwrite="no"
		;;
	*)
		overwrite="yes"
		;;
esac

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

echo -------------------------
echo --- CFC anaysis ended ---
echo -------------------------
