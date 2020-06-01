# CFC
Calibrador Fotométrico de CAFOS

Conjunto de subrutinas de bash y python que permiten realizar una calibración fotométrica automática de las imágenes del intrumento CAFOS del telescopio de 2.2m del observatorio de Calar Alto (CAHA)

# Requisitos

-*[SExtractor][3]* (Source-Extractor) versión 2.25.0 o mayor.

-*[PSFEx][3]* (PSF Extractor) versión 3.17.1 o mayor

-*[Astroquery][3]*.

-*[Astropy][4]*.

[1]: https://github.com/astromatic/sextractor
[2]: https://www.astromatic.net/software/psfex
[3]: https://astroquery.readthedocs.io/en/latest/
[4]: https://www.astropy.org/

# Funcionamiento

El programa funciona depositando el fichero CFC.sh, la carpeta CFC_configuration y una carpeta con las imágenes a calibrar (por ejemplo "files") y utilizando el siguiente comando en la terminal de UNIX "sh CFC.sh files".

Las imágenes a calibrar deben estar reducidas y calibradas astrometricamente, para ello se recomienda utilizar *[filabres][4]*.

[5]: https://github.com/nicocardiel/filabres

## Procesado de las imágenes

-Estimación de las dimensiones de la imagen y aplicación de una máscara e identificación del filtro utilizado (de momento unicamente funciona correctamente para filtros SDSS).

-Utilización de SExtractor y PSFEx para la identificación de los objetos en la imagen y la estimación de sus parámetros.

-Aplicación de diferentes criterios de slección para la identificación de objetos reales y descarte de artefactos.

-Selección de los objetos de mayor calidad fotométrica y comparación con los catálogos SDSS DR12 o APASS DR9.

-Elaboración de un catálogo asociado a cada imagen con los objetos calibrados en magnituds y diferentes parámetros de calidad.

# Resultados

-Catálogos de objetos calibrados en magnitud asociados a cada imagen de CAFOS
