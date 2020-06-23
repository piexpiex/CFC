Operation
*********

Operation modes
===============

CFC has two operation modes, one for images obtained using filabres and other for images obtained by own reduction and astrometric calibration (it also works with filabres images) deposited in a folder. 

filabres operation mode
-----------------------

Calibration of images obtained by filabres requires to deposite the CFC_configuration and CFC.sh files in the same directory as the resulting filabres folders and a csv file with the information of the selected images and their location in the filabres estructure (for this example it has the name "filabres_tree.csv")

The operating command is:

.. code-block:: bash 

   sh CFC.sh filabres_tree.csv

Own calibrated images operation mode
------------------------------------

Calibration of images with own reduction and astrometric calibration only requires to deposite the CFC_configuration and CFC.sh files in the same directory as the folder with the images (for this example it has the name "files").

The operating command is:

.. code-block:: bash 

   sh CFC.sh files

