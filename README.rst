VerMeerKAT
==========

RATT/RARG MeerKAT continuum self-calibration pipeline.

Usage
-----

Download and image observations in the last 3 days
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    vermeerkat

Specify a different configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    vermeerkat -c myconfig.conf

Download and image a specific observation file, using a custom configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    vermeerkat -f 123456789.h5 -c myconfig.conf

The Astronomer, by Vermeer
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: https://upload.wikimedia.org/wikipedia/commons/0/0e/Johannes_Vermeer_-_The_Astronomer_-_WGA24685.jpg
    :alt: The Astronomer
    :width: 500px
    :height: 500px
    :align: center

The latest version of the pipeline is depicted here. Unimplemented steps are shown in red:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: https://github.com/ska-sa/vermeerkat/blob/master/misc/Vermeerkat_flow.png
   :alt: Pipeline
