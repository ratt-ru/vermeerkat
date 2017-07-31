**DEPRECATED PROTOTYPE: CHANGES TO BE INCORPORATED IN HI PIPELINE**
====================================================================

VerMeerKAT
==========

RATT/RARG MeerKAT continuum self-calibration pipeline.

Source Repository
-----------------

VerMeerKAT is available at `<https://github.com/ska-sa/vermeerkat>`_.

Installation
------------

Clone VerMeerKAT from github and install it in a virtual environment:

.. code:: bash

    $ cd /src
    $ git clone git@github.com:ska-sa/vermeerkat.git
    $ virtualenv ~/venv/vermeerkat
    $ source ~/venv/vermeerkat/bin/activate
    (vermeerkat) $ pip install pip -U
    (vermeerkat) $ pip install setuptools -U
    (vermeerkat) $ pip install /src/vermeerkat

VerMeerKAT depends on stimela_, which uses radio astronomy packages
installed in docker containers to perform a reduction.
The containers must built and installed.
Note that this step will involve downloading GB of docker images.

.. code:: bash

    (vermeerkat) $ stimela build

Usage
-----

In general, it is necessary for VerMeerKAT to access the SDP pysolr server
in order to query the latest observations and obtain the location of
observation hdf5 files for download.
As such it initially needs to run within the SKA SA network.
Once the observation metadata and the hdf5 file have been downloaded,
it is possible to rerun VerMeerKAT on an observation in an ``offline mode``
without access to these servers.
This mode requires specifying the specific observation file, as mentioned
in the use case below.

Initial setup
~~~~~~~~~~~~~

Make a directory in which you wish to perform your reduction and copy the default
configuration file into the ``input`` subdirectory.

.. code:: bash

    (vermeerkat) $ mkdir -p ~/pipelines/input
    (vermeerkat) $ cd ~/pipelines
    (vermeerkat) $ cp /src/vermeerkat/vermeerkat/conf/default.conf myconfig.conf

At present, three files referenced in ``myconfig.conf`` need to be present in the ``input`` folder

* A MeerKAT 4096 channel rfi mask. (``rfi_mask.pickle``)
* A firstpass aoflagger flagging strategy file (``29Dec2016_firstpass_strategy.rfis``)
* A secondpass aoflagger flagging strategy file (``06Jan2017_secondpass_strategy.rfis``)

Configuration File
~~~~~~~~~~~~~~~~~~

``default.conf`` is a YAML file containing many options for configuring
the various tasks within the pipeline.
Many of them correspond to options for the CASA tasks that they call.
As they are legion, they are not documented here, and you should have a
familiarity with CASA when modifying them.

Download and image observations in the last 3 days
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    (vermeerkat) $ vermeerkat

Specify a user configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    (vermeerkat) $ vermeerkat -c myconfig.conf

Download and image a specific observation file, as well as a custom configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    (vermeerkat) $ vermeerkat -f 123456789.h5 -c myconfig.conf

Other useful command line options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


-b, --bandpass-calibrator   Manually set the bandpass calibrator used
                            used to estimate the flux present in the observation.

-g, --gain-calibrator       Manually specify the gain calibrator for estimating the
                            gains during the observation.




The latest version of the pipeline is depicted here. Unimplemented steps are shown in red:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: https://github.com/ska-sa/vermeerkat/blob/master/misc/Vermeerkat_flow.png
   :alt: Pipeline

The Astronomer, by Vermeer
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: https://upload.wikimedia.org/wikipedia/commons/0/0e/Johannes_Vermeer_-_The_Astronomer_-_WGA24685.jpg
    :alt: The Astronomer
    :width: 500px
    :height: 500px
    :align: center

.. _stimela: https://github.com/SpheMakh/Stimela
.. _vermeerkat: https://github.com/ska-sa/vermeerkat
