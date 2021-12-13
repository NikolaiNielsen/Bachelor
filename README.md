# Bachelor - Visualization of Concepts in Condensed Matter Physics
This repository stores my bachelors project in Physics from the University of Copenhagen, under the supervision of Mark Spencer Rudner.

The source code is located in the cmp folder, along with a rawtext installation guide (reproduced below) and a list of packages used (packages.txt).


## Usage
To use this package, you need a python installation with the packages and versions specified in the packages.txt document. This unfortunately means that you will not be able to use ERDA, for example, since it does not support the PyQT package which we use to create the graphical user interface. You will need a local installation of Python.

The easiest way to do this is, in my opinion, to use Anaconda.

Simply download and install Anaconda from the following link:
https://www.anaconda.com/distribution/#download-section

When prompted if you want to "prepend Anaconda to $PATH" or something along those lines, I urge you to press yes. What it does is make it so when you type in "python" in the command-line or terminal, it opens up the python interpreter from Anaconda. This will also enable you to use the "conda" commands in the command-line/terminal, which will be needed for the next step.

If you choose "no" on the above prompt, you will need to either append or prepend Anaconda yourself, so the conda commands will run.

When the distribution has been installed, open up a command-line / terminal, navigate to the folder with the scripts you have downloaded and type/paste in the following command:

    conda create -n cmp --file=packages.txt

This will create a new python "environment" called "cmp", with the correct version of both python and the packages needed. Next we activate the environment using the following command:

    conda activate cmp

or

    source activate cmp

or what-ever Anaconda prompts you to write. This will make it so when you use python, it uses the correct versions of the packages. To revert back to the base environment (the one installed with anaconda), you can either restart your terminal or write

    conda activate base

or what-ever the correct command was before (source activate base, etc).

Now you're ready to use the scripts! Simply type in the following command to run the scripts:

    python cmp.py


## Abstract
In this thesis we explore key concepts in the field of condensed matter physics and translate these into visualisations using the computer programming language Python. Concepts explored include the geometry of crystal structures, families of lattice planes, neutron scattering and the band structure of two dimensional materials. Four distinct programs are written to visualize each of these concepts, and they are bundled into one open-source software package available at https://github.com/NikolaiNielsen/Bachelor, along with documentation on its use.
