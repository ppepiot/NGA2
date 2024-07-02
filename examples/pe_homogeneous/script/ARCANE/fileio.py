""" General I/O file management

Contains all generic routines to compare, check, create and remove
files and folders

"""

import os
import shutil
from pathlib import Path

import ARCANE.display as display

logger = display.Logger()
logger.set_log('logCompute')


def create_dir(myfolder, overwrite=False):
    """Clean up and create folder where all simulations will be stored. Check if directory `myfolder` exists,
    if not create it. If it exists, stop, unless optional flag allows overwriting

    **Note**: if `myfolder` is a file rather than a folder, save it with a different name first.

    Parameters
    ----------
    myfolder : str
        name of folder to store all simulations
    overwrite : bool
        if True, existing folder will be overwritten (Default value = False)

    Returns
    -------
    True if successful


    Created: 11/12/17 [PP]

    Last modified: 11/12/17 [PP]


    TODO
    ----
    - Write unit test for function `create_dir()`
    - can we use only `os` instead of `os` and `pathlib`? (or only `pathlib`)

    """

    if Path(myfolder).exists() and not Path(myfolder).is_dir():
        logger.info("Warning creating case database: ', myfolder, ' is a file, not a folder! ")
        logger.info("... Moved to ', myfolder.strip(), '_old")
        movefile = myfolder + '_old'
        os.rename(myfolder, movefile)

    try:
        if not Path(myfolder).exists():
            # Folder does not exist, just create it
            os.mkdir(myfolder)
            return True
        elif Path(myfolder).exists() and overwrite is True:
            # Folder exists, but should be overwritten
            shutil.rmtree(myfolder)
            os.mkdir(myfolder)
            return True
        else:
            # If folder already exists and is not to be overwritten
            raise TypeError

    except TypeError:
        logger.debug("Folder {0} already exists and will not be overwritten: " \
                     "data will be re-used.".format(myfolder))


def cp_file(myfile, myfolder):
    """Shorthand notation to copy files

    Parameters
    ----------
    myfile : str
        file to be copied
    myfolder : str
        folder to copy to

    Returns
    -------
    True if successful

    TODO
    ----
    - Write unit test for function `cp_file()`
    """

    # Copy
    try:
        shutil.copy(myfile, myfolder)
    except SystemError:
        logger.error("{0} was not successfully copied to folder {1}".format(myfile, myfolder))
    return True


# Test functions directly from module
if __name__ == "__main__":
    create_dir('dummy_folder', overwrite=False)
