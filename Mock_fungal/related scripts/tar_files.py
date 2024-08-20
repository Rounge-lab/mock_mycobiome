# -*- coding: utf-8 -*-
"""

Tar all fastq files and remove unarchived ones

"""

import tarfile as tf
import os

wdir='/FULL_PATH_TO/Mock_fungal/data/simdata/'
files=os.listdir(wdir)
files=[f for f in files if '.tar.gz' not in f]

for f in files:
    if os.path.isfile(os.path.join(wdir,f+'.tar.gz')):
       os.remove(os.path.join(wdir,f))
    else:
        with tf.open(os.path.join(wdir,f+'.tar.gz'),'w:gz') as tar:
            tar.add(os.path.join(wdir,f),arcname=f)
            os.remove(os.path.join(wdir,f))
          

