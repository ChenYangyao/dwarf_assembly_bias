# Copyright (C) 2025 Yangyao Chen (yangyaochen.astro@foxmail.com) - All Rights 
# Reserved
# 
# You may use, distribute and modify this code under the MIT license. We kindly
# request you to give credit to the original author(s) of this code, and cite 
# the following paper(s) if you use this code in your research: 
# - Zhang Z. et al. 2025. Nature ???, ???.

from __future__ import annotations
import numpy as np
from pathlib import Path
from pyhipp.core import abc, DataTable
from pyhipp.io import h5

class GalaxySample(abc.HasLog, DataTable):
    '''
    A simple galaxy sample.
    
    @data: attributes of the galaxies, e.g. {'ra': ..., 'dec': ...}.
    @verbose: whether to print log messages.
    @copy: whether to copy the data.
    
    Attributes:
    - n_objs: number of objects in the sample.
    '''
    def __init__(self, data: dict[str, np.ndarray], verbose=True, copy=True):

        keys = tuple(data.keys())
        n_objs = len(data[keys[0]])
        for k, v in data.items():
            assert len(v) == n_objs, f"Size of {k} != {n_objs}"
        
        super().__init__(data=data, verbose=verbose, copy=copy)
        
        self.n_objs = n_objs
        self.log(f'GalaxySample: {n_objs=}, {keys=}')

    @classmethod
    def from_file(cls, path: Path | str, **init_kw):
        '''
        Create the sample from a file. The file should be in HDF5 format,
        with datasets containing the galaxy attributes.
        
        @path: path to the file.
        @init_kw: additional keyword arguments passed to __init__().
        '''
        data = h5.File.load_from(path)
        return cls(data, **init_kw)