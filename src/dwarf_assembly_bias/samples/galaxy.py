# Copyright (C) 2025 Yangyao Chen (yangyaochen.astro@foxmail.com) - All Rights 
# Reserved
# 
# You may use, distribute and modify this code under the MIT license. We kindly
# request you to give credit to the original author(s) of this code, and cite 
# the following paper(s) if you use this code in your research: 
# - Zhang Z. et al. 2025. Nature 642, 47-52.

from __future__ import annotations
from typing import Self
import numpy as np
from pathlib import Path
from pyhipp.core import abc, DataTable
from pyhipp.io import h5

class GalaxySample(abc.HasLog, abc.HasDictRepr):
    
    repr_attr_keys = ('n_objs',)
    
    '''
    A simple galaxy sample.
    
    @data: attributes of the galaxies, e.g. {'ra': ..., 'dec': ...}.
    @verbose: whether to print log messages.
    @copy: whether to copy the data.
    '''
    def __init__(self, data: dict[str, np.ndarray], verbose=True, copy=True):

        keys = tuple(data.keys())
        n_objs = len(data[keys[0]])
        for k, v in data.items():
            assert len(v) == n_objs, f"Size of {k} != {n_objs}"
        
        super().__init__(verbose=verbose)
        self.data = DataTable(data=data, copy=copy)
        self.log(f'GalaxySample: {n_objs=}, {keys=}')
        
    def __getitem__(self, key: str | tuple[str, ...]):
        '''
        Get the attribute of the galaxy sample.
        
        @key: key of the attribute.
        '''
        return self.data[key]
    
    @property
    def n_objs(self):
        key = next(iter(self.keys()))    
        val = self[key]
        return len(val)
    
    def keys(self):
        return self.data.keys()
    
    def values(self):
        return self.data.values()
    
    def items(self):
        return self.data.items()
    
    def __contains__(self, key: str):
        return key in self.data
    
    def subset(self, args: np.ndarray | slice, copy=True) -> Self:
        '''
        Row-wise subsetting by `args`, applied to each value. Return a new 
        (copied) sample.
        
        @copy: if False, values are not guaranteed to be copied.
        '''
        return GalaxySample({
            k: v[args] for k, v in self.data.items()    
        }, copy=copy, verbose=self.verbose)
        
    def subset_by_val(self, key: str, lo=None, hi=None, eq=None) -> Self:
        sel = np.ones(self.n_objs, dtype=bool)
        val = self[key]
        if lo is not None:
            sel &= val >= lo
        if hi is not None:
            sel &= val < hi
        if eq is not None:
            sel &= val == eq
        return self.subset(sel, copy=False)
    
    def subset_if(self, key: str) -> Self:
        '''
        Equivalent to subset_by_val(key, eq=True).
        '''
        return self.subset_by_val(key, eq=True)
    
    def subset_if_not(self, key: str) -> Self:
        '''
        Equivalent to subset_by_val(key, eq=False).
        '''
        return self.subset_by_val(key, eq=False)
    
    def subset_by_p(self, key: str, p_lo=None, p_hi=None) -> Self:
        val = self[key]
        lo, hi = None, None
        if p_lo is not None:
            lo = np.quantile(val, p_lo)
        if p_hi is not None:
            hi = np.quantile(val, p_hi)
        return self.subset_by_val(key, lo=lo, hi=hi)    
    
    @classmethod
    def from_file(cls, path: Path | str, **init_kw):
        '''
        Create the sample from a file. The file should be in HDF5 format,
        with datasets containing the galaxy attributes.
        
        @path: path to the file.
        @init_kw: additional keyword arguments passed to __init__().
        '''
        data = h5.File.load_from(path)
        return cls(data, copy=False, **init_kw)
    
    def to_file(self, path: Path | str, f_flag='x', **dump_kw):
        '''
        @f_flag, dump_kw: passed to h5.File.dump_to().
        '''
        h5.File.dump_to(path, data=self.data, f_flag=f_flag, **dump_kw)