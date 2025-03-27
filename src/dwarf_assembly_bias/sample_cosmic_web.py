from __future__ import annotations
import typing
from typing import Tuple, Any, Self
from pyhipp.core import abc, DataDict, Num, DataTable
from pyhipp.field.cubic_box import TidalClassifier
from .sample import SubhaloSet

# not used
class CosmicWebSampling(abc.HasDictRepr):

    repr_attr_keys = ['parent_samp', 'tid_clsf']

    def __init__(self, parent_samp: SubhaloSet,
                 tid_clsf: TidalClassifier) -> None:

        self.parent_samp = parent_samp
        self.tid_clsf = tid_clsf

    def sample_of_web_type(self, web_type: str):
        parent_samp, tid_clsf = self.parent_samp, self.tid_clsf

        x = parent_samp['x']
        sel = tid_clsf.web_type_at(x) == tid_clsf.web_types[web_type]

        return parent_samp.subset(sel)
