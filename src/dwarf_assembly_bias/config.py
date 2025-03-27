from __future__ import annotations
import os
from pyhipp_sims import sims
from pathlib import Path
from pyhipp import plot
from pyhipp.core import DataDict
from pyhipp.io import h5

plot.runtime_config.use_stylesheets('mathtext-it')


class ColorSets:
    dark2 = plot.ColorSeq.predefined('dark2').get_rgba()
    tab10 = plot.ColorSeq.predefined('dark2').get_rgba()
    set1 = plot.ColorSeq.predefined('set1').get_rgba()
    set2 = plot.ColorSeq.predefined('set2').get_rgba()


class ProjPaths:
    proj_dir = Path(os.environ['MAHGIC_WORK_DIR']
                    ).resolve() / 'workspaces/obs/dwarf_bias/data'
    sims_dir = proj_dir / 'sims'
    obs_dir = proj_dir / 'obs'
    models_dir = proj_dir / 'models'
    figs_dir = proj_dir / 'figures'

    @staticmethod
    def sim_dir_of(sim_info: sims.SimInfo) -> Path:
        return ProjPaths.sims_dir / sim_info.name
    
    @staticmethod
    def man_sim_dir(sim_info: sims.SimInfo) -> SimDir:
        return SimDir(sim_info)

    @staticmethod
    def save_fig(file_name: str, **savefig_kw):
        plot.savefig(ProjPaths.figs_dir / file_name, **savefig_kw)


class SimDir:
    def __init__(self, sim_info: sims.SimInfo):

        base_dir = ProjPaths.sim_dir_of(sim_info)

        self.sim_info = sim_info
        self.base_dir = base_dir
        self.sample_dir = base_dir / 'samples'

    def dump_sample(self, obj, file_name: str, f_flag='ac', dump_flag='ac'):
        h5.File.dump_to(self.sample_dir/file_name, obj,
                        f_flag=f_flag, dump_flag=dump_flag)

    def load_sample(self, file_name: str):
        return h5.File.load_from(self.sample_dir/file_name)
