# clpipe/file_types.py
"""
CLPipe-specific file types, subclassing base types from ceci.
SACCFile and FiducialCosmology adapted from LSSTDESC/TXPipe (file_types.py).
CosmosisFile and PythonFile are CLPipe additions, not present in core ceci.
"""
from ceci.file_types import DataFile, YamlFile
import yaml


class SACCFile(DataFile):
    """
    Adapted from TXPipe.
    """
    suffix = "sacc"

    @classmethod
    def open(cls, path, mode, **kwargs):
        import sacc
        if mode == "w":
            raise ValueError(
                "Do not use the open_output method to write sacc files. Use sacc.write_fits"
            )
        return sacc.Sacc.load_fits(path)

    def read_provenance(self):
        meta = self.file.metadata
        return {
            "uuid": meta.get("provenance/uuid", "UNKNOWN"),
            "creation": meta.get("provenance/creation", "UNKNOWN"),
            "domain": meta.get("provenance/domain", "UNKNOWN"),
            "username": meta.get("provenance/username", "UNKNOWN"),
        }

    def close(self):
        pass


class CosmosisFile(DataFile):
    """A CosmoSIS .ini file. Not present in core ceci."""
    suffix = "ini"


class PythonFile(DataFile):
    """A plain-text Python source file. Not present in core ceci."""
    suffix = "py"


class FiducialCosmology(YamlFile):
    """
    Adapted from TXPipe (LSSTDESC/TXPipe, file_types.py).
    TODO replace when CCL has more complete serialization tools.
    """
    def to_ccl(self, **kwargs):
        import pyccl as ccl

        with open(self.path, "r") as fp:
            params = yaml.load(fp, Loader=yaml.Loader)

        inits = dict(
            Omega_c=params["Omega_c"], Omega_b=params["Omega_b"], h=params["h"],
            n_s=params["n_s"],
            sigma8=None if params["sigma8"] == "nan" else params["sigma8"],
            A_s=None if params["A_s"] == "nan" else params["A_s"],
            Omega_k=params["Omega_k"], Neff=params["Neff"],
            w0=params["w0"], wa=params["wa"],
        )
        if ccl.__version__[0] == "2":
            inits.update(dict(
                bcm_log10Mc=params["bcm_log10Mc"], bcm_etab=params["bcm_etab"],
                bcm_ks=params["bcm_ks"], mu_0=params["mu_0"], sigma_0=params["sigma_0"],
            ))
        if "z_mg" in params:
            inits["z_mg"] = params["z_mg"]
            inits["df_mg"] = params["df_mg"]
        if "m_nu" in params:
            inits["m_nu"] = params["m_nu"]
            inits["m_nu_type"] = "list"

        inits.update(kwargs)
        return ccl.Cosmology(**inits)
