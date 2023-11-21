import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Kpoints

import logging
logger = logging.getLogger(__name__)

def generate_cell_parameters(qe, params):
    if qe.struct.use_ibrav:
        data = None
    else:
        data = qe.struct.structure.lattice.matrix
    return {
        'key': 'CELL_PARAMETERS',
        'option': 'angstrom',
        'data': data,
    }

def generate_atomic_species(qe, params):
    data = []

    if qe.pp_list is not None:
        for aname, ename in zip(qe.struct.atom_names, qe.struct.elem_names):
            pp = qe.pp_list.at[ename, "pseudopotential"]
            ppfile = "{}.{}{}.UPF".format(ename, ("rel-" if qe.is_soc else ""), pp)
            data += [[aname, Element(ename).atomic_mass, ppfile]]
    else:
        logger.error("generate_atomic_species: pp_file is not specified")
        raise RuntimeError("generate_atomic_species: pp_file not specified")
    return {
        'key': 'ATOMIC_SPECIES',
        'option': None,
        'data': data
    }

def generate_atomic_positions(qe, params):
    data = []
    option = params.get("option", "crystal")

    ign = None
    if "ignore_species" in params:
        ign = params["ignore_species"]
        if np.isscalar(ign):
            ign = [ign]

    for atom, coords in zip(qe.struct.atoms, qe.struct.structure.frac_coords):
        idx = qe.struct.atom_types.index(atom)
        atom_name = qe.struct.atom_names[idx]
        if ign is not None:
            data += [[atom_name, *coords, *([0,0,0] if atom_name in ign else []) ]]
        else:
            data += [[atom_name, *coords]]

    return {
        'key': 'ATOMIC_POSITIONS',
        'option': option,
        'data': data
    }

def generate_k_points(qe, params):

    def _get_kmesh(kresolution):
        nn = 2 * np.pi / kresolution
        kmesh = [int(np.ceil(nn/lp)) for lp in qe.struct.structure.lattice.abc]
        return kmesh

    def _get_kmesh_by_vol(vol_density):
        return Kpoints.automatic_density_by_vol(qe.struct.structure, vol_density, force_gamma=True).kpts[0]

    data = []
    option = params.get("option", "automatic")

    if option == "gamma":
        pass
    elif option in ["automatic", "crystal"]:
        # kgrid: [n1,n2,n3]
        # kgrid: n or [n] -> [n,n,n]
        # kgrid: [n1,n2]  -> [n1,n2,1]
        if "grid" in params:
            kgrid = params["grid"]
            if isinstance(kgrid, (int, float, str)):
                kgrid = [ int(kgrid) ] * 3
            elif isinstance(kgrid, list) and len(kgrid) > 0:
                if len(kgrid) == 1:
                    kgrid = kgrid * 3
                elif len(kgrid) == 2:
                    kgrid = kgrid + [1]
                elif len(kgrid) == 3:
                    pass
                else:
                    kgrid = kgrid[0:3]
            else:
                raise ValueError("generage_k_points: invalid grid {}".format(kgrid))
        else:
            vol_density = params.get("vol_density", 0)
            if vol_density > 0:
                kgrid = _get_kmesh_by_vol(vol_density)
            else:
                kres = params.get("k_resolution", 0.15)
                kgrid = _get_kmesh(kres)

        if option == "automatic":
            kshift = params.get("kshifts", [0, 0, 0])
            data = [ kgrid + kshift ]
        else:
            ngrid = kgrid[0] * kgrid[1] * kgrid[2]

            def _kval(kpt):
                return kpt - np.round(kpt + 1.e-10)

            data = [ [ ngrid ] ]
            data += [ [_kval(kx/kgrid[0]), _kval(ky/kgrid[1]), _kval(kz/kgrid[2]), 1.0/ngrid]
                      for kx in range(kgrid[0])
                      for ky in range(kgrid[1])
                      for kz in range(kgrid[2]) ]
    else:
        logger.error("generate_k_points: unknown option {}".format(option))
        return {}

    return {
        'key': 'K_POINTS',
        'option': option,
        'data': data
    }

