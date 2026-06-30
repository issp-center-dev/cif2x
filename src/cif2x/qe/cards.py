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
            ppfile = qe._pp_filename(ename)
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

    option = params.get("option", "automatic")
    mode = getattr(qe, "mode", None)
    band_options = ("crystal_b", "tpiba_b")

    # A bands calculation must use a band-path option; reject a silent uniform
    # mesh (wrong physics) instead of generating one.
    if mode == "bands" and option not in band_options:
        raise ValueError(
            "mode 'bands' requires a band-path K_POINTS option such as "
            "'crystal_b' (got option '{}'). Set 'option: crystal_b' in the "
            "K_POINTS content block.".format(option)
        )

    data = []

    if option in band_options:
        # honour a hand-written block verbatim (manual escape hatch)
        if params.get("data"):
            return {'key': 'K_POINTS', 'option': option, 'data': params["data"]}
        if option == "tpiba_b":
            raise ValueError(
                "K_POINTS option 'tpiba_b' is not generated automatically; "
                "supply explicit data or use 'crystal_b'."
            )
        data = _generate_band_path(qe, params)
    elif option == "gamma":
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
        # unknown option: keep the card as-is rather than dropping it, so a
        # hand-written K_POINTS block in a template survives.
        logger.warning("generate_k_points: passing through unrecognized option '{}'".format(option))
        return {'key': 'K_POINTS', 'option': option, 'data': params.get("data", [])}

    return {
        'key': 'K_POINTS',
        'option': option,
        'data': data
    }


def _generate_band_path(qe, params):
    """Build K_POINTS crystal_b rows from pymatgen's high-symmetry path.

    Each high-symmetry point becomes a row ``[kx, ky, kz, nk]`` in fractional
    coordinates; ``nk`` (the QE per-line integer = points to the next k-point)
    is ``line_npoints`` for interior points and 0 at the terminal point of each
    disjoint segment so bands.x does not interpolate across a break.
    """
    from pymatgen.symmetry.bandstructure import HighSymmKpath

    line_npoints = int(params.get("line_npoints", 20))
    kpath = HighSymmKpath(qe.struct.structure).kpath
    label_coords = kpath["kpoints"]
    segments = params.get("path") or kpath["path"]

    if any(isinstance(seg, str) for seg in segments):
        raise ValueError(
            "K_POINTS 'path' must be a list of label sequences "
            "(e.g. [['\\Gamma', 'X', 'M']]), not a flat list of labels."
        )

    rows = []
    for segment in segments:
        for i, label in enumerate(segment):
            if label not in label_coords:
                raise ValueError(
                    "unknown high-symmetry label '{}'. Available: {}.".format(
                        label, ", ".join(sorted(label_coords))
                    )
                )
            coord = label_coords[label]
            nk = 0 if i == len(segment) - 1 else line_npoints
            rows.append([float(coord[0]), float(coord[1]), float(coord[2]), nk])
    return [[len(rows)]] + rows

